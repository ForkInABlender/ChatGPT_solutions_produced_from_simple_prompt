# Dylan Kenneth Eliot & GPT-4-Plugins (Beta Edition)

"""
This is mostly for parts to use in modeling an actual AI. But with the bits for emulating the brain.







https://chat.openai.com/share/2a9f0d74-5687-4eb8-a4ed-8195b8270516

"""
from pybrain3.structure import SoftmaxLayer, LinearLayer
from pybrain3.datasets import SupervisedDataSet
from pybrain3.structure import FeedForwardNetwork, FullConnection
from pybrain3.supervised.trainers import BackpropTrainer
from pybrain3.structure.modules.neuronlayer import NeuronLayer
from brian2 import *
from sympy import symbols, diff, lambdify
from rdkit import Chem
from rdkit.Chem import Descriptors
import numpy as np


class FeedForwardLayer(NeuronLayer):
		def __init__(self, indim, outdim):
				super().__init__(indim, outdim)
				self.weights = np.random.randn(indim, outdim)
				self.bias = np.random.randn(outdim)
		def _forwardImplementation(self, inbuf, outbuf):
				outbuf[:] = inbuf @ self.weights + self.bias
		def _backwardImplementation(self, outerr, inerr, outbuf, inbuf):
				self.weights -= inbuf[:, np.newaxis] @ outerr[np.newaxis, :]
				self.bias -= outerr
				inerr[:] = outerr @ self.weights.T

class EmbeddingLayer(LinearLayer):
		def __init__(self, vocab_size, embedding_dim):
				super().__init__(embedding_dim)
				self.embeddings = np.random.randn(vocab_size, embedding_dim)
		def _forwardImplementation(self, inbuf, outbuf):
				self.token_idx = np.argmax(inbuf)
				if self.token_idx <= len(self.embeddings):
						outbuf[:] = self.embeddings[self.token_idx]
				else:
						self.token_idx = len(self.embeddings)-1
						outbuf[:] = self.embeddings[self.token_idx]
		def _backwardImplementation(self, outerr, inerr, outbuf, inbuf):
				gradient = np.zeros_like(self.embeddings)
				gradient[self.token_idx] = outerr
				self.embeddings -= gradient
				inerr[:] = self.embeddings.T @ outerr

class GeLULayer(NeuronLayer):
		def __init__(self, dim):
				super().__init__(dim, dim)
		def _forwardImplementation(self, inbuf, outbuf):
				outbuf[:] = 0.5 * inbuf * (
						1 + np.tanh(np.sqrt(2 / np.pi) * (inbuf + 0.044715 * inbuf**3)))
		def _backwardImplementation(self, outerr, inerr, outbuf, inbuf):
				cdf = 0.5 * (1 + np.tanh(np.sqrt(2 / np.pi) * (inbuf + 0.044715 * inbuf**3)))
				pd = np.sqrt(2 / np.pi) * (1 + 0.134145 * inbuf**2) * (1 / np.cosh(np.sqrt(2 / np.pi) * (inbuf + 0.044715 * inbuf**3)))**2
				inerr[:] = outerr * (cdf + inbuf * pd)

class AttentionLayer(NeuronLayer):
		def __init__(self, indim, outdim):
				super().__init__(indim, outdim)
				self.attention_weights = np.random.rand(indim, outdim)
		def _forwardImplementation(self, inbuf, outbuf):
				outbuf[:] = np.dot(inbuf, self.attention_weights)
		def _backwardImplementation(self, outerr, inerr, outbuf, inbuf):
				gradient = np.dot(inbuf.T, outerr)
				self.attention_weights -= gradient
				inerr[:] = np.dot(outerr, self.attention_weights.T)

def softmax(x, axis=-1):
		"""Compute softmax values for each set of scores in x."""
		e_x = np.exp(x - np.max(x, axis=axis, keepdims=True))
		return e_x / e_x.sum(axis=axis, keepdims=True)

class MultiHeadSelfAttention(NeuronLayer):
		def __init__(self, indim, outdim, num_heads):
				super().__init__(indim, outdim)
				self.num_heads = num_heads
				self.depth = indim // num_heads
				self.W_q = np.random.rand(indim, indim)
				self.W_k = np.random.rand(indim, indim)
				self.W_v = np.random.rand(indim, indim)
				self.W_o = np.random.rand(indim, outdim)

		def scaled_dot_product_attention(self, Q, K, V):
				matmul_qk = np.dot(Q, K.T)
				d_k = Q.shape[-1]
				scaled_attention_logits = matmul_qk / np.sqrt(d_k)
				attention_weights = softmax(scaled_attention_logits, axis=-1)
				return np.dot(attention_weights, V)

		def _forwardImplementation(self, inbuf, outbuf):
				if len(inbuf.shape) == 1:
						inbuf = inbuf[np.newaxis, :]
				Q = np.dot(inbuf, self.W_q)
				K = np.dot(inbuf, self.W_k)
				V = np.dot(inbuf, self.W_v)
				Q = np.split(Q, self.num_heads, axis=1)
				K = np.split(K, self.num_heads, axis=1)
				V = np.split(V, self.num_heads, axis=1)

				attention_heads = []
				for i in range(self.num_heads):
						attention_head = self.scaled_dot_product_attention(Q[i], K[i], V[i])
						attention_heads.append(attention_head)

				# Concatenate attention heads and pass through the final linear layer
				concatenated_heads = np.concatenate(attention_heads, axis=1)
				outbuf[:] = np.dot(concatenated_heads, self.W_o)

		def _backwardImplementation(self, outerr, inerr, outbuf, inbuf):
				d_concat_heads = np.dot(outerr, self.W_o.T)
				d_attention_heads = np.split(d_concat_heads, self.num_heads, axis=1)
				dQ_total, dK_total, dV_total = 0, 0, 0
				for i in range(self.num_heads):
						d_out = d_attention_heads[i]
						dV = np.dot(self.attention_weights[i].T, d_out)
						d_attention_weights = np.dot(d_out, self.V[i].T)
						dQK = d_attention_weights * (1 - self.attention_weights[i]) * self.attention_weights[i]
						dQ = np.dot(dQK, self.K[i])
						dK = np.dot(dQK, self.Q[i])
						dQ_total += dQ
						dK_total += dK
						dV_total += dV
				inerr[:] = np.dot(dQ_total, self.W_q.T) + np.dot(dK_total, self.W_k.T) + np.dot(dV_total, self.W_v.T)
				self.W_q -= np.dot(inbuf.T, dQ_total)
				self.W_k -= np.dot(inbuf.T, dK_total)
				self.W_v -= np.dot(inbuf.T, dV_total)
				self.W_o -= np.dot(self.concatenated_heads.T, outerr)

class LayerNorm(NeuronLayer):
		def __init__(self, size, eps=1e-6):
				super().__init__(size, size)
				self.gamma = np.ones(size)
				self.beta = np.zeros(size)
				self.eps = eps
		def _forwardImplementation(self, inbuf, outbuf):
				mean = np.mean(inbuf)
				std = np.std(inbuf)
				outbuf[:] = self.gamma * (inbuf - mean) / (std + self.eps) + self.beta
		def _backwardImplementation(self, outerr, inerr, outbuf, inbuf):
				N = inbuf.size
				dbeta = np.sum(outerr)
				dgamma = np.sum((inbuf - np.mean(inbuf)) / (np.std(inbuf) + self.eps) * outerr)
				dinbuf = self.gamma / (N * np.std(inbuf) + self.eps) * (N * outerr - np.sum(outerr) - (inbuf - np.mean(inbuf)) / (np.std(inbuf) + self.eps) * np.sum(outerr * (inbuf - np.mean(inbuf))))
				inerr[:] = dinbuf

VOCAB_SIZE = 50257 # GPT-3 & 4 use that value --> 50257; otherwise, use 120 to interface with this model until retooled.
D_MODEL = 128
NUM_BLOCKS = 96 # 94 to 96 for GPT 3 & 4; otherwise set to 4.
NUM_HEADS = 128
FFN_DIM = 128
net = FeedForwardNetwork()
inLayer = LinearLayer(VOCAB_SIZE)
net.addInputModule(inLayer)
embedding = EmbeddingLayer(VOCAB_SIZE, D_MODEL)
net.addModule(embedding)
net.addConnection(FullConnection(inLayer, embedding))
attention = MultiHeadSelfAttention(D_MODEL, D_MODEL, NUM_HEADS)
net.addModule(attention)
net.addConnection(FullConnection(embedding, attention))
prev_layer = attention
for _ in range(NUM_BLOCKS):
		norm1 = LayerNorm(D_MODEL)
		net.addModule(norm1)
		net.addConnection(FullConnection(prev_layer, norm1))
		ffn1 = LinearLayer(D_MODEL, FFN_DIM)
		net.addModule(ffn1)
		net.addConnection(FullConnection(norm1, ffn1))
		gelu = GeLULayer(FFN_DIM)
		net.addModule(gelu)
		net.addConnection(FullConnection(ffn1, gelu))
		ffn2 = LinearLayer(FFN_DIM, D_MODEL)
		net.addModule(ffn2)
		net.addConnection(FullConnection(gelu, ffn2))
		norm2 = LayerNorm(D_MODEL)
		net.addModule(norm2)
		net.addConnection(FullConnection(ffn2, norm2))    
		prev_layer = norm2
outLayer = SoftmaxLayer(VOCAB_SIZE)
net.addOutputModule(outLayer)
net.addConnection(FullConnection(prev_layer, outLayer))

#from pybrain3.tools.xml.networkwriter import NetworkWriter
#from pybrain3.tools.xml.networkreader import NetworkReader
#NetworkWriter.writeToFile(net, 'my_model.xml')
#loaded_net = NetworkReader.readFrom('my_model.xml')
