# Dylan Kenneth Eliot & GPT-4-Plugins (Beta Edition)

"""
This is mostly for parts to use in modeling an actual AI. But with the bits for emulating the brain.

05/05/2026 @ 6:10 to 6:33 pm update --- this file has been updated so it also runs on android cooies of the pybrain3,
 scipy, and numpy libraries. due keep in mind that this is designed to handle 50'257 input and output of data





https://chat.openai.com/share/2a9f0d74-5687-4eb8-a4ed-8195b8270516

"""


from pybrain3.structure import SoftmaxLayer, LinearLayer
from pybrain3.datasets import SupervisedDataSet
from pybrain3.structure import FeedForwardNetwork, FullConnection
from pybrain3.supervised.trainers import BackpropTrainer
from pybrain3.structure.modules.neuronlayer import NeuronLayer
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
				self.embeddings[self.token_idx] -= outerr
				inerr[:] = 0

class GeLULayer(NeuronLayer):
		def __init__(self, dim):
				super().__init__(dim, dim)
		def _forwardImplementation(self, inbuf, outbuf):
				outbuf[:] = 0.5 * inbuf * (
						1 + np.tanh(np.sqrt(2 / np.pi) * (inbuf + 0.044715 * inbuf**3)))
		def _backwardImplementation(self, outerr, inerr, outbuf, inbuf):
				x = np.clip(inbuf, -10, 10)
				cdf = 0.5 * (1 + np.tanh(np.sqrt(2 / np.pi) * (x + 0.044715 * x**3)))
				pd = np.sqrt(2 / np.pi) * (1 + 0.134145 * x**2) * (1 / np.cosh(np.sqrt(2 / np.pi) * (x + 0.044715 * x**3)))**2
				inerr[:] = outerr * (cdf + x * pd)

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
				self._inbuf = inbuf
				Q = np.dot(inbuf, self.W_q)
				K = np.dot(inbuf, self.W_k)
				V = np.dot(inbuf, self.W_v)
				self.Q = np.split(Q, self.num_heads, axis=1)
				self.K = np.split(K, self.num_heads, axis=1)
				self.V = np.split(V, self.num_heads, axis=1)

				attention_heads = []
				self.attention_weights = []
				for i in range(self.num_heads):
						matmul_qk = np.dot(self.Q[i], self.K[i].T)
						aw = softmax(matmul_qk / np.sqrt(self.depth), axis=-1)
						self.attention_weights.append(aw)
						attention_heads.append(np.dot(aw, self.V[i]))

				self.concatenated_heads = np.concatenate(attention_heads, axis=1)
				outbuf[:] = np.dot(self.concatenated_heads, self.W_o)

		def _backwardImplementation(self, outerr, inerr, outbuf, inbuf):
				seq_len = self._inbuf.shape[0]
				outerr_2d = outerr.reshape(seq_len, -1)
				d_concat_heads = np.dot(outerr_2d, self.W_o.T)
				d_attention_heads = np.split(d_concat_heads, self.num_heads, axis=1)
				dQ_heads, dK_heads, dV_heads = [], [], []
				for i in range(self.num_heads):
						d_out = d_attention_heads[i]
						dV_heads.append(np.dot(self.attention_weights[i].T, d_out))
						d_aw = np.dot(d_out, self.V[i].T)
						dQK = d_aw * self.attention_weights[i] * (1 - self.attention_weights[i])
						dQ_heads.append(np.dot(dQK, self.K[i]))
						dK_heads.append(np.dot(dQK, self.Q[i]))
				dQ_total = np.concatenate(dQ_heads, axis=1)
				dK_total = np.concatenate(dK_heads, axis=1)
				dV_total = np.concatenate(dV_heads, axis=1)
				dinput = np.dot(dQ_total, self.W_q.T) + np.dot(dK_total, self.W_k.T) + np.dot(dV_total, self.W_v.T)
				inerr[:] = dinput.reshape(-1)[:inerr.size]
				self.W_q -= np.dot(self._inbuf.T, dQ_total)
				self.W_k -= np.dot(self._inbuf.T, dK_total)
				self.W_v -= np.dot(self._inbuf.T, dV_total)
				self.W_o -= np.dot(self.concatenated_heads.T, outerr_2d)

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
D_MODEL = 64
NUM_BLOCKS = 48 # 94 to 96 for GPT 3 & 4; otherwise set to 4.
NUM_HEADS = 16
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
net.sortModules()


from pybrain3.tools.xml.networkwriter import NetworkWriter
from pybrain3.tools.xml.networkreader import NetworkReader
NetworkWriter.writeToFile(net, 'my_model.xml')
#loaded_net = NetworkReader.readFrom('my_model.xml')
