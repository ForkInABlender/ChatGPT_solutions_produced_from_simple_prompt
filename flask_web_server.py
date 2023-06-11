# Dylan Kenneth Eliot & GPT


from flask import Flask, render_template
from flask_socketio import SocketIO, emit
import torch
import numpy as np
from transformers import GPT2LMHeadModel, GPT2Tokenizer
import anvil.server

app = Flask(__name__)
app.config['SECRET_KEY'] = 'your-secret-key'
socketio = SocketIO(app)

def get_gpt2_model_and_tokenizer(model_name):
    model = GPT2LMHeadModel.from_pretrained(model_name)
    tokenizer = GPT2Tokenizer.from_pretrained(model_name)
    return model, tokenizer

def train_gpt2_on_offset(gpt2_model, input_ids, offset, num_iterations):
    optimizer = torch.optim.Adam(gpt2_model.parameters())
    loss_fn = torch.nn.CrossEntropyLoss()

    for _ in range(num_iterations):
        optimizer.zero_grad()
        input_ids_offset = input_ids[:, :offset + 1]
        output = gpt2_model(input_ids_offset)
        loss = loss_fn(output.logits[:, :-1].flatten(0, 1), input_ids_offset[:, 1:].flatten())
        loss.backward()
        optimizer.step()

def generate_text(gpt2_model, gpt2_tokenizer, input_text, max_length=100, num_return_sequences=1):
    input_ids = gpt2_tokenizer.encode(input_text, return_tensors="pt")
    with torch.no_grad():
        output = gpt2_model.generate(input_ids, max_length=max_length, num_return_sequences=num_return_sequences)
    output_text = gpt2_tokenizer.decode(output[0], skip_special_tokens=True)
    return output_text

@app.route('/')
def index():
    return render_template('index.html')

@socketio.on('generate_modified_text')
def handle_generate_modified_text(input_text):
    output_text_modified = input_text
    for offset in range(num_offsets):
        train_gpt2_on_offset(gpt2_model, input_ids, offset, num_iterations)
        input_text_offset = generate_text(gpt2_model, gpt2_tokenizer, output_text_modified)
        output_text_modified = generate_text(gpt2_model, gpt2_tokenizer, input_text_offset)
    emit('modified_text', output_text_modified)

if __name__ == '__main__':
    model_name = "gpt2-medium"
    num_iterations = 10
    num_offsets = 113333

    gpt2_model, gpt2_tokenizer = get_gpt2_model_and_tokenizer(model_name)
    input_ids = gpt2_tokenizer.encode("Hello", return_tensors="pt")

    # Anvil Uplink integration
    anvil.server.connect("<ANVIL_UPLINK_KEY>")
    @anvil.server.callable
    def generate_modified_text(input_text):
        output_text_modified = input_text
        for offset in range(num_offsets):
            train_gpt2_on_offset(gpt2_model, input_ids, offset, num_iterations)
            input_text_offset = generate_text(gpt2_model, gpt2_tokenizer, output_text_modified)
            output_text_modified = generate_text(gpt2_model, gpt2_tokenizer, input_text_offset)
        return output_text_modified

    socketio.run(app)

