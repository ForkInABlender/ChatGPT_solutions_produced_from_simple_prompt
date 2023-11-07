# Dylan Kenneth Eliot & GPT-4-plugins (Alpha Edition)


"""
Why have I integrated GPT-Neo (also known as GPT-3) with GPT-2 tokenizers?

The reason is that at their core, both GPT-3 and GPT-4 operate on similar principles. By leveraging PyBrain3 neural networks and applying some development operations finesse, we can emulate the functionality of the latest models with the use of custom plugins. It's important to note that this setup is quite rudimentary. However, it serves as an excellent starting point.

Key features of this setup include:

Execution of functions
Capability to engage in dialogue, although initially without "comprehension" – this requires separate implementation of attention mechanisms
Compact enough to be deployed on any UTD (up-to-date) Android device
Customizable to specific requirements by design
An ideal template for developing an offline version of GPT-3, akin to GPT-Neo
To begin, you'll need to install the necessary packages using the following commands:
```
pip install torch==1.13.1+cpu torchvision=0.16.0+cpu -f https://download.pytorch.org/whl/cpu/torch_stable.html
pip install transformers
```
After running the script, it will automatically download the chosen model. Depending on the model, the size can range from over 530MB to as large as 10.7GB. For testing purposes, I opted for the smallest model to conserve time and computational resources. This approach also simplifies the integration of brain emulation techniques using libraries like NumPy and PyBrain.
"""


from transformers import GPTNeoForCausalLM, GPT2Tokenizer

# Load the model and tokenizer
model_name = "EleutherAI/gpt-neo-125M"  # or another size like "2.7B"
model = GPTNeoForCausalLM.from_pretrained(model_name)
tokenizer = GPT2Tokenizer.from_pretrained(model_name)

# Define some functions that can be called by the user
def add_numbers(a, b):
    return a + b

def subtract_numbers(a, b):
    return a - b

# Function to generate text based on a prompt
def generate_text(prompt):
    input_ids = tokenizer(prompt, return_tensors="pt").input_ids
    gen_tokens = model.generate(input_ids, do_sample=True, attention_mask=21, temperature=0.9, max_length=100)
    gen_text = tokenizer.batch_decode(gen_tokens)[0]
    return gen_text

# Function to call a specific function based on the user's request
def call_function_by_request(request):
    # A dictionary mapping keywords to functions
    functions = {
        'add': add_numbers,
        'subtract': subtract_numbers,
    }
    
    # Split the request to find the keyword and arguments
    parts = request.split()
    if parts[0] in functions:
        # Assume the rest of the parts are arguments
        args = [int(arg) for arg in parts[1:]]
        # Call the function with the arguments
        return functions[parts[0]](*args)
    else:
        # If the first word is not a keyword, generate text
        return generate_text(request)

# Example usage
user_request = "add 5 10"
result = call_function_by_request(user_request)
print(result)  # Output will be 15

# For generating text
user_request = "Once upon a time"
result = call_function_by_request(user_request)
print(result)  # Output will be generated text
