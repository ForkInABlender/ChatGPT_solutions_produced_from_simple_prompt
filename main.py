# Dylan Kenneth Eliot & GPT-4-plugins (Alpha Edition)

"""

Why have I laid gpt-neo aka GPT3 out on the gpt-2 tokenizers? 

Well, this, with what you see below is the basic way GPT-3 & GPT-4 work at the heart of them. Now this, with the pybrain3 neural networks and 
 a little dev ops polish, It can work like the up to date models using custom plugins. Do note that this is a basic layout. But, it is a good
	starting point.

* runs the functions
* responds and continues conversing like normal but with no "comprehension". needs attention heads added separately
* small enough that it fits on every UTD android device.
* can be tailored for needs by design
* good point to template from for an offline version of GPT3 as GPT-NEO 


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
