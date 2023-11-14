# Dylan Kenneth Eliot & GPT-4-plugins (Alpha Edition)

"""


Now GPT could make use of this to come closer to being more than a token parser. :)

This script integrates several libraries and functionalities, primarily focusing on chemical data processing, neural network training, and interaction 
   with OpenAI's API. Here's a breakdown of its components and functionalities:
   
1. **Import Libraries**: The script imports necessary libraries including `openai` for accessing OpenAI's API, `json` for data serialization, `numpy` for
   numerical operations, `rdkit` for chemical informatics, `pybrain3` for neural network operations, and `sympy` for symbolic mathematics.
   
2. **OpenAI API Key Setup**: It sets up the OpenAI API key, which is necessary for making requests to OpenAI's services.

3. **Chemical Data Processing with RDKit**: 
    - It uses RDKit to process a chemical compound represented by a SMILES (Simplified Molecular Input Line Entry System) string.
    - It calculates the molecular weight (`mol_weight`) and the logarithm of the partition coefficient (`logp`) of the molecule.

4. **Symbolic Calculations with SymPy**: 
    - It defines a symbolic expression using SymPy and simplifies it.
    - The expression is then evaluated using the molecular weight and logP values calculated earlier.

5. **Neural Network Setup with PyBrain3**:
    - It prepares a dataset for a neural network, using the features (molecular weight and logP) and the target (evaluated expression).
    - A neural network is built and trained using backpropagation. The network has 2 input neurons, 3 hidden neurons, and 1 output neuron.

6. **Neural Network Processing Function**:
    - It defines a function `process_nn` that takes input data, processes it through the neural network, and returns the prediction in JSON format.

7. **OpenAI Function Definition**:
    - It defines a function for OpenAI to call (`process_nn`), including its parameters and description.

8. **OpenAI API Interaction**:
    - The script makes an API call to OpenAI's ChatCompletion model, sending a user message and the defined function.
    - It then checks for a function call in the response from OpenAI. If there's a function call, it processes the input data through the neural network
         and appends the response to the message chain.
    - A second API call is made to OpenAI with the updated messages, including the neural network's output.

9. **Output**:
    - Finally, the script prints the content of the response from the second OpenAI API call.

This script is a sophisticated integration of chemical informatics, symbolic computation, neural network processing, and interaction with OpenAI's API, 
  demonstrating a complex workflow that spans multiple domains.

This should give you an idea of how to adapt the ELIZA project imitation that is GPT-3 & GPT-4 to your model.

Happy Templating!


 


"""



import openai
import json
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors
from pybrain3.datasets import SupervisedDataSet
from pybrain3.supervised.trainers import BackpropTrainer
from pybrain3.tools.shortcuts import buildNetwork
import sympy as sp

# Set up OpenAI API key
openai.api_key = 'your-api-key'

# RDKit: Compute molecular properties
smiles = "[Zn](OOC)(OOC)C[N+](C)C(C)C(C)[N+](C)C(O)=O.[Mg](OOC)(OOC)C[N+](C)C(C)C(C)[N+](C)C(O)=O.[SO4]"
molecule = Chem.MolFromSmiles(smiles)
mol_weight = Descriptors.MolWt(molecule)
logp = Descriptors.MolLogP(molecule)

# SymPy for symbolic calculations
x, y = sp.symbols('x y')
expression = sp.simplify(x**2 + y**2)
evaluated_expression = expression.subs({x: mol_weight, y: logp})

# Prepare dataset for PyBrain3 using NumPy
features = np.array([mol_weight, logp])
targets = np.array([evaluated_expression])
ds = SupervisedDataSet(2, 1)
ds.addSample(features, targets)

# Build and train the neural network
net = buildNetwork(2, 3, 1)
trainer = BackpropTrainer(net, ds)
trainer.trainUntilConvergence()

# Define a function to process input through the neural network
def process_nn(input_data):
    prediction = net.activate(input_data)
    return json.dumps({"prediction": prediction.tolist()})

# Define the function for OpenAI to call
functions = [{
    "name": "process_nn",
    "description": "Process input through the neural network",
    "parameters": {
        "type": "object",
        "properties": {
            "input_data": {
                "type": "array",
                "description": "Input data for the neural network",
            },
        },
        "required": ["input_data"],
    },
}]

# Example OpenAI API call
messages = [{"role": "user", "content": "Process this data through the neural network"}]
response = openai.ChatCompletion.create(
    model="gpt-3.5-turbo-0613",
    messages=messages,
    functions=functions,
    function_call="auto",
)

# Check if there's a function call in the response
response_message = response["choices"][0]["message"]
if response_message.get("function_call"):
    available_functions = {
        "process_nn": process_nn,
    }
    function_name = response_message["function_call"]["name"]
    function_to_call = available_functions[function_name]
    function_args = json.loads(response_message["function_call"]["arguments"])
    function_response = function_to_call(
        input_data=function_args.get("input_data"),
    )
    messages.append(response_message)
    messages.append({
        "role": "function",
        "name": function_name,
        "content": function_response,
    })
    second_response = openai.ChatCompletion.create(
        model="gpt-3.5-turbo-0613",
        messages=messages,
    )
    print(second_response["choices"][0]["message"]['content'])
