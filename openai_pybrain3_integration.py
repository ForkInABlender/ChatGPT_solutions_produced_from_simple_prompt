# Dylan Kenneth Eliot & GPT-4-plugins (Alpha Edition)

"""


Now GPT could make use of this to come closer to being more than a token parser. :)

This should give you an idea of how to adapt the ELIXA project imitation that is GPT-3 & GPT-4.

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
