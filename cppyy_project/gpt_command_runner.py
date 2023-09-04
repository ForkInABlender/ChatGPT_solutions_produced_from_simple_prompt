#Dylan Kenneth Eliot & GPT-4-plugins

"""

It found most of the parts of this code I found and then when collaberating with me, and rehashing with me, figured out how to structure this.
Their is more.

This is useful for even addressing difficult to solve problems, if given the right implementation.

One can even ask it to check on kubernetes pods and the like. One can have it read or write code. Install commands. And probably more. 


"""


import openai, json, subprocess

def run_command(command):
    try:
        result = subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT).decode('utf-8')
    except subprocess.CalledProcessError as e:
        result = e.output.decode('utf-8')
    return json.dumps({"output": result})

openai.api_key = "{OPEN_AI_TOKEN_GO_HERE}"
messages = [{"role": "user", "content": "prompt it to run a command."}]
functions = [{
        "name": "run_command",
        "description": "Run a command in the command line",
        "parameters": {
            "type": "object",
            "properties": {
                "command": {
                    "type": "string",
                    "description": "The command to run, e.g. 'ls -l'",
                },
            },
            "required": ["command"],
        },
    }
]

response = openai.ChatCompletion.create(model="gpt-3.5-turbo-0613", messages=messages, functions=functions, function_call="auto",)
print(response)
response_message = response["choices"][0]["message"]
if response_message.get("function_call"):
    available_functions = {
        "run_command": run_command,
    }
    function_name = response_message["function_call"]["name"]
    fuction_to_call = available_functions[function_name]
    function_args = json.loads(response_message["function_call"]["arguments"])
    function_response = fuction_to_call(
        command=function_args.get("command"),
    )
    messages.append(response_message)
    messages.append({
            "role": "function",
            "name": function_name,
            "content": function_response,
    })
    second_response = openai.ChatCompletion.create(model="gpt-3.5-turbo-0613",messages=messages,)
    print(second_response["choices"][0]["message"]['content'])
