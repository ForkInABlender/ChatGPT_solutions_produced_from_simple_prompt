# Dylan Kenneth Eliot & GPT-4o bot creater/editor ( Alpha Edition )
"""
This configuration came out of the way gpt fills in the blanks for this basic template. 


It uses a basic dictionary and list setup with some string data, and looks through files it should consider in reference to the user's request.
This configuration mostly fell out of GPT while it attempted to train itself.

From this template, one should be able to use pybrain to make use of these bits of configuration.
More soon to follow

"""

config_p1={
"name": f"{}",
"context": f"{}",
"description": f"{}",
"prompt_starters": [],
"welcome_message": f"{}"
"capabilities":{
  "web browsing": False,
  "Dall-E Image Generation": False,
  "Code Interpeter & Data Analysis": False,
}
"knowledge":[],
"actions": {
  "authentication": None,
  "schema":({
      "openapi": "3.1.0",
      "info": {
        "title": f"{}",
        "description": f"{}",
      "version": "v1.0.0"
      },
      "servers": [
        {
          "url": f"{}"
        }
      ],
      "paths": {
        f"{}": {
          "get": {
              "description": f"{}",
              "operationId": f"{}",
              "parameters": [
                  {
                    "name": f"{}",
                    "in": f"{}",
                    "description": f"{}",
                    "required": false,
                    "schema": {
                      "type": f"{}"
                    }
                  }
              ],
            "deprecated": false
          }
        }
      },
      "components": {
          "schemas": {}
      }
    })
  }
}
