# Dylan Kenneth Eliot

"""
This code does function passing from the server to the client; the python interpreter being the client in this example.

One could also use this to go the other way in bridging. 

"""


import subprocess
import time
import xmlrpc.client

# Step 1: Start the Ruby XML-RPC server as a background process
ruby_server_code = """
require 'xmlrpc/server'

server = XMLRPC::Server.new(8080)

class RubyService
  def say_hello(name)
    "Hello, #{name} from Ruby!"
  end

  def add_numbers(a, b)
    a + b
  end
end

server.add_handler("ruby_service", RubyService.new)

# Keep the server running
trap("INT") { exit }
server.serve
"""

# Write the Ruby code to a temporary file
with open("ruby_server.rb", "w") as file:
    file.write(ruby_server_code)

# Start the Ruby server as a subprocess
ruby_process = subprocess.Popen(["ruby", "ruby_server.rb"])

# Give the server a moment to start
time.sleep(1)

try:
    # Step 2: Connect to the Ruby XML-RPC server and call methods
    server = xmlrpc.client.ServerProxy("http://localhost:8080")

    # Call the say_hello method
    name = "Python"
    response = server.ruby_service.say_hello(name)
    print("Response from Ruby (say_hello):", response)

    # Call the add_numbers method
    a, b = 5, 7
    result = server.ruby_service.add_numbers(a, b)
    print("Response from Ruby (add_numbers):", result)

finally:
    # Clean up: Terminate the Ruby server process
    ruby_process.terminate()
    ruby_process.wait()
