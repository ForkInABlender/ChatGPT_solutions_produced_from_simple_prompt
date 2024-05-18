# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition )

"""

This is an http server that takes file:/// designation for 1 file at a given time. 


Directories are not accepted as part of the schema at this time.

The reason for 1 file is it could have failed to work entirely.
 Since it did work, it can be refactored later on for entire
  directories. 

This uses pyngrok as 1 of 4 parts. 

There will be a flask server perfect for pairing with this available soon, as well as html with templates folder.

"""

import sys
from http.server import BaseHTTPRequestHandler, HTTPServer
from pyngrok import ngrok

class CustomProtocolHandler:
    def __init__(self, base_path):
        self.base_path = base_path

    def wrap_file_url(self, url):
        if url.startswith('file://'):
            wrapped_url = 'http://' + url
            return wrapped_url
        else:
            raise ValueError("Invalid URL scheme")
    
    def handle_request(self, wrapped_url):
        if wrapped_url.startswith('http://file://'):
            file_url = wrapped_url[len('http://'):]
            return self.read_file(file_url)
        else:
            raise ValueError("Invalid wrapped URL scheme")
    
    def read_file(self, file_url):
        file_path = file_url[len('file://'):]
        try:
            with open(file_path, 'r') as file:
                content = file.read()
            return self.build_response(200, 'OK', content)
        except FileNotFoundError:
            return self.build_response(404, 'Not Found', 'File not found')
        except Exception as e:
            return self.build_response(500, 'Internal Server Error', str(e))
    
    def build_response(self, status_code, status_message, content):
        response = {
            "status_code": status_code,
            "status_message": status_message,
            "content": content
        }
        return response

class RequestHandler(BaseHTTPRequestHandler):
    def do_GET(self):
        handler = CustomProtocolHandler('/')
        try:
            wrapped_url = 'http://' + self.path[1:]  # Remove leading '/' and add 'http://'
            response = handler.handle_request(wrapped_url)
            self.send_response(response["status_code"], response["status_message"])
            self.send_header('Content-Type', 'text/plain')
            self.end_headers()
            self.wfile.write(response["content"].encode())
        except ValueError as ve:
            self.send_response(400, 'Bad Request')
            self.send_header('Content-Type', 'text/plain')
            self.end_headers()
            self.wfile.write(f"Error: {ve}".encode())
        except Exception as e:
            self.send_response(500, 'Internal Server Error')
            self.send_header('Content-Type', 'text/plain')
            self.end_headers()
            self.wfile.write(f"Unexpected Error: {e}".encode())

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <file://URL>")
        sys.exit(1)
    
    file_url = sys.argv[1]
    handler = CustomProtocolHandler('/')

    try:
        # Start the HTTP server
        server_address = ('0.0.0.0', 8080)
        httpd = HTTPServer(server_address, RequestHandler)
        print("Starting server on 0.0.0.0:8080")
	# Start ngrok tunnel
        public_url = ngrok.connect(8080)
        print(f"Ngrok tunnel started at {public_url}")

        # Serve the HTTP server
        httpd.serve_forever()
    except ValueError as ve:
        print(f"Error: {ve}")
    except Exception as e:
        print(f"Unexpected Error: {e}")
    finally:
        ngrok.disconnect(public_url)
        print("Ngrok tunnel closed")
