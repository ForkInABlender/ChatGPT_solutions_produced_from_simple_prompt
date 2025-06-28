# Dylan Kenneth Eliot

"""
Create a flask app that gets its localtunnel "password" (ip assigned to the tunnel) then setup the tunnel given a string for a localtunnel subdomain...

"""



from flask import Flask, render_template, send_from_directory
import subprocess
import threading
import os
from time import sleep


app = Flask(__name__)

@app.route('/')
def index():
	return render_template('index.html')
#
def get_tunnel_auth():
	try:
		result = subprocess.check_output(['wget', '-q', '-O', '-', 'https://loca.lt/mytunnelpassword'])
		auth_token = result.decode().strip()
		print(f"[+] Got auth token: {auth_token}")
		return auth_token
	except subprocess.CalledProcessError as e:
		print("[-] Failed to get tunnel auth:", e)
		return None
#
def start_localtunnel(port, auth_token):
	print(f"[+] Starting LocalTunnel on port {port}")
	cmd = ['lt', '--subdomain',	'wakwinjnachal', '--port', str(port)]
	tunnel = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	for line in iter(tunnel.stdout.readline, b''):
		decoded = line.decode().strip()
		print("[localtunnel] " + decoded)
		if decoded.count("your url is:") == 1:
			break
	return tunnel
#
if __name__ == '__main__':
	threading.Thread(target=start_localtunnel, args=(8000, get_tunnel_auth())).start()
	#app.run(host="127.0.0.1", port=8000, debug=True)
	#sleep(5)


	def f(*z):
		x = dict(*z)
		app.run(x['host'], x['port'],debug=x['debug'])
		#app.run(z['host'], z['port'],debug=z['debug'])
	
	f({'host': '127.0.0.1', 'port': int(8000), 'debug': True})
