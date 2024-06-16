# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition )

"""
This allows for imitation of dockerd.
 As long as types are properly marshelled and unmarshelled, seemless interaction with docker as if it is dockerd itself...

For more information, consult the documentation: https://docs.docker.com/engine/api/v1.45/

Product class: 
     rating: kepler/nebulus; otherwise considered safe for development
     potential: limit unknown

The purpose was to explore an option that was more conformant the the logic of docker and not so much as be stuck to architectural limitations.
 But instead have a clear path that allowed for running even on a cellphone running python inside android. With that being said, dockerd's remaining logic
  could be implemented in python with golang & ctypes interfacing now that I've figured out the hard part.

The next steps is it needs unicorn-engine and kernel emulation. Plus more documentation surfing. Imagine if you will the bash operator from hell had a 
 protege. Only it wears the personality of the character "tony stark". "that guy's playing galiga" 

"""


from flask import Flask, jsonify, request

app = Flask(__name__)

@app.route('/_ping', methods=['GET'])
def ping():
    return '', 200

@app.route('/v1.24/containers/json', methods=['GET'])
def list_containers():
    # Fake container data
    container_info = [{
        "Id": "1234567890abcdef",
        "Names": ["/fake_container"],
        "Image": "fake_image:latest",
        "ImageID": "sha256:1234567890abcdef",
        "Command": "/bin/sh -c 'while true; do echo hello world; sleep 1; done'",
        "Created": 1625158000,
        "Ports": [],
        "Labels": {},
        "State": "running",
        "Status": "Up 5 minutes",
        "HostConfig": {
            "NetworkMode": "default"
        },
        "NetworkSettings": {
            "Networks": {
                "bridge": {
                    "IPAMConfig": None,
                    "Links": None,
                    "Aliases": None,
                    "NetworkID": "1234567890abcdef",
                    "EndpointID": "1234567890abcdef",
                    "Gateway": "172.17.0.1",
                    "IPAddress": "172.17.0.2",
                    "IPPrefixLen": 16,
                    "IPv6Gateway": "",
                    "GlobalIPv6Address": "",
                    "GlobalIPv6PrefixLen": 0,
                    "MacAddress": "02:42:ac:11:00:02",
                    "DriverOpts": None
                }
            }
        },
        "Mounts": []
    }]
    return jsonify(container_info), 200

@app.errorhandler(404)
def page_not_found(e):
    return '', 404

if __name__ == '__main__':
    app.run(host='127.0.0.1', port=2375)
