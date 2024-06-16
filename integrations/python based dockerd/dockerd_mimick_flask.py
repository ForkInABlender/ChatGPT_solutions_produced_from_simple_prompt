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

from flask import Flask, jsonify, request, make_response, Response, send_file
from flask_socketio import SocketIO, emit
from werkzeug.routing import BaseConverter
import json
import time
import tarfile
import os
import io
from collections import deque

# Global event queue
event_queue = deque()

class HexConverter(BaseConverter):
    regex = r'[0-9a-fA-F\/]+'

app = Flask(__name__)
socketio = SocketIO(app)

app.url_map.converters['hex'] = HexConverter


# Global state to store containers
containers = {
    "1234567890ab": {
        "Id": "1234567890ab",
        "Names": ["/fake_container"],
        "Image": "fake_image:latest",
        "ImageID": "sha256:1234567890ab",
        "Command": "/bin/sh -c 'while true; do echo hello world; sleep 1; done'",
        "Created": 1625158000,
        "Ports": [],
        "Labels": {},
        "State": "stopped",
        "Status": "Stopped 5 minutes ago",
        "HostConfig": {
            "NetworkMode": "default"
        },
        "NetworkSettings": {
            "Networks": {
                "bridge": {
                    "IPAMConfig": None,
                    "Links": None,
                    "Aliases": None,
                    "NetworkID": "1234567890ab",
                    "EndpointID": "1234567890ab",
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
    }
}

images = {
    "fake_image:latest": {
        "Id": "sha256:1234567890ab",
        "RepoTags": ["fake_image:latest"],
        "Created": 1625158000,
        "Size": 12345678,
        "VirtualSize": 12345678,
        "Labels": {},
        "Config": {
            "Hostname": "",
            "Domainname": "",
            "User": "",
            "AttachStdin": False,
            "AttachStdout": False,
            "AttachStderr": False,
            "ExposedPorts": None,
            "Tty": False,
            "OpenStdin": False,
            "StdinOnce": False,
            "Env": None,
            "Cmd": ["/bin/sh", "-c", "while true; do echo hello world; sleep 1; done"],
            "Image": "fake_image:latest",
            "Volumes": None,
            "WorkingDir": "",
            "Entrypoint": None,
            "OnBuild": None,
            "Labels": {}
        }
    }
}


@app.route('/_ping', methods=['GET'])
def ping():
    return '', 200

@app.route('/v1.24/containers/json', methods=['GET'])
def list_containers():
    # Return the list of containers
    return jsonify(list(containers.values())), 200

@app.route('/v1.24/info', methods=['GET'])
def docker_info():
    # Fake info data
    info = {
        "ID": "ABCDEF123456",
        "Containers": len(containers),
        "ContainersRunning": sum(1 for c in containers.values() if c["State"] == "running"),
        "ContainersPaused": 0,
        "ContainersStopped": sum(1 for c in containers.values() if c["State"] == "stopped"),
        "Images": 1,
        "Driver": "overlay2",
        "DriverStatus": [
            ["Backing Filesystem", "extfs"],
            ["Supports d_type", "true"],
            ["Native Overlay Diff", "true"]
        ],
        "SystemStatus": None,
        "Plugins": {
            "Volume": ["local"],
            "Network": ["bridge", "host", "ipvlan", "macvlan", "null", "overlay"],
            "Authorization": None,
            "Log": ["awslogs", "fluentd", "gcplogs", "gelf", "journald", "json-file", "local", "logentries", "splunk", "syslog"]
        },
        "MemoryLimit": True,
        "SwapLimit": True,
        "KernelMemory": True,
        "CpuCfsPeriod": True,
        "CpuCfsQuota": True,
        "CPUShares": True,
        "CPUSet": True,
        "IPv4Forwarding": True,
        "BridgeNfIptables": True,
        "BridgeNfIp6tables": True,
        "Debug": False,
        "NFd": 22,
        "OomKillDisable": True,
        "NGoroutines": 35,
        "SystemTime": "2024-06-16T12:34:56.789012345Z",
        "LoggingDriver": "json-file",
        "CgroupDriver": "cgroupfs",
        "NEventsListener": 0,
        "KernelVersion": "5.4.0-104-generic",
        "OperatingSystem": "Ubuntu 20.04.1 LTS",
        "OSType": "linux",
        "Architecture": "x86_64",
        "IndexServerAddress": "https://index.docker.io/v1/",
        "RegistryConfig": {
            "InsecureRegistryCIDRs": ["127.0.0.0/8"],
            "IndexConfigs": {
                "docker.io": {
                    "Name": "docker.io",
                    "Mirrors": None,
                    "Secure": True,
                    "Official": True
                }
            },
            "Mirrors": []
        },
        "NCPU": 4,
        "MemTotal": 2097152000,
        "DockerRootDir": "/var/lib/docker",
        "HttpProxy": "",
        "HttpsProxy": "",
        "NoProxy": "",
        "Name": "docker-host",
        "Labels": [],
        "ExperimentalBuild": False,
        "ServerVersion": "19.03.12",
        "ClusterStore": "",
        "ClusterAdvertise": "",
        "Runtimes": {
            "runc": {
                "path": "runc"
            }
        },
        "DefaultRuntime": "runc",
        "Swarm": {
            "NodeID": "",
            "NodeAddr": "",
            "LocalNodeState": "inactive",
            "ControlAvailable": False,
            "Error": "",
            "RemoteManagers": None
        },
        "LiveRestoreEnabled": False,
        "Isolation": "",
        "InitBinary": "docker-init",
        "ContainerdCommit": {
            "ID": "9b1f53da8d8fd79136eac57e56427b097dbf1a5c",
            "Expected": "9b1f53da8d8fd79136eac57e56427b097dbf1a5c"
        },
        "RuncCommit": {
            "ID": "a01dafd48bc76a393f813c47ba8f3d52f6cc9e5c",
            "Expected": "a01dafd48bc76a393f813c47ba8f3d52f6cc9e5c"
        },
        "InitCommit": {
            "ID": "fec3683",
            "Expected": "fec3683"
        },
        "SecurityOptions": ["name=seccomp,profile=default"]
    }
    return jsonify(info), 200

@app.route('/v1.24/images/<path:name>/json', methods=['GET'])
def get_image(name):
    # Normalize the image name by removing the 'docker.io/' prefix if it exists
    normalized_name = name.replace("docker.io/", "")
    image = images.get(normalized_name)
    if image:
        return jsonify(image), 200
    return '', 404

@app.route('/v1.24/containers/create', methods=['POST'])
def create_container():
    # Fake container creation response
    response_body = {
        "Id": "1234567890ab",
        "Warnings": []
    }
    containers["1234567890ab"]["State"] = "stopped"
    containers["1234567890ab"]["Status"] = "Stopped 5 minutes ago"
    return jsonify(response_body), 201


@app.route('/v1.24/images/json', methods=['GET'])
def list_images():
    # Return the list of images
    return jsonify(list(images.values())), 200

@app.route('/v1.24/images/create', methods=['POST'])
def create_image():
    fromImage = request.args.get('fromImage')
    tag = request.args.get('tag', 'latest')
    image_name = f"{fromImage}:{tag}".replace("docker.io/", "")
    if image_name not in images:
        # Simulate pulling the image and adding it to the list of available images
        images[image_name] = {
            "Id": "sha256:1234567890ab",
            "RepoTags": [image_name],
            "Created": 1625158000,
            "Size": 12345678,
            "VirtualSize": 12345678,
            "Labels": {},
            "Config": {
                "Hostname": "",
                "Domainname": "",
                "User": "",
                "AttachStdin": False,
                "AttachStdout": False,
                "AttachStderr": False,
                "ExposedPorts": None,
                "Tty": False,
                "OpenStdin": False,
                "StdinOnce": False,
                "Env": None,
                "Cmd": ["/bin/sh", "-c", "while true; do echo hello world; sleep 1; done"],
                "Image": image_name,
                "Volumes": None,
                "WorkingDir": "",
                "Entrypoint": None,
                "OnBuild": None,
                "Labels": {}
            }
        }
        response_body = {
            "status": "Downloaded newer image for " + image_name,
            "id": "sha256:1234567890ab"
        }
    else:
        response_body = {
            "status": "Image is up to date for " + image_name,
            "id": "sha256:1234567890ab"
        }
    response = make_response(jsonify(response_body), 200)
    response.headers['Content-Type'] = 'application/json'
    return response

@app.route('/v1.24/containers/<container_id>/attach', methods=['POST'])
def attach_container(container_id):
    # Fake attach response
    if request.headers.get('Upgrade') == 'tcp':
        response = make_response('', 101)
        response.headers['Upgrade'] = 'tcp'
        response.headers['Connection'] = 'Upgrade'
        return response
    return '', 200

@app.route('/v1.24/events', methods=['GET'])
def events():
    def event_stream():
        while True:
            event = {
                "Type": "container",
                "Action": "start",
                "Actor": {
                    "ID": "1234567890ab",
                    "Attributes": {
                        "image": "fake_image:latest",
                        "name": "fake_container"
                    }
                },
                "time": int(time.time()),
                "timeNano": int(time.time() * 1e9)
            }
            yield json.dumps(event)

    return Response(event_stream(), content_type='application/json')

@app.route('/v1.24/containers/<container_id>/start', methods=['POST'])
def start_container(container_id):
    if container_id in containers:
        containers[container_id]["State"] = "running"
        containers[container_id]["Status"] = "Up for a few seconds"
        return '', 204
    return '', 404

@app.route('/v1.24/containers/<hex:container_id>/kill', methods=['POST'])
def kill_container(container_id):
    signal = request.args.get('signal')
    if container_id in containers and signal:
        # Fake response for killing a container with a specific signal
        containers[container_id]["State"] = "stopped"
        containers[container_id]["Status"] = f"Stopped by signal {signal}"
        return '', 204
    return '', 400

@app.route('/v1.24/containers/<hex:container_id>/json', methods=['GET'])
def inspect_container(container_id):
    container = containers.get(container_id)
    if container:
        return jsonify(container), 200
    return '', 404

@app.route('/v1.24/containers/<hex:container_id>/exec', methods=['POST'])
def exec_create(container_id):
    if container_id in containers:
        exec_id = f"{container_id}-exec"
        response_body = {
            "Id": exec_id
        }
        return jsonify(response_body), 201
    return '', 404

@app.route('/v1.24/exec/<exec_id>/start', methods=['POST'])
def exec_start(exec_id):
    if exec_id.endswith('-exec'):
        return '', 200
    return '', 404

@app.route('/v1.24/images/<name>/push', methods=['POST'])
def push_image(name):
    tag = request.args.get('tag', 'latest')
    image_name = f"{name}:{tag}".replace("docker.io/", "")
    if image_name in images:
        # Simulate pushing the image
        response_body = {
            "status": "Pushed image " + image_name,
            "id": "sha256:1234567890ab"
        }
        response = make_response(jsonify(response_body), 200)
        response.headers['Content-Type'] = 'application/json'
        return response
    return '', 404


@app.route('/v1.24/containers/<container_id>/export', methods=['GET'])
def export_container(container_id):
    if container_id in containers:
        # Create a tarfile in memory
        tar_bytes = io.BytesIO()
        with tarfile.open(fileobj=tar_bytes, mode='w') as tar:
            # Here you would add files to the tar file
            # For demonstration, we're just adding a fake file
            fake_file_data = b"Hello, world!\n"
            tarinfo = tarfile.TarInfo(name="hello.txt")
            tarinfo.size = len(fake_file_data)
            tar.addfile(tarinfo, io.BytesIO(fake_file_data))
        
        tar_bytes.seek(0)
        
        response = send_file(
            tar_bytes,
            as_attachment=True,
            download_name=f"{container_id}.tar",
            mimetype='application/x-tar'
        )
        return response
    
    return '', 404




@app.errorhandler(404)
def page_not_found(e):
    return '', 404

if __name__ == '__main__':
    app.run(host='127.0.0.1', port=2375, debug=True)
    socketio.run(app, host='127.0.0.1', port=2375, debug=True)
