# Dylan Kenneth Eliot & GPT-4-plugins

"""


This is a custom flask app I put together using code already produced from GPT-4-plugins by OpenAI. The app itself renders a 3-d cube, and
 gives it a set of coordinates to rotate about the plane. Thus far, all it does is rotate a cube with a single image mapped onto it.

Why use brython.js, THREE.js, and flask?
 Because it is simpler than the alternative. It also saves time when the only thing to worry about is functional and if it runs. 



"""


from flask import Flask
from flask_cors import CORS, cross_origin
from requests import get
app = Flask(__name__)
CORS(app)
@app.route("/*", methods=["OPTIONS"])
def star_OPTIONS():
	return "", 200

@app.route("/", methods=["GET"])
@cross_origin()
def index():
    return """<!DOCTYPE html>
<html>
<head>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/brython/3.11.3/brython.min.js"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/brython/3.11.3/brython_stdlib.min.js"></script>
</head>
<body onload="brython()">
    <script type="text/python">
        from browser import document, window, timer
        THREE = window.THREE
        scene = THREE.Scene.new()
        loader = THREE.TextureLoader.new()
        camera = THREE.PerspectiveCamera.new(75, window.innerWidth / window.innerHeight, 0.1, 1000)
        renderer = THREE.WebGLRenderer.new()
        renderer.setSize(window.innerWidth, window.innerHeight)
        document <= renderer.domElement
        geometry = THREE.BoxGeometry.new(1, 1, 1)
        texture = loader.load("https://cors-anywhere.herokuapp.com/https://illustoon.com/photo/7251.png")
        material = THREE.MeshBasicMaterial.new({"map": texture})
        cube = THREE.Mesh.new(geometry, material)
        scene.add(cube)
        camera.position.z = 5
        def animate(i):
            cube.rotation.x += 0.01
            cube.rotation.y += 0.01
            cube.rotation.z += 0.01
            renderer.render(scene, camera)
            timer.request_animation_frame(animate)
        animate(0)
    </script>
</body>
</html>"""
app.run(host="0.0.0.0")
