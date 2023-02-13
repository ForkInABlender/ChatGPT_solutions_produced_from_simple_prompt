"""
This is a flask application built with the purposes of imitating the kubernetes HTTP requests it relies on to function properly. Below is a example of how to approach
 if you're looking to imitate some of the kubernetes backend. The point being, if you try hard enough, anything imitated stops being flattery and more of a puzzle. "Am I
  or the mask not adequately solving the imitation-of-agency problem? Am I doing the wrong level of overthinking in any one area leading me to not see?"

This answer only provides a way to proof of how to do such. This probably also is a proof that python, nginx and openssl with a docker-in-docker container solution is
 all the developer realistically needs. This is the problem of setting up the honey-pot in reverse. Where the tool gets mock tested and the developer ends up exposing
  there entire hand as to what they're trying to do to the hacker. As a result, this also opens the door to the problem even moreso of "spy versus spy" problem and
 anything like it. This in turn leaves every user potentially exposed to any hacker's self hosted website using the same domain name. The spoof is unavoidable. It can
  not be patched. Even if you think rewriting the code would prevent that.

Energy[in] = Energy[out]
Garbage[in] = Garbage[out]

Play stupid games, win stupid prices. Being terminally stuck on stupid is not exclusionary for any human that plays stupid as well. The simplest way to win is to
 simply not play the game.
"""

from flask import Flask, jsonify, request

app = Flask(__name__)

pods = []

@app.route("/api/v1/pods", methods=["GET", "POST"])
def pods_handler():
    if request.method == "GET":
        return jsonify(pods)
    elif request.method == "POST":
        pod = request.get_json()
        pods.append(pod)
        return jsonify(pod), 201

@app.route("/api/v1/pods/<name>", methods=["GET", "DELETE"])
def pod_handler(name):
    for i, pod in enumerate(pods):
        if pod["metadata"]["name"] == name:
            if request.method == "GET":
                return jsonify(pod)
            elif request.method == "DELETE":
                del pods[i]
                return "", 204
    return "Pod not found", 404

if __name__ == "__main__":
    app.run(debug=True)
