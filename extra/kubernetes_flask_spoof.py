#Dylan Kenneth Eliot

""""
This file spoofs kubernetes functionality.


in this way, it both breaks the default convention of kubernetes development, making it able to equally function on Android.

due to type marshalling of information, it matches just enough of what kubectl expects in normal kubernetes environment minimally.


This file was created on 3/3/2026 @ 11:49 am. And retested on android via termux for quality assurance.

Professional development will continue within 6 hrs from 2:06 on 3/4/2026;;; code patched 11:51 pm

#### ADDENDUM ####

Some corrections and additions were made by claude.ai & https://perchance.org/ai-text-generator as together on freemium compute they can break 
 the kubernetes stack tools. If fleshed out further, it could mean kubernetes becomes obsolete. Like Docker engine aka dockerd, if it cannot run
  on a cellphone, it is not good enough for production release.

Further development of this server will continue until full server spoof is complete.
""""



# fake_kube_cluster_fixed.py
from flask import Flask, request, jsonify, Response
import uuid, copy, json
from datetime import datetime
import time, random

app = Flask(__name__)

# -----------------------------
# Global state
# -----------------------------
store = {
    "namespaces": {},
    "pods": {},
    "services": {},
    "deployments": {},
    "replicasets": {},
    "configmaps": {},
    "secrets": {},
    "roles": {},
    "rolebindings": {},
    "clusterroles": {},
    "clusterrolebindings": {},
    "nodes": {},
    "storageclasses": {}
}

GLOBAL_RV = 0

# -----------------------------
# Utilities
# -----------------------------
def next_rv():
    global GLOBAL_RV
    GLOBAL_RV += 1
    return str(GLOBAL_RV)

def now():
    return datetime.utcnow().isoformat() + "Z"

def assign_metadata(obj, namespace=None):
    meta = obj.setdefault("metadata", {})
    meta.setdefault("uid", str(uuid.uuid4()))
    meta.setdefault("creationTimestamp", now())
    if namespace:
        meta["namespace"] = namespace
    meta["resourceVersion"] = next_rv()
    meta.setdefault("generation", 1)
    return obj

def make_list(kind, items):
    return {
        "kind": f"{kind}List",
        "apiVersion": "v1" if kind in ["Pod","Service","Namespace","ConfigMap","Secret","Node","StorageClass"] else "apps/v1",
        "metadata": {"resourceVersion": next_rv()},
        "items": items
    }

def deep_merge(original, patch):
    for k, v in patch.items():
        if isinstance(v, dict) and isinstance(original.get(k), dict):
            deep_merge(original[k], v)
        else:
            original[k] = v
    return original

# -----------------------------
# Discovery Endpoints
# -----------------------------
@app.route("/version")
def version():
    return jsonify({"major":"1","minor":"28","gitVersion":"v1.28.0-fake"})

@app.route("/api")
def api_root():
    return jsonify({"kind":"APIVersions","versions":["v1"],"serverAddressByClientCIDRs":[]})

@app.route("/apis")
def apis_root():
    return jsonify({"kind":"APIGroupList","groups":[
        {"name":"apps","versions":[{"groupVersion":"apps/v1","version":"v1"}],"preferredVersion":{"groupVersion":"apps/v1","version":"v1"}},
        {"name":"rbac.authorization.k8s.io","versions":[{"groupVersion":"rbac.authorization.k8s.io/v1","version":"v1"}],"preferredVersion":{"groupVersion":"rbac.authorization.k8s.io/v1","version":"v1"}}
    ]})

@app.route("/api/v1")
def api_v1():
    return jsonify({"kind":"APIResourceList","groupVersion":"v1","resources":[
        {"name":"pods","namespaced":True,"kind":"Pod"},
        {"name":"services","namespaced":True,"kind":"Service"},
        {"name":"namespaces","namespaced":False,"kind":"Namespace"},
        {"name":"configmaps","namespaced":True,"kind":"ConfigMap"},
        {"name":"secrets","namespaced":True,"kind":"Secret"},
        {"name":"nodes","namespaced":False,"kind":"Node"},
        {"name":"storageclasses","namespaced":False,"kind":"StorageClass"}
    ]})

@app.route("/apis/apps/v1")
def apps_v1():
    return jsonify({"kind":"APIResourceList","groupVersion":"apps/v1","resources":[
        {"name":"deployments","namespaced":True,"kind":"Deployment"},
        {"name":"replicasets","namespaced":True,"kind":"ReplicaSet"}
    ]})

@app.route("/apis/rbac.authorization.k8s.io/v1")
def rbac_v1():
    return jsonify({"kind":"APIResourceList","groupVersion":"rbac.authorization.k8s.io/v1","resources":[
        {"name":"roles","namespaced":True,"kind":"Role"},
        {"name":"rolebindings","namespaced":True,"kind":"RoleBinding"},
        {"name":"clusterroles","namespaced":False,"kind":"ClusterRole"},
        {"name":"clusterrolebindings","namespaced":False,"kind":"ClusterRoleBinding"}
    ]})

@app.route("/openapi/v2")
def openapi():
    return jsonify({})  # stub

# -----------------------------
# Namespace endpoints
# -----------------------------
def ensure_namespace(name):
    if name not in store["namespaces"]:
        ns = {"metadata":{"name":name}}
        assign_metadata(ns)
        store["namespaces"][name] = ns

@app.route("/api/v1/namespaces", methods=["GET","POST"])
def namespaces():
    if request.method=="GET":
        return jsonify(make_list("Namespace", list(store["namespaces"].values())))
    data=request.json
    assign_metadata(data)
    name=data["metadata"]["name"]
    store["namespaces"][name]=data
    return jsonify(data)

@app.route("/api/v1/namespaces/<name>", methods=["GET","DELETE"])
def namespace_detail(name):
    if request.method=="GET":
        return jsonify(store["namespaces"].get(name, {}))
    store["namespaces"].pop(name,None)
    return jsonify({"status":"Success"})

# -----------------------------
# Enhanced Status Handling
# -----------------------------
def simulate_pod_status(pod_name):
    """Generates realistic pod status with simulated conditions"""
    return {
        "phase": "Running",
        "conditions": [
            {"type": "Initialized", "status": "True"},
            {"type": "Ready", "status": "True"},
            {"type": "ContainersReady", "status": "True"},
            {"type": "PodScheduled", "status": "True"}
        ],
        "containerStatuses": [{
            "name": "main",
            "state": {"running": {"startedAt": now()}},
            "ready": True,
            "restartCount": 0,
            "image": "alpine:latest",
            "imageID": f"docker-pullable://alpine@sha256:{str(uuid.uuid4())[:32]}",
            "containerID": f"docker://{str(uuid.uuid4())}"
        }],
        "startTime": now(),
        "podIP": f"10.244.{random.randint(0,255)}.{random.randint(1,254)}",
        "qosClass": "Burstable"
    }

# -----------------------------
# Enhanced Pod Creation
# -----------------------------
@app.route("/api/v1/namespaces/<ns>/pods", methods=["POST"])
def create_pod(ns):
    ensure_namespace(ns)
    data = request.json
    assign_metadata(data, ns)
    
    # Auto-populate required fields if missing
    data.setdefault("spec", {})
    data["spec"].setdefault("containers", [{
        "name": "main",
        "image": "alpine:latest",
        "imagePullPolicy": "IfNotPresent"
    }])
    
    # Generate realistic status
    data["status"] = simulate_pod_status(data["metadata"]["name"])
    
    key = f"{ns}/{data['metadata']['name']}"
    store["pods"][key] = data
    
    # Simulate scheduling delay
    time.sleep(0.5)
    
    return jsonify(data)

# -----------------------------
# Enhanced Watch Endpoint
# -----------------------------
@app.route("/api/v1/watch/<path:path>")
def watch_resource(path):
    def event_generator():
        last_rv = request.args.get('resourceVersion', '0')
        
        while True:
            current_rv = GLOBAL_RV
            if int(last_rv) < current_rv:
                # Find changed resources
                changed = []
                for resource_type in store:
                    for key, obj in store[resource_type].items():
                        if obj.get("metadata", {}).get("resourceVersion", "0") > last_rv:
                            changed.append(("MODIFIED", obj))
                
                # Sort by resourceVersion
                changed.sort(key=lambda x: x[1]["metadata"]["resourceVersion"])
                
                for event_type, obj in changed:
                    yield f"data: {json.dumps({'type': event_type, 'object': obj})}\n\n"
                
                last_rv = str(current_rv)
            
            time.sleep(1)
    
    return Response(event_generator(), mimetype="text/event-stream")

# -----------------------------
# Metrics Endpoint
# -----------------------------
@app.route("/metrics")
def metrics():
    metrics = [
        "# HELP kubernetes_mock_resources_total Total number of resources in mock cluster",
        "# TYPE kubernetes_mock_resources_total gauge"
    ]
    
    for resource_type, resources in store.items():
        metrics.append(f'kubernetes_mock_resources_total{{resource_type="{resource_type}"}} {len(resources)}')
    
    return Response("\n".join(metrics), mimetype="text/plain")

# -----------------------------
# Enhanced Health Checks
# -----------------------------
@app.route("/healthz")
def health_check():
    return jsonify({
        "status": "healthy",
        "components": {
            "api": True,
            "etcd": False,  # Mock doesn't use etcd
            "controller": True,
            "scheduler": True
        }
    })

# -----------------------------
# Enhanced Error Handling
# -----------------------------
@app.errorhandler(404)
def not_found(e):
    return jsonify({
        "kind": "Status",
        "apiVersion": "v1",
        "metadata": {},
        "status": "Failure",
        "message": str(e),
        "reason": "NotFound",
        "details": {},
        "code": 404
    }), 404

# -----------------------------
# Generic Namespaced Resource Handler
# -----------------------------
def create_ns_handler(resource_key, kind, default_status=None):
    def handler(ns):
        ensure_namespace(ns)
        if request.method=="GET":
            items=[v for k,v in store[resource_key].items() if k.startswith(ns+"/")]
            return jsonify(make_list(kind, items))
        data=request.json
        assign_metadata(data, ns)
        if default_status:
            data.setdefault("status", copy.deepcopy(default_status))
        key=f"{ns}/{data['metadata']['name']}"
        store[resource_key][key]=data
        return jsonify(data)
    handler.__name__ = f"{resource_key}_list"
    return handler

def create_ns_detail_handler(resource_key):
    def handler(ns,name):
        key=f"{ns}/{name}"
        if request.method=="GET":
            return jsonify(store[resource_key].get(key, {}))
        if request.method=="DELETE":
            store[resource_key].pop(key,None)
            return jsonify({"status":"Success"})
        if request.method in ["PUT","PATCH"]:
            obj=store[resource_key].get(key)
            if not obj:
                return jsonify({"error":"Not found"}),404
            patch=request.json
            deep_merge(obj,patch)
            obj["metadata"]["resourceVersion"]=next_rv()
            return jsonify(obj)
    handler.__name__ = f"{resource_key}_detail"
    return handler

# -----------------------------
# Pod and Service endpoints
# -----------------------------
pod_status={"phase":"Running","conditions":[{"type":"Ready","status":"True"}],"containerStatuses":[]}
app.add_url_rule("/api/v1/namespaces/<ns>/pods", view_func=create_ns_handler("pods","Pod",pod_status), methods=["GET","POST"], endpoint="pods_list")
app.add_url_rule("/api/v1/namespaces/<ns>/pods/<name>", view_func=create_ns_detail_handler("pods"), methods=["GET","PUT","PATCH","DELETE"], endpoint="pods_detail")
app.add_url_rule("/api/v1/namespaces/<ns>/services", view_func=create_ns_handler("services","Service"), methods=["GET","POST"], endpoint="services_list")
app.add_url_rule("/api/v1/namespaces/<ns>/services/<name>", view_func=create_ns_detail_handler("services"), methods=["GET","PUT","PATCH","DELETE"], endpoint="services_detail")
app.add_url_rule("/api/v1/namespaces/<ns>/configmaps", view_func=create_ns_handler("configmaps","ConfigMap"), methods=["GET","POST"], endpoint="configmaps_list")
app.add_url_rule("/api/v1/namespaces/<ns>/configmaps/<name>", view_func=create_ns_detail_handler("configmaps"), methods=["GET","PUT","PATCH","DELETE"], endpoint="configmaps_detail")
app.add_url_rule("/api/v1/namespaces/<ns>/secrets", view_func=create_ns_handler("secrets","Secret"), methods=["GET","POST"], endpoint="secrets_list")
app.add_url_rule("/api/v1/namespaces/<ns>/secrets/<name>", view_func=create_ns_detail_handler("secrets"), methods=["GET","PUT","PATCH","DELETE"], endpoint="secrets_detail")

# -----------------------------
# Deployment and ReplicaSet endpoints
# -----------------------------
deployment_status={"replicas":1,"readyReplicas":1,"availableReplicas":1,"observedGeneration":1}
app.add_url_rule("/apis/apps/v1/namespaces/<ns>/deployments", view_func=create_ns_handler("deployments","Deployment",deployment_status), methods=["GET","POST"], endpoint="deployments_list")
app.add_url_rule("/apis/apps/v1/namespaces/<ns>/deployments/<name>", view_func=create_ns_detail_handler("deployments"), methods=["GET","PATCH","DELETE"], endpoint="deployments_detail")
app.add_url_rule("/apis/apps/v1/namespaces/<ns>/replicasets", view_func=create_ns_handler("replicasets","ReplicaSet",deployment_status), methods=["GET","POST"], endpoint="replicasets_list")
app.add_url_rule("/apis/apps/v1/namespaces/<ns>/replicasets/<name>", view_func=create_ns_detail_handler("replicasets"), methods=["GET","PATCH","DELETE"], endpoint="replicasets_detail")

# -----------------------------
# RBAC endpoints
# -----------------------------
app.add_url_rule("/apis/rbac.authorization.k8s.io/v1/namespaces/<ns>/roles", view_func=create_ns_handler("roles","Role"), methods=["GET","POST"], endpoint="roles_list")
app.add_url_rule("/apis/rbac.authorization.k8s.io/v1/namespaces/<ns>/roles/<name>", view_func=create_ns_detail_handler("roles"), methods=["GET","PUT","PATCH","DELETE"], endpoint="roles_detail")
app.add_url_rule("/apis/rbac.authorization.k8s.io/v1/namespaces/<ns>/rolebindings", view_func=create_ns_handler("rolebindings","RoleBinding"), methods=["GET","POST"], endpoint="rolebindings_list")
app.add_url_rule("/apis/rbac.authorization.k8s.io/v1/namespaces/<ns>/rolebindings/<name>", view_func=create_ns_detail_handler("rolebindings"), methods=["GET","PUT","PATCH","DELETE"], endpoint="rolebindings_detail")
app.add_url_rule("/apis/rbac.authorization.k8s.io/v1/clusterroles", view_func=create_ns_handler("clusterroles","ClusterRole"), methods=["GET","POST"], endpoint="clusterroles_list")
app.add_url_rule("/apis/rbac.authorization.k8s.io/v1/clusterroles/<name>", view_func=create_ns_detail_handler("clusterroles"), methods=["GET","PUT","PATCH","DELETE"], endpoint="clusterroles_detail")
app.add_url_rule("/apis/rbac.authorization.k8s.io/v1/clusterrolebindings", view_func=create_ns_handler("clusterrolebindings","ClusterRoleBinding"), methods=["GET","POST"], endpoint="clusterrolebindings_list")
app.add_url_rule("/apis/rbac.authorization.k8s.io/v1/clusterrolebindings/<name>", view_func=create_ns_detail_handler("clusterrolebindings"), methods=["GET","PUT","PATCH","DELETE"], endpoint="clusterrolebindings_detail")

# -----------------------------
# Watch endpoint stub
# -----------------------------
@app.route("/api/v1/namespaces/<ns>/pods/watch")
def watch_pods(ns):
    def generate():
        for v in store["pods"].values():
            yield json.dumps({"type":"ADDED","object":v}) + "\n"
    return Response(generate(), mimetype="application/json")

# -----------------------------
# Cluster-level stubs
# -----------------------------
# kube-system namespace
ensure_namespace("kube-system")

# Fake control-plane pods
for pod_name in ["kube-apiserver","kube-controller-manager","kube-scheduler","coredns"]:
    key=f"kube-system/{pod_name}"
    if key not in store["pods"]:
        pod={"metadata":{"name":pod_name},"spec":{},"status":{"phase":"Running"}}
        assign_metadata(pod,"kube-system")
        store["pods"][key]=pod

# kube-dns service
if "kube-system/kube-dns" not in store["services"]:
    svc={"metadata":{"name":"kube-dns"},"spec":{"clusterIP":"10.96.0.10"}}
    assign_metadata(svc,"kube-system")
    store["services"]["kube-system/kube-dns"]=svc

# Fake node
if "fake-node" not in store["nodes"]:
    store["nodes"]["fake-node"]={"metadata":{"name":"fake-node"},"status":{"conditions":[{"type":"Ready","status":"True"}]}}

# ClusterRole/Binding
if "cluster-admin" not in store["clusterroles"]:
    cr={"metadata":{"name":"cluster-admin"}}
    assign_metadata(cr)
    store["clusterroles"]["cluster-admin"]=cr
if "cluster-admin-binding" not in store["clusterrolebindings"]:
    crb={"metadata":{"name":"cluster-admin-binding"},"roleRef":{"kind":"ClusterRole","name":"cluster-admin","apiGroup":"rbac.authorization.k8s.io"}}
    assign_metadata(crb)
    store["clusterrolebindings"]["cluster-admin-binding"]=crb

# Minimal storageclass
if "standard" not in store["storageclasses"]:
    sc={"metadata":{"name":"standard"}}
    assign_metadata(sc)
    store["storageclasses"]["standard"]=sc
@app.route("/apis/storage.k8s.io/v1/storageclasses")
def storage_classes():
    return jsonify(make_list("StorageClass", list(store["storageclasses"].values())))

# Fake nodes endpoint
@app.route("/api/v1/nodes")
def nodes():
    return jsonify(make_list("Node", list(store["nodes"].values())))

# -----------------------------
# Main
# -----------------------------
if __name__=="__main__":
    ensure_namespace("default")
    ensure_namespace("kube-system")
    # Seed an initial pod directly into the store — do NOT call create_pod()
    # here because it is a Flask route handler and requires an active HTTP
    # request context (request.json would raise RuntimeError otherwise).
    _seed_pod_name = "default-pod"
    _seed_key = f"default/{_seed_pod_name}"
    if _seed_key not in store["pods"]:
        _seed_pod = {
            "metadata": {"name": _seed_pod_name},
            "spec": {
                "containers": [{
                    "name": "main",
                    "image": "alpine:latest",
                    "imagePullPolicy": "IfNotPresent"
                }]
            }
        }
        assign_metadata(_seed_pod, "default")
        _seed_pod["status"] = simulate_pod_status(_seed_pod_name)
        store["pods"][_seed_key] = _seed_pod
    app.run(host="0.0.0.0", port=8080, debug=True)
