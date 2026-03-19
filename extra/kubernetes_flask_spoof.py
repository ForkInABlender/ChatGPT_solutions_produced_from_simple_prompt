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

Product is not useable for those aiming to 'fake it to make it'.

#### END OF ADDENDUM ####

updated 9:51:34 am on 03/05/2026

#### ADDENDUM ####

Kubernetes is now flexibly broken, moreover than before.

These patches were added by "cursor.ai", a even nummer than life toy thing of the banana.

Now kubernetes is flattened, and given ``kubectl --server="http://localhost:8080"`` with this server running in another process, root permissions
 no longer matter....
#### END OF ADDENDUM ####

updated 07:30:11 pm on 03/18/2026

""""

#!/usr/bin/env python3
"""
Expanded Kubernetes mock server (Flask).

This is intended for local/testing scenarios where you want a lightweight API
surface that looks like Kubernetes enough for basic kubectl/client-go usage.
"""

from __future__ import annotations

import copy
import json
import random
import threading
import time
import uuid
from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Any, Dict, Iterable, List, Optional, Tuple

from flask import Flask, Response, jsonify, request

app = Flask(__name__)


def _utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


GLOBAL_RV_LOCK = threading.Lock()
GLOBAL_RV = 0


def next_rv() -> str:
    global GLOBAL_RV
    with GLOBAL_RV_LOCK:
        GLOBAL_RV += 1
        return str(GLOBAL_RV)


def make_uid() -> str:
    return str(uuid.uuid4())


def deep_merge(original: Dict[str, Any], patch: Dict[str, Any]) -> Dict[str, Any]:
    """
    JSON-merge style deep merge for dicts.

    This is not a full Kubernetes strategic/JSON patch implementation, but is
    good enough for a local mock.
    """

    for k, v in patch.items():
        if isinstance(v, dict) and isinstance(original.get(k), dict):
            deep_merge(original[k], v)
        else:
            original[k] = v
    return original


def _metadata_get_name(obj: Dict[str, Any]) -> Optional[str]:
    meta = obj.get("metadata") or {}
    return meta.get("name")


def _metadata_get_generate_name(obj: Dict[str, Any]) -> Optional[str]:
    meta = obj.get("metadata") or {}
    return meta.get("generateName")


def assign_metadata(
    obj: Dict[str, Any],
    *,
    namespace: Optional[str],
    name: str,
) -> Dict[str, Any]:
    meta = obj.setdefault("metadata", {})
    meta.setdefault("uid", make_uid())
    meta.setdefault("creationTimestamp", _utc_now_iso())
    meta["name"] = name
    if namespace is not None:
        meta["namespace"] = namespace
    meta["resourceVersion"] = next_rv()
    meta.setdefault("generation", 1)
    return obj


def make_list(kind: str, api_version: str, items: List[Dict[str, Any]]) -> Dict[str, Any]:
    # Kubernetes list objects have kind "<Kind>List"
    # and metadata.resourceVersion generally reflects the state at list time.
    return {
        "kind": f"{kind}List",
        "apiVersion": api_version,
        "metadata": {"resourceVersion": str(GLOBAL_RV)},
        "items": items,
    }


def sse_encode(event: Dict[str, Any]) -> str:
    # Standard Kubernetes watch response uses SSE with "data: <json>\n\n"
    return f"data: {json.dumps(event, separators=(',', ':'))}\n\n"


@dataclass(frozen=True)
class ResourceSpec:
    kind: str
    api_version: str  # e.g. "v1" or "apps/v1"
    plural: str
    namespaced: bool


CORE = "v1"

RESOURCE_SPECS: Dict[str, ResourceSpec] = {
    # Core group
    "Pod": ResourceSpec(kind="Pod", api_version="v1", plural="pods", namespaced=True),
    "Service": ResourceSpec(kind="Service", api_version="v1", plural="services", namespaced=True),
    "Namespace": ResourceSpec(kind="Namespace", api_version="v1", plural="namespaces", namespaced=False),
    "ConfigMap": ResourceSpec(kind="ConfigMap", api_version="v1", plural="configmaps", namespaced=True),
    "Secret": ResourceSpec(kind="Secret", api_version="v1", plural="secrets", namespaced=True),
    "Node": ResourceSpec(kind="Node", api_version="v1", plural="nodes", namespaced=False),
    "StorageClass": ResourceSpec(kind="StorageClass", api_version="v1", plural="storageclasses", namespaced=False),
    # apps group
    "Deployment": ResourceSpec(kind="Deployment", api_version="apps/v1", plural="deployments", namespaced=True),
    "ReplicaSet": ResourceSpec(kind="ReplicaSet", api_version="apps/v1", plural="replicasets", namespaced=True),
    # rbac group
    "Role": ResourceSpec(kind="Role", api_version="rbac.authorization.k8s.io/v1", plural="roles", namespaced=True),
    "RoleBinding": ResourceSpec(
        kind="RoleBinding", api_version="rbac.authorization.k8s.io/v1", plural="rolebindings", namespaced=True
    ),
    "ClusterRole": ResourceSpec(
        kind="ClusterRole", api_version="rbac.authorization.k8s.io/v1", plural="clusterroles", namespaced=False
    ),
    "ClusterRoleBinding": ResourceSpec(
        kind="ClusterRoleBinding",
        api_version="rbac.authorization.k8s.io/v1",
        plural="clusterrolebindings",
        namespaced=False,
    ),
}


def api_base_for(api_version: str) -> str:
    # core: "v1" => /api/v1
    if api_version == "v1":
        return "/api/v1"
    group, version = api_version.split("/", 1)
    return f"/apis/{group}/{version}"


def sse_watch_snapshot(items: List[Dict[str, Any]], rv: str) -> None:
    # Not used; kept for clarity if you later want to add "ADDED" baseline.
    _ = (items, rv)


# -----------------------------
# In-memory "cluster" state
# -----------------------------
store: Dict[str, Dict[str, Dict[str, Any]]] = {
    kind: {} for kind in RESOURCE_SPECS.keys()
}
namespaces: Dict[str, Dict[str, Any]] = {}
event_log: List[Dict[str, Any]] = []
EVENT_LOCK = threading.Lock()


def _key_for(spec: ResourceSpec, *, namespace: Optional[str], name: str) -> str:
    if spec.namespaced:
        if namespace is None:
            raise ValueError(f"namespace required for {spec.kind}")
        return f"{namespace}/{name}"
    return name


def list_resources(spec: ResourceSpec, *, namespace: Optional[str] = None) -> List[Dict[str, Any]]:
    if spec.namespaced:
        if namespace is None:
            return [obj for obj in store[spec.kind].values()]
        prefix = f"{namespace}/"
        return [obj for k, obj in store[spec.kind].items() if k.startswith(prefix)]
    return list(store[spec.kind].values())


def get_resource(spec: ResourceSpec, *, namespace: Optional[str], name: str) -> Optional[Dict[str, Any]]:
    key = _key_for(spec, namespace=namespace, name=name)
    return store[spec.kind].get(key)


def push_event(event_type: str, obj: Dict[str, Any], rv: Optional[str] = None) -> None:
    if rv is None:
        rv = str(GLOBAL_RV)
    with EVENT_LOCK:
        event_log.append({"type": event_type, "rv": rv, "object": copy.deepcopy(obj)})


def simulate_pod_status(pod_name: str) -> Dict[str, Any]:
    def rand_ip() -> str:
        return f"10.244.{random.randint(0,255)}.{random.randint(1,254)}"

    return {
        "phase": "Running",
        "conditions": [
            {"type": "Initialized", "status": "True"},
            {"type": "Ready", "status": "True"},
            {"type": "ContainersReady", "status": "True"},
            {"type": "PodScheduled", "status": "True"},
        ],
        "containerStatuses": [
            {
                "name": "main",
                "state": {"running": {"startedAt": _utc_now_iso()}},
                "ready": True,
                "restartCount": 0,
                "image": "alpine:latest",
                "imageID": f"docker-pullable://alpine@sha256:{uuid.uuid4().hex[:32]}",
                "containerID": f"docker://{uuid.uuid4().hex}",
            }
        ],
        "startTime": _utc_now_iso(),
        "podIP": rand_ip(),
        "qosClass": "Burstable",
    }


def simulate_deployment_status() -> Dict[str, Any]:
    return {
        "replicas": 1,
        "readyReplicas": 1,
        "availableReplicas": 1,
        "observedGeneration": 1,
    }


def ensure_namespace(ns: str) -> None:
    if ns in namespaces:
        return
    obj = {"metadata": {"name": ns}}
    obj = assign_metadata(obj, namespace=None, name=ns)
    namespaces[ns] = obj
    # Namespace itself is also a resource of kind Namespace
    store["Namespace"][_key_for(RESOURCE_SPECS["Namespace"], namespace=None, name=ns)] = obj
    push_event("ADDED", obj, rv=str(GLOBAL_RV))


# -----------------------------
# Discovery endpoints
# -----------------------------
@app.route("/version")
def version():
    return jsonify({"major": "1", "minor": "28", "gitVersion": "v1.28.0-fake"})


@app.route("/api")
def api_root():
    return jsonify({"kind": "APIVersions", "versions": ["v1"], "serverAddressByClientCIDRs": []})


@app.route("/apis")
def apis_root():
    return jsonify(
        {
            "kind": "APIGroupList",
            "groups": [
                {
                    "name": "apps",
                    "versions": [{"groupVersion": "apps/v1", "version": "v1"}],
                    "preferredVersion": {"groupVersion": "apps/v1", "version": "v1"},
                },
                {
                    "name": "rbac.authorization.k8s.io",
                    "versions": [{"groupVersion": "rbac.authorization.k8s.io/v1", "version": "v1"}],
                    "preferredVersion": {"groupVersion": "rbac.authorization.k8s.io/v1", "version": "v1"},
                },
                {
                    "name": "storage.k8s.io",
                    "versions": [{"groupVersion": "storage.k8s.io/v1", "version": "v1"}],
                    "preferredVersion": {"groupVersion": "storage.k8s.io/v1", "version": "v1"},
                },
            ],
        }
    )


@app.route("/api/v1")
def core_resources():
    # Kubectl's "all" shortcut relies on resource discovery, including which verbs
    # are supported. If `verbs` is omitted, kubectl may treat matching resources
    # as not listable and refuse shortcuts like `kubectl get all`.
    common_get_list_watch = ["get", "list", "watch"]
    common_write = ["create", "update", "patch", "delete"]
    core_verbs = {
        "pods": common_get_list_watch + common_write,
        "services": common_get_list_watch + common_write,
        "namespaces": common_get_list_watch + common_write,
        "configmaps": common_get_list_watch + common_write,
        "secrets": common_get_list_watch + common_write,
        "nodes": common_get_list_watch,
        "storageclasses": common_get_list_watch,
        # Pseudo-resource used so kubectl's `get all` can succeed even though
        # our mock doesn't implement kubectl's full client-side shortcut logic.
        "all": common_get_list_watch,
    }
    return jsonify(
        {
            "kind": "APIResourceList",
            "groupVersion": "v1",
            "resources": [
                {"name": "all", "namespaced": True, "kind": "All", "verbs": core_verbs["all"]},
                {"name": "pods", "namespaced": True, "kind": "Pod", "verbs": core_verbs["pods"]},
                {
                    "name": "services",
                    "namespaced": True,
                    "kind": "Service",
                    "verbs": core_verbs["services"],
                },
                {
                    "name": "namespaces",
                    "namespaced": False,
                    "kind": "Namespace",
                    "verbs": core_verbs["namespaces"],
                },
                {
                    "name": "configmaps",
                    "namespaced": True,
                    "kind": "ConfigMap",
                    "verbs": core_verbs["configmaps"],
                },
                {"name": "secrets", "namespaced": True, "kind": "Secret", "verbs": core_verbs["secrets"]},
                {"name": "nodes", "namespaced": False, "kind": "Node", "verbs": core_verbs["nodes"]},
                {
                    "name": "storageclasses",
                    "namespaced": False,
                    "kind": "StorageClass",
                    "verbs": core_verbs["storageclasses"],
                },
            ],
        }
    )


def _all_items_for_namespace(namespace: Optional[str]) -> List[Dict[str, Any]]:
    items: List[Dict[str, Any]] = []
    # This is a pragmatic approximation of `kubectl get all`'s intent.
    wanted = ["Pod", "Service", "Deployment", "ReplicaSet"]
    for kind in wanted:
        spec_kind_items = store[kind]
        if namespace is None:
            items.extend(copy.deepcopy(list(spec_kind_items.values())))
        else:
            prefix = f"{namespace}/"
            for k, obj in spec_kind_items.items():
                if k.startswith(prefix):
                    items.append(copy.deepcopy(obj))
    return items


@app.route("/api/v1/all", methods=["GET"])
def all_resources_all_namespaces():
    # Cluster-wide list of the pseudo-resource `all`.
    items = _all_items_for_namespace(namespace=None)
    return jsonify(make_list("All", "v1", items))


@app.route("/api/v1/namespaces/<ns>/all", methods=["GET"])
def all_resources_in_namespace(ns: str):
    items = _all_items_for_namespace(namespace=ns)
    return jsonify(make_list("All", "v1", items))


@app.route("/apis/apps/v1")
def apps_resources():
    common_get_list_watch = ["get", "list", "watch"]
    common_write = ["create", "update", "patch", "delete"]
    return jsonify(
        {
            "kind": "APIResourceList",
            "groupVersion": "apps/v1",
            "resources": [
                {
                    "name": "deployments",
                    "namespaced": True,
                    "kind": "Deployment",
                    "verbs": common_get_list_watch + common_write,
                },
                {
                    "name": "replicasets",
                    "namespaced": True,
                    "kind": "ReplicaSet",
                    "verbs": common_get_list_watch + common_write,
                },
            ],
        }
    )


@app.route("/apis/rbac.authorization.k8s.io/v1")
def rbac_resources():
    common_get_list_watch = ["get", "list", "watch"]
    common_write = ["create", "update", "patch", "delete"]
    return jsonify(
        {
            "kind": "APIResourceList",
            "groupVersion": "rbac.authorization.k8s.io/v1",
            "resources": [
                {"name": "roles", "namespaced": True, "kind": "Role", "verbs": common_get_list_watch + common_write},
                {
                    "name": "rolebindings",
                    "namespaced": True,
                    "kind": "RoleBinding",
                    "verbs": common_get_list_watch + common_write,
                },
                {
                    "name": "clusterroles",
                    "namespaced": False,
                    "kind": "ClusterRole",
                    "verbs": common_get_list_watch + common_write,
                },
                {
                    "name": "clusterrolebindings",
                    "namespaced": False,
                    "kind": "ClusterRoleBinding",
                    "verbs": common_get_list_watch + common_write,
                },
            ],
        }
    )


@app.route("/apis/storage.k8s.io/v1")
def storage_resources():
    common_get_list_watch = ["get", "list", "watch"]
    # Some kubectl flows query the storage group for StorageClass discovery.
    return jsonify(
        {
            "kind": "APIResourceList",
            "groupVersion": "storage.k8s.io/v1",
            "resources": [
                {"name": "storageclasses", "namespaced": False, "kind": "StorageClass", "verbs": common_get_list_watch}
            ],
        }
    )


@app.route("/openapi/v2")
def openapi_stub():
    # Many kubectl flows don't require this.
    return jsonify({})


# -----------------------------
# Error handling
# -----------------------------
@app.errorhandler(404)
def not_found(e):
    return jsonify(
        {
            "kind": "Status",
            "apiVersion": "v1",
            "metadata": {},
            "status": "Failure",
            "message": str(e),
            "reason": "NotFound",
            "details": {},
            "code": 404,
        }
    ), 404


# -----------------------------
# Metrics + health
# -----------------------------
@app.route("/healthz")
def healthz():
    return jsonify(
        {
            "status": "healthy",
            "components": {"api": True, "etcd": False, "controller": True, "scheduler": True},
        }
    )


@app.route("/metrics")
def metrics():
    lines = [
        "# HELP kubernetes_mock_resources_total Total number of resources in mock cluster",
        "# TYPE kubernetes_mock_resources_total gauge",
    ]
    for kind, objs in store.items():
        lines.append(f'kubernetes_mock_resources_total{{resource_type="{kind}"}} {len(objs)}')
    return Response("\n".join(lines), mimetype="text/plain")


# -----------------------------
# Watch implementation
# -----------------------------
def watch_stream(*, spec: ResourceSpec, namespace: Optional[str], resource_version: str) -> Iterable[str]:
    # Kubernetes uses resourceVersion as a monotonic integer-ish string.
    try:
        start_rv = int(resource_version)
    except Exception:
        start_rv = 0

    timeout_s = request.args.get("timeoutSeconds", type=float)
    deadline = time.time() + timeout_s if timeout_s is not None else None

    cursor = 0
    with EVENT_LOCK:
        # Move cursor to first event with rv > start_rv
        while cursor < len(event_log) and int(event_log[cursor]["rv"]) <= start_rv:
            cursor += 1

    while True:
        if deadline is not None and time.time() > deadline:
            return

        # Grab newly appended events
        with EVENT_LOCK:
            new_events = event_log[cursor:]
            cursor = len(event_log)

        for ev in new_events:
            obj = ev["object"]
            obj_ns = obj.get("metadata", {}).get("namespace")
            if spec.namespaced and namespace is not None and obj_ns != namespace:
                continue
            if (not spec.namespaced) and namespace is not None:
                # Should not happen.
                continue
            yield sse_encode({"type": ev["type"], "object": obj})

        # Poll interval: keep lightweight.
        time.sleep(0.25)


def _get_namespace_from_request(spec: ResourceSpec) -> Optional[str]:
    if not spec.namespaced:
        return None
    # For list/watch endpoints, kubectl typically uses URL namespace, not body.
    # We'll try URL first (view functions pass it in via closure) and fallback to body.
    return None


# -----------------------------
# CRUD handlers (generic)
# -----------------------------
def _ensure_correct_kind_fields(spec: ResourceSpec, obj: Dict[str, Any]) -> Dict[str, Any]:
    # Kubernetes clients commonly expect these fields.
    obj.setdefault("kind", spec.kind)
    obj["apiVersion"] = spec.api_version
    return obj


def _default_object_for_create(spec: ResourceSpec, obj: Dict[str, Any], namespace: Optional[str], name: str) -> Dict[str, Any]:
    obj = _ensure_correct_kind_fields(spec, obj)
    if spec.kind == "Pod":
        obj.setdefault("spec", {})
        obj["spec"].setdefault(
            "containers",
            [{"name": "main", "image": "alpine:latest", "imagePullPolicy": "IfNotPresent"}],
        )
        obj.setdefault("status", simulate_pod_status(name))
        # Mark a sensible nodeName, etc. (optional)
    elif spec.kind == "Service":
        obj.setdefault("spec", {})
        # Basic clusterIP behavior: if absent, assign an RFC1918-ish one.
        obj["spec"].setdefault("clusterIP", f"10.96.0.{random.randint(2, 254)}")
        obj.setdefault("status", {"loadBalancer": {}})
    elif spec.kind == "Deployment":
        obj.setdefault("spec", {})
        obj.setdefault("status", simulate_deployment_status())
    elif spec.kind == "ReplicaSet":
        obj.setdefault("spec", {})
        obj.setdefault("status", simulate_deployment_status())
    elif spec.kind == "Namespace":
        obj.setdefault("status", {})
    elif spec.kind == "Node":
        obj.setdefault("status", {"conditions": [{"type": "Ready", "status": "True"}]})
    elif spec.kind == "StorageClass":
        obj.setdefault("provisioner", "kubernetes.io/no-provisioner")
        obj.setdefault("reclaimPolicy", "Delete")
        obj.setdefault("volumeBindingMode", "Immediate")

    # Ensure metadata namespace/name exists.
    assign_metadata(obj, namespace=namespace, name=name)
    return obj


def handle_list_or_create(spec: ResourceSpec, *, namespace: Optional[str]):
    api_base = api_base_for(spec.api_version)
    plural = spec.plural

    def view():
        watch = request.args.get("watch", default="false").lower() in ("1", "true", "yes", "y", "on")
        resource_version = request.args.get("resourceVersion", "0")

        if watch:
            # Namespaced resources can be watched per-namespace; cluster-scoped resources ignore namespace.
            return Response(
                watch_stream(spec=spec, namespace=namespace, resource_version=resource_version),
                mimetype="text/event-stream",
            )

        if request.method == "GET":
            items = list_resources(spec, namespace=namespace)
            return jsonify(make_list(spec.kind, spec.api_version, items))

        # POST create
        data = request.get_json(force=True, silent=True) or {}

        meta = data.get("metadata") or {}
        name = meta.get("name")
        if not name:
            gen = meta.get("generateName")
            if gen:
                name = f"{gen}{uuid.uuid4().hex[:8]}"
            else:
                return jsonify({"message": "metadata.name or metadata.generateName is required"}), 400

        if spec.namespaced and namespace is None:
            # If you hit the cluster-wide list/create route, require namespace in body.
            namespace_from_body = meta.get("namespace")
            if not namespace_from_body:
                return jsonify({"message": "namespace is required for namespaced resources"}), 400
            ns = namespace_from_body
        else:
            ns = namespace

        obj = _default_object_for_create(spec, data, ns, name)
        store[spec.kind][_key_for(spec, namespace=ns, name=name)] = obj
        push_event("ADDED", obj)
        return jsonify(obj)

    view.__name__ = f"list_or_create_{spec.kind}_{'ns' if spec.namespaced else 'cluster'}"
    return view


def handle_detail_update_delete(spec: ResourceSpec, *, namespace: Optional[str]):
    def view(name: str):
        if request.method == "GET":
            obj = get_resource(spec, namespace=namespace, name=name)
            if not obj:
                return jsonify({}), 404
            return jsonify(obj)

        if request.method in ("PUT", "PATCH"):
            existing = get_resource(spec, namespace=namespace, name=name)
            if not existing:
                return jsonify({"message": "Not found"}), 404

            patch = request.get_json(force=True, silent=True) or {}
            # For PATCH/PUT in this mock, we do a merge. Real Kubernetes behavior differs.
            merged = deep_merge(existing, patch)
            merged = _ensure_correct_kind_fields(spec, merged)
            assign_metadata(merged, namespace=namespace, name=name)
            # Keep spec/status objects if provided.
            store[spec.kind][_key_for(spec, namespace=namespace, name=name)] = merged
            push_event("MODIFIED", merged)
            return jsonify(merged)

        if request.method == "DELETE":
            existing = get_resource(spec, namespace=namespace, name=name)
            if not existing:
                return jsonify({"status": "Success"})
            store[spec.kind].pop(_key_for(spec, namespace=namespace, name=name), None)
            # Use a tombstone object for the event.
            tomb = copy.deepcopy(existing)
            tomb.setdefault("metadata", {})
            tomb["metadata"]["resourceVersion"] = next_rv()
            push_event("DELETED", tomb, rv=tomb["metadata"]["resourceVersion"])
            return jsonify({"status": "Success"})

        return jsonify({"message": "Unsupported method"}), 405

    view.__name__ = f"detail_update_delete_{spec.kind}"
    return view


def handle_watch_subresource(spec: ResourceSpec, *, namespace: Optional[str]):
    def view():
        resource_version = request.args.get("resourceVersion", "0")
        return Response(
            watch_stream(spec=spec, namespace=namespace, resource_version=resource_version),
            mimetype="text/event-stream",
        )

    view.__name__ = f"watch_{spec.kind}"
    return view


def handle_cluster_list(spec: ResourceSpec):
    # For namespaced resources, kubectl often lists with -A by calling the cluster-scoped plural.
    def view():
        watch = request.args.get("watch", default="false").lower() in ("1", "true", "yes", "y", "on")
        resource_version = request.args.get("resourceVersion", "0")
        if watch:
            return Response(
                watch_stream(spec=spec, namespace=None, resource_version=resource_version),
                mimetype="text/event-stream",
            )
        items = list_resources(spec, namespace=None)
        return jsonify(make_list(spec.kind, spec.api_version, items))

    view.__name__ = f"cluster_list_{spec.kind}"
    return view


def register_routes() -> None:
    for kind, spec in RESOURCE_SPECS.items():
        base = api_base_for(spec.api_version)
        plural = spec.plural

        if spec.kind == "Namespace":
            # Namespace is core but cluster-scoped. Create/list at /api/v1/namespaces.
            app.add_url_rule(
                f"{base}/{plural}",
                view_func=handle_list_or_create(spec, namespace=None),
                methods=["GET", "POST"],
                endpoint=f"ns_list_create_{kind}",
            )
            app.add_url_rule(
                f"{base}/{plural}/<name>",
                view_func=handle_detail_update_delete(spec, namespace=None),
                methods=["GET", "PUT", "PATCH", "DELETE"],
                endpoint=f"ns_detail_{kind}",
            )
            app.add_url_rule(
                f"{base}/{plural}/watch",
                view_func=handle_watch_subresource(spec, namespace=None),
                methods=["GET"],
                endpoint=f"ns_watch_{kind}",
            )
            continue

        if spec.namespaced:
            # Namespaced list/create
            app.add_url_rule(
                f"{base}/namespaces/<ns>/{plural}",
                view_func=lambda ns, s=spec: handle_list_or_create(s, namespace=ns)(),  # type: ignore[misc]
                methods=["GET", "POST"],
                endpoint=f"{spec.kind}_list_create_ns",
            )
            # Namespaced detail
            app.add_url_rule(
                f"{base}/namespaces/<ns>/{plural}/<name>",
                view_func=lambda ns, name, s=spec: handle_detail_update_delete(s, namespace=ns)(name),  # type: ignore[misc]
                methods=["GET", "PUT", "PATCH", "DELETE"],
                endpoint=f"{spec.kind}_detail_ns",
            )
            # Namespaced watch subresource
            app.add_url_rule(
                f"{base}/namespaces/<ns>/{plural}/watch",
                view_func=lambda ns, s=spec: handle_watch_subresource(s, namespace=ns)(),  # type: ignore[misc]
                methods=["GET"],
                endpoint=f"{spec.kind}_watch_ns",
            )

            # Cluster-wide list for -A
            app.add_url_rule(
                f"{base}/{plural}",
                view_func=handle_cluster_list(spec),
                methods=["GET"],
                endpoint=f"{spec.kind}_list_all",
            )
        else:
            # Cluster-scoped list/create
            app.add_url_rule(
                f"{base}/{plural}",
                view_func=handle_list_or_create(spec, namespace=None),
                methods=["GET", "POST"],
                endpoint=f"{spec.kind}_list_create_cluster",
            )
            # Cluster-scoped detail
            app.add_url_rule(
                f"{base}/{plural}/<name>",
                view_func=handle_detail_update_delete(spec, namespace=None),
                methods=["GET", "PUT", "PATCH", "DELETE"],
                endpoint=f"{spec.kind}_detail_cluster",
            )
            # Cluster-scoped watch
            app.add_url_rule(
                f"{base}/{plural}/watch",
                view_func=handle_watch_subresource(spec, namespace=None),
                methods=["GET"],
                endpoint=f"{spec.kind}_watch_cluster",
            )


register_routes()


# -----------------------------
# Seed some objects so kubectl "feels" alive
# -----------------------------
def seed_defaults() -> None:
    for ns in ["default", "kube-system"]:
        ensure_namespace(ns)

    # Fake node(s)
    node_spec = RESOURCE_SPECS["Node"]
    if not get_resource(node_spec, namespace=None, name="fake-node"):
        obj: Dict[str, Any] = {"metadata": {"name": "fake-node"}}
        obj = _default_object_for_create(node_spec, obj, None, "fake-node")
        store["Node"][_key_for(node_spec, namespace=None, name="fake-node")] = obj
        push_event("ADDED", obj)

    # kube-system control-plane-ish pods
    for pod_name in ["kube-apiserver", "kube-controller-manager", "kube-scheduler", "coredns"]:
        pod_spec = RESOURCE_SPECS["Pod"]
        ns = "kube-system"
        if not get_resource(pod_spec, namespace=ns, name=pod_name):
            obj = {"metadata": {"name": pod_name}, "spec": {}}
            obj = _default_object_for_create(pod_spec, obj, ns, pod_name)
            store["Pod"][_key_for(pod_spec, namespace=ns, name=pod_name)] = obj
            push_event("ADDED", obj)

    # kube-dns service
    svc_spec = RESOURCE_SPECS["Service"]
    ns = "kube-system"
    if not get_resource(svc_spec, namespace=ns, name="kube-dns"):
        obj = {"metadata": {"name": "kube-dns"}, "spec": {"ports": [{"port": 53}]}}
        obj = _default_object_for_create(svc_spec, obj, ns, "kube-dns")
        store["Service"][_key_for(svc_spec, namespace=ns, name="kube-dns")] = obj
        push_event("ADDED", obj)

    # Minimal storageclass
    sc_spec = RESOURCE_SPECS["StorageClass"]
    if not get_resource(sc_spec, namespace=None, name="standard"):
        obj = {"metadata": {"name": "standard"}}
        obj = _default_object_for_create(sc_spec, obj, None, "standard")
        store["StorageClass"][_key_for(sc_spec, namespace=None, name="standard")] = obj
        push_event("ADDED", obj)

    # ClusterRole/Binding (for basic auth flows)
    cr_spec = RESOURCE_SPECS["ClusterRole"]
    if not get_resource(cr_spec, namespace=None, name="cluster-admin"):
        obj = {"metadata": {"name": "cluster-admin"}, "rules": []}
        obj = _default_object_for_create(cr_spec, obj, None, "cluster-admin")
        store["ClusterRole"][_key_for(cr_spec, namespace=None, name="cluster-admin")] = obj
        push_event("ADDED", obj)

    crb_spec = RESOURCE_SPECS["ClusterRoleBinding"]
    if not get_resource(crb_spec, namespace=None, name="cluster-admin-binding"):
        obj = {
            "metadata": {"name": "cluster-admin-binding"},
            "roleRef": {"kind": "ClusterRole", "name": "cluster-admin", "apiGroup": "rbac.authorization.k8s.io"},
            "subjects": [],
        }
        obj = _default_object_for_create(crb_spec, obj, None, "cluster-admin-binding")
        store["ClusterRoleBinding"][_key_for(crb_spec, namespace=None, name="cluster-admin-binding")] = obj
        push_event("ADDED", obj)

    # Default pod in default namespace
    pod_spec = RESOURCE_SPECS["Pod"]
    ns = "default"
    if not get_resource(pod_spec, namespace=ns, name="default-pod"):
        obj = {"metadata": {"name": "default-pod"}, "spec": {}}
        obj = _default_object_for_create(pod_spec, obj, ns, "default-pod")
        store["Pod"][_key_for(pod_spec, namespace=ns, name="default-pod")] = obj
        push_event("ADDED", obj)


seed_defaults()


def main() -> None:
    ensure_namespace("default")
    ensure_namespace("kube-system")
    app.run(host="0.0.0.0", port=8080, debug=True)


if __name__ == "__main__":
    main()
