import ctypes
import threading
from typing import List, Optional, Dict
from dataclasses import dataclass, field
from datetime import datetime

class Fs(ctypes.Structure):
    _fields_ = [
        ("rw_lock", ctypes.py_object),
        ("root", ctypes.c_char_p)
    ]
    
    def __init__(self, root: str):
        super().__init__()
        self.rw_lock = threading.RLock()
        self.root = root.encode()

class RootFS(ctypes.Structure):
    _fields_ = [
        ("type", ctypes.c_char_p),
        ("diff_ids", ctypes.POINTER(ctypes.c_char_p))
    ]
    
    def __init__(self, type: str, diff_ids: List[str]):
        super().__init__()
        self.type = type.encode()
        self.diff_ids = (ctypes.c_char_p * len(diff_ids))(*[d.encode() for d in diff_ids])

class LocalImageCache(ctypes.Structure):
    _fields_ = [
        ("store", ctypes.py_object)
    ]
    
    def __init__(self, store):
        super().__init__()
        self.store = store

class ImageCache(ctypes.Structure):
    _fields_ = [
        ("sources", ctypes.POINTER(ctypes.py_object)),
        ("store", ctypes.py_object),
        ("local_image_cache", ctypes.POINTER(LocalImageCache))
    ]
    
    def __init__(self, sources: List, store, local_image_cache: Optional[LocalImageCache]):
        super().__init__()
        self.sources = (ctypes.py_object * len(sources))(*sources)
        self.store = store
        self.local_image_cache = ctypes.pointer(local_image_cache) if local_image_cache else None

class ImageDescriptor(ctypes.Structure):
    _fields_ = [
        ("refs", ctypes.POINTER(ctypes.c_char_p)),
        ("layers", ctypes.POINTER(ctypes.c_char_p)),
        ("image", ctypes.py_object),
        ("layer_ref", ctypes.py_object)
    ]
    
    def __init__(self, refs: List[str], layers: List[str], image, layer_ref):
        super().__init__()
        self.refs = (ctypes.c_char_p * len(refs))(*[r.encode() for r in refs])
        self.layers = (ctypes.c_char_p * len(layers))(*[l.encode() for l in layers])
        self.image = image
        self.layer_ref = layer_ref

class SaveSession(ctypes.Structure):
    _fields_ = [
        ("tarexporter", ctypes.py_object),
        ("out_dir", ctypes.c_char_p),
        ("images", ctypes.py_object),
        ("saved_layers", ctypes.py_object),
        ("saved_configs", ctypes.py_object)
    ]
    
    def __init__(self, tarexporter, out_dir: str, images: Dict, saved_layers: Dict, saved_configs: Dict):
        super().__init__()
        self.tarexporter = tarexporter
        self.out_dir = out_dir.encode()
        self.images = images
        self.saved_layers = saved_layers
        self.saved_configs = saved_configs

class ParentLink(ctypes.Structure):
    _fields_ = [
        ("id", ctypes.c_char_p),
        ("parent_id", ctypes.c_char_p)
    ]
    
    def __init__(self, id: str, parent_id: str):
        super().__init__()
        self.id = id.encode()
        self.parent_id = parent_id.encode()

class ManifestItem(ctypes.Structure):
    _fields_ = [
        ("config", ctypes.c_char_p),
        ("repo_tags", ctypes.POINTER(ctypes.c_char_p)),
        ("layers", ctypes.POINTER(ctypes.c_char_p)),
        ("parent", ctypes.c_char_p),
        ("layer_sources", ctypes.py_object)
    ]
    
    def __init__(self, config: str, repo_tags: List[str], layers: List[str], parent: Optional[str] = None, layer_sources: Optional[Dict] = None):
        super().__init__()
        self.config = config.encode()
        self.repo_tags = (ctypes.c_char_p * len(repo_tags))(*[rt.encode() for rt in repo_tags])
        self.layers = (ctypes.c_char_p * len(layers))(*[l.encode() for l in layers])
        self.parent = parent.encode() if parent else None
        self.layer_sources = layer_sources

class TarExporter(ctypes.Structure):
    _fields_ = [
        ("is_store", ctypes.py_object),
        ("lss", ctypes.py_object),
        ("rs", ctypes.py_object),
        ("logger_img_event", ctypes.py_object)
    ]
    
    def __init__(self, is_store, lss, rs, logger_img_event):
        super().__init__()
        self.is_store = is_store
        self.lss = lss
        self.rs = rs
        self.logger_img_event = logger_img_event

class TCase(ctypes.Structure):
    _fields_ = [
        ("input", ctypes.POINTER(ctypes.c_ubyte)),
        ("expected", ctypes.c_char_p)
    ]
    
    def __init__(self, input: bytes, expected: str):
        super().__init__()
        self.input = (ctypes.c_ubyte * len(input))(*input)
        self.expected = expected.encode()

class V1Image(ctypes.Structure):
    _fields_ = [
        ("id", ctypes.c_char_p),
        ("parent", ctypes.c_char_p),
        ("comment", ctypes.c_char_p),
        ("created", ctypes.POINTER(ctypes.c_long)),
        ("container", ctypes.c_char_p),
        ("container_config", ctypes.py_object),
        ("docker_version", ctypes.c_char_p),
        ("author", ctypes.c_char_p),
        ("config", ctypes.py_object),
        ("architecture", ctypes.c_char_p),
        ("variant", ctypes.c_char_p),
        ("os", ctypes.c_char_p),
        ("size", ctypes.c_long)
    ]
    
    def __init__(self, id: Optional[str] = None, parent: Optional[str] = None, comment: Optional[str] = None, created: Optional[datetime] = None, container: Optional[str] = None, container_config = None, docker_version: Optional[str] = None, author: Optional[str] = None, config = None, architecture: Optional[str] = None, variant: Optional[str] = None, os: Optional[str] = None, size: Optional[int] = None):
        super().__init__()
        self.id = id.encode() if id else None
        self.parent = parent.encode() if parent else None
        self.comment = comment.encode() if comment else None
        self.created = ctypes.pointer(ctypes.c_long(int(created.timestamp()))) if created else None
        self.container = container.encode() if container else None
        self.container_config = container_config
        self.docker_version = docker_version.encode() if docker_version else None
        self.author = author.encode() if author else None
        self.config = config
        self.architecture = architecture.encode() if architecture else None
        self.variant = variant.encode() if variant else None
        self.os = os.encode() if os else None
        self.size = size

class Image(ctypes.Structure):
    _fields_ = [
        ("v1_image", V1Image),
        ("parent", ctypes.c_char_p),
        ("root_fs", ctypes.py_object),
        ("history", ctypes.POINTER(ctypes.py_object)),
        ("os_version", ctypes.c_char_p),
        ("os_features", ctypes.POINTER(ctypes.c_char_p)),
        ("raw_json", ctypes.POINTER(ctypes.c_ubyte)),
        ("computed_id", ctypes.c_char_p),
        ("details", ctypes.py_object)
    ]
    
    def __init__(self, v1_image: V1Image, parent: Optional[str] = None, root_fs = None, history: List = [], os_version: Optional[str] = None, os_features: List[str] = [], raw_json: Optional[bytes] = None, computed_id: Optional[str] = None, details = None):
        super().__init__()
        self.v1_image = v1_image
        self.parent = parent.encode() if parent else None
        self.root_fs = root_fs
        self.history = (ctypes.py_object * len(history))(*history)
        self.os_version = os_version.encode() if os_version else None
        self.os_features = (ctypes.c_char_p * len(os_features))(*[f.encode() for f in os_features])
        self.raw_json = (ctypes.c_ubyte * len(raw_json))(*raw_json) if raw_json else None
        self.computed_id = computed_id.encode() if computed_id else None
        self.details = details

class Details(ctypes.Structure):
    _fields_ = [
        ("references", ctypes.POINTER(ctypes.c_char_p)),
        ("size", ctypes.c_long),
        ("metadata", ctypes.py_object),
        ("driver", ctypes.c_char_p),
        ("last_updated", ctypes.c_long)
    ]
    
    def __init__(self, references: List[str], size: int, metadata: Dict[str, str], driver: str, last_updated: datetime):
        super().__init__()
        self.references = (ctypes.c_char_p * len(references))(*[r.encode() for r in references])
        self.size = size
        self.metadata = metadata
        self.driver = driver.encode()
        self.last_updated = int(last_updated.timestamp())

class ChildConfig(ctypes.Structure):
    _fields_ = [
        ("container_id", ctypes.c_char_p),
        ("author", ctypes.c_char_p),
        ("comment", ctypes.c_char_p),
        ("diff_id", ctypes.c_char_p),
        ("container_config", ctypes.py_object),
        ("config", ctypes.py_object)
    ]
    
    def __init__(self, container_id: str, author: str, comment: str, diff_id: str, container_config, config):
        super().__init__()
        self.container_id = container_id.encode()
        self.author = author.encode()
        self.comment = comment.encode()
        self.diff_id = diff_id.encode()
        self.container_config = container_config
        self.config = config

class ImageMeta(ctypes.Structure):
    _fields_ = [
        ("layer", ctypes.py_object),
        ("children", ctypes.py_object)
    ]
    
    def __init__(self, layer, children: Dict):
        super().__init__()
        self.layer = layer
        self.children = children

class Store(ctypes.Structure):
    _fields_ = [
        ("rw_lock", ctypes.py_object),
        ("lss", ctypes.py_object),
        ("images", ctypes.py_object),
        ("fs", ctypes.py_object),
        ("digest_set", ctypes.py_object)
    ]
    
    def __init__(self, lss, images: Dict, fs, digest_set):
        super().__init__()
        self.rw_lock = threading.RLock()
        self.lss = lss
        self.images = images
        self.fs = fs
        self.digest_set = digest_set

class UnpackSizeCounter(ctypes.Structure):
    _fields_ = [
        ("unpacker", ctypes.py_object),
        ("size", ctypes.c_long)
    ]
    
    def __init__(self, unpacker, size: int):
        super().__init__()
        self.unpacker = unpacker
        self.size = size

class PackSizeCounter(ctypes.Structure):
    _fields_ = [
        ("packer", ctypes.py_object),
        ("size", ctypes.c_long)
    ]
    
    def __init__(self, packer, size: int):
        super().__init__()
        self.packer = packer
        self.size = size

class Metadata(ctypes.Structure):
    _fields_ = [
        ("chain_id", ctypes.c_char_p),
        ("diff_id", ctypes.c_char_p),
        ("size", ctypes.c_long),
        ("diff_size", ctypes.c_long)
    ]
    
    def __init__(self, chain_id: str, diff_id: str, size: int, diff_size: int):
        super().__init__()
        self.chain_id = chain_id.encode()
        self.diff_id = diff_id.encode()
        self.size = size
        self.diff_size = diff_size

class CreateRWLayerOpts(ctypes.Structure):
    _fields_ = [
        ("mount_label", ctypes.c_char_p),
        ("init_func", ctypes.py_object),
        ("storage_opt", ctypes.py_object)
    ]
    
    def __init__(self, mount_label: str, init_func, storage_opt: Dict):
        super().__init__()
        self.mount_label = mount_label.encode()
        self.init_func = init_func
        self.storage_opt = storage_opt

class ChangeSorter(ctypes.Structure):
    _fields_ = [
        ("changes", ctypes.py_object)
    ]
    
    def __init__(self, changes: List):
        super().__init__()
        self.changes = changes

class LayerStore(ctypes.Structure):
    _fields_ = [
        ("store", ctypes.py_object),
        ("driver", ctypes.py_object),
        ("use_tar_split", ctypes.c_bool),
        ("layer_map", ctypes.py_object),
        ("layer_lock", ctypes.py_object),
        ("mounts", ctypes.py_object),
        ("mount_lock", ctypes.py_object),
        ("locker", ctypes.py_object)
    ]
    
    def __init__(self, store, driver, use_tar_split: bool, layer_map: Dict, mounts: Dict, locker):
        super().__init__()
        self.store = store
        self.driver = driver
        self.use_tar_split = use_tar_split
        self.layer_map = layer_map
        self.layer_lock = threading.Lock()
        self.mounts = mounts
        self.mount_lock = threading.Lock()
        self.locker = locker

class StoreOptions(ctypes.Structure):
    _fields_ = [
        ("root", ctypes.c_char_p),
        ("metadata_store_path_template", ctypes.c_char_p),
        ("graph_driver", ctypes.c_char_p),
        ("graph_driver_options", ctypes.POINTER(ctypes.c_char_p)),
        ("id_mapping", ctypes.py_object),
        ("plugin_getter", ctypes.py_object),
        ("experimental_enabled", ctypes.c_bool)
    ]
    
    def __init__(self, root: str, metadata_store_path_template: str, graph_driver: str, graph_driver_options: List[str], id_mapping, plugin_getter, experimental_enabled: bool):
        super().__init__()
        self.root = root.encode()
        self.metadata_store_path_template = metadata_store_path_template.encode()
        self.graph_driver = graph_driver.encode()
        self.graph_driver_options = (ctypes.c_char_p * len(graph_driver_options))(*[opt.encode() for opt in graph_driver_options])
        self.id_mapping = id_mapping
        self.plugin_getter = plugin_getter
        self.experimental_enabled = experimental_enabled

class NaiveDiffPathDriver(ctypes.Structure):
    _fields_ = [
        ("driver", ctypes.py_object)
    ]
    
    def __init__(self, driver):
        super().__init__()
        self.driver = driver

class FileGetPutter(ctypes.Structure):
    _fields_ = [
        ("storage_file_getter", ctypes.py_object),
        ("driver", ctypes.py_object),
        ("id", ctypes.c_char_p)
    ]
    
    def __init__(self, storage_file_getter, driver, id: str):
        super().__init__()
        self.storage_file_getter = storage_file_getter
        self.driver = driver
        self.id = id.encode()

class TestFile(ctypes.Structure):
    _fields_ = [
        ("name", ctypes.c_char_p),
        ("content", ctypes.POINTER(ctypes.c_ubyte)),
        ("permission", ctypes.c_uint)
    ]
    
    def __init__(self, name: str, content: bytes, permission: int):
        super().__init__()
        self.name = name.encode()
        self.content = (ctypes.c_ubyte * len(content))(*content)
        self.permission = permission

class FileMetadataStore(ctypes.Structure):
    _fields_ = [
        ("root", ctypes.c_char_p)
    ]
    
    def __init__(self, root: str):
        super().__init__()
        self.root = root.encode()

class FileMetadataTransaction(ctypes.Structure):
    _fields_ = [
        ("store", ctypes.py_object),
        ("ws", ctypes.py_object)
    ]
    
    def __init__(self, store: FileMetadataStore, ws):
        super().__init__()
        self.store = store
        self.ws = ws

class MountedLayer(ctypes.Structure):
    _fields_ = [
        ("name", ctypes.c_char_p),
        ("mount_id", ctypes.c_char_p),
        ("init_id", ctypes.c_char_p),
        ("parent", ctypes.py_object),
        ("layer_store", ctypes.py_object),
        ("lock", ctypes.py_object),
        ("references", ctypes.py_object)
    ]
    
    def __init__(self, name: str, mount_id: str, init_id: str, parent, layer_store, references: Dict):
        super().__init__()
        self.name = name.encode()
        self.mount_id = mount_id.encode()
        self.init_id = init_id.encode()
        self.parent = parent
        self.layer_store = layer_store
        self.lock = threading.Lock()
        self.references = references

class ReferencedRWLayer(ctypes.Structure):
    _fields_ = [
        ("mounted_layer", ctypes.py_object)
    ]
    
    def __init__(self, mounted_layer: MountedLayer):
        super().__init__()
        self.mounted_layer = mounted_layer

class RoLayer(ctypes.Structure):
    _fields_ = [
        ("chain_id", ctypes.c_char_p),
        ("diff_id", ctypes.c_char_p),
        ("parent", ctypes.py_object),
        ("cache_id", ctypes.c_char_p),
        ("size", ctypes.c_long),
        ("layer_store", ctypes.py_object),
        ("descriptor", ctypes.py_object),
        ("reference_count", ctypes.c_int),
        ("references", ctypes.py_object)
    ]
    
    def __init__(self, chain_id: str, diff_id: str, parent, cache_id: str, size: int, layer_store, descriptor, reference_count: int, references: Dict):
        super().__init__()
        self.chain_id = chain_id.encode()
        self.diff_id = diff_id.encode()
        self.parent = parent
        self.cache_id = cache_id.encode()
        self.size = size
        self.layer_store = layer_store
        self.descriptor = descriptor
        self.reference_count = reference_count
        self.references = references

class ReferencedCacheLayer(ctypes.Structure):
    _fields_ = [
        ("ro_layer", ctypes.py_object)
    ]
    
    def __init__(self, ro_layer: RoLayer):
        super().__init__()
        self.ro_layer = ro_layer

class VerifiedReadCloser(ctypes.Structure):
    _fields_ = [
        ("rc", ctypes.py_object),
        ("digest", ctypes.c_char_p),
        ("verifier", ctypes.py_object)
    ]
    
    def __init__(self, rc, digest: str, verifier):
        super().__init__()
        self.rc = rc
        self.digest = digest.encode()
        self.verifier = verifier

class ImageSaveManifestEntry(ctypes.Structure):
    _fields_ = [
        ("config", ctypes.c_char_p),
        ("repo_tags", ctypes.POINTER(ctypes.c_char_p)),
        ("layers", ctypes.POINTER(ctypes.c_char_p))
    ]
    
    def __init__(self, config: str, repo_tags: List[str], layers: List[str]):
        super().__init__()
        self.config = config.encode()
        self.repo_tags = (ctypes.c_char_p * len(repo_tags))(*[rt.encode() for rt in repo_tags])
        self.layers = (ctypes.c_char_p * len(layers))(*[l.encode() for l in layers])

class TestCase(ctypes.Structure):
    _fields_ = [
        ("image", ctypes.c_char_p),
        ("expected_oci_ref", ctypes.c_char_p),
        ("expected_containerd_ref", ctypes.c_char_p)
    ]
    
    def __init__(self, image: str, expected_oci_ref: str, expected_containerd_ref: str):
        super().__init__()
        self.image = image.encode()
        self.expected_oci_ref = expected_oci_ref.encode()
        self.expected_containerd_ref = expected_containerd_ref.encode()

class BuildLine(ctypes.Structure):
    _fields_ = [
        ("stream", ctypes.c_char_p),
        ("aux", ctypes.py_object)
    ]
    
    def __init__(self, stream: str, aux: Dict[str, str]):
        super().__init__()
        self.stream = stream.encode()
        self.aux = aux

class TestWriter(ctypes.Structure):
    _fields_ = [
        ("test", ctypes.py_object)
    ]
    
    def __init__(self, test):
        super().__init__()
        self.test = test

class ExecResult(ctypes.Structure):
    _fields_ = [
        ("exit_code", ctypes.c_int),
        ("out_buffer", ctypes.py_object),
        ("err_buffer", ctypes.py_object)
    ]
    
    def __init__(self, exit_code: int, out_buffer, err_buffer):
        super().__init__()
        self.exit_code = exit_code
        self.out_buffer = out_buffer
        self.err_buffer = err_buffer

class TestContainerConfig(ctypes.Structure):
    _fields_ = [
        ("name", ctypes.c_char_p),
        ("config", ctypes.py_object),
        ("host_config", ctypes.py_object),
        ("networking_config", ctypes.py_object),
        ("platform", ctypes.py_object)
    ]
    
    def __init__(self, name: str, config, host_config, networking_config, platform):
        super().__init__()
        self.name = name.encode()
        self.config = config
        self.host_config = host_config
        self.networking_config = networking_config
        self.platform = platform

class RunResult(ctypes.Structure):
    _fields_ = [
        ("container_id", ctypes.c_char_p),
        ("exit_code", ctypes.c_int),
        ("stdout", ctypes.py_object),
        ("stderr", ctypes.py_object)
    ]
    
    def __init__(self, container_id: str, exit_code: int, stdout, stderr):
        super().__init__()
        self.container_id = container_id.encode()
        self.exit_code = exit_code
        self.stdout = stdout
        self.stderr = stderr

class Streams(ctypes.Structure):
    _fields_ = [
        ("stdout", ctypes.py_object),
        ("stderr", ctypes.py_object)
    ]
    
    def __init__(self, stdout, stderr):
        super().__init__()
        self.stdout = stdout
        self.stderr = stderr

class ContainerOutput(ctypes.Structure):
    _fields_ = [
        ("stdout", ctypes.c_char_p),
        ("stderr", ctypes.c_char_p)
    ]
    
    def __init__(self, stdout: str, stderr: str):
        super().__init__()
        self.stdout = stdout.encode()
        self.stderr = stderr.encode()

class StringHandler(ctypes.Structure):
    _fields_ = [
        ("ansi_event_handler", ctypes.py_object),
        ("cursor", ctypes.c_int),
        ("buffer", ctypes.POINTER(ctypes.c_ubyte))
    ]
    
    def __init__(self, ansi_event_handler, cursor: int, buffer: bytes):
        super().__init__()
        self.ansi_event_handler = ansi_event_handler
        self.cursor = cursor
        self.buffer = (ctypes.c_ubyte * len(buffer))(*buffer)

class BridgesOpts(ctypes.Structure):
    _fields_ = [
        ("bridge1_opts", ctypes.POINTER(ctypes.py_object)),
        ("bridge2_opts", ctypes.POINTER(ctypes.py_object))
    ]
    
    def __init__(self, bridge1_opts: List[callable], bridge2_opts: List[callable]):
        super().__init__()
        self.bridge1_opts = (ctypes.py_object * len(bridge1_opts))(*bridge1_opts)
        self.bridge2_opts = (ctypes.py_object * len(bridge2_opts))(*bridge2_opts)

class TestStep(ctypes.Structure):
    _fields_ = [
        ("step_name", ctypes.c_char_p),
        ("fixed_cidr_v6", ctypes.c_char_p),
        ("expected_addrs", ctypes.POINTER(ctypes.c_char_p))
    ]
    
    def __init__(self, step_name: str, fixed_cidr_v6: str, expected_addrs: List[str]):
        super().__init__()
        self.step_name = step_name.encode()
        self.fixed_cidr_v6 = fixed_cidr_v6.encode()
        self.expected_addrs = (ctypes.c_char_p * len(expected_addrs))(*[addr.encode() for addr in expected_addrs])

class RestartTestCase(ctypes.Structure):
    _fields_ = [
        ("desc", ctypes.c_char_p),
        ("restart_policy", ctypes.py_object),
        ("x_running", ctypes.c_bool),
        ("x_running_live_restore", ctypes.c_bool),
        ("x_start", ctypes.c_bool),
        ("x_health_check", ctypes.c_bool)
    ]
    
    def __init__(self, desc: str, restart_policy, x_running: bool, x_running_live_restore: bool, x_start: bool, x_health_check: bool):
        super().__init__()
        self.desc = desc.encode()
        self.restart_policy = restart_policy
        self.x_running = x_running
        self.x_running_live_restore = x_running_live_restore
        self.x_start = x_start
        self.x_health_check = x_health_check

class Config(ctypes.Structure):
    _fields_ = [
        ("runtimes", ctypes.py_object)
    ]
    
    def __init__(self, runtimes: Dict[str, 'Runtime']):
        super().__init__()
        self.runtimes = runtimes

class AuthorizationController(ctypes.Structure):
    _fields_ = [
        ("req_res", ctypes.py_object),
        ("res_res", ctypes.py_object),
        ("version_req_count", ctypes.c_int),
        ("version_res_count", ctypes.c_int),
        ("requests_uris", ctypes.POINTER(ctypes.c_char_p)),
        ("req_user", ctypes.c_char_p),
        ("res_user", ctypes.c_char_p)
    ]
    
    def __init__(self, req_res, res_res, version_req_count: int, version_res_count: int, requests_uris: List[str], req_user: str, res_user: str):
        super().__init__()
        self.req_res = req_res
        self.res_res = res_res
        self.version_req_count = version_req_count
        self.version_res_count = version_res_count
        self.requests_uris = (ctypes.c_char_p * len(requests_uris))(*[uri.encode() for uri in requests_uris])
        self.req_user = req_user.encode()
        self.res_user = res_user.encode()

class StartLoggingRequest(ctypes.Structure):
    _fields_ = [
        ("file", ctypes.c_char_p)
    ]
    
    def __init__(self, file: str):
        super().__init__()
        self.file = file.encode()

class CapabilitiesResponse(ctypes.Structure):
    _fields_ = [
        ("cap", ctypes.py_object)
    ]
    
    def __init__(self, cap: Dict[str, bool]):
        super().__init__()
        self.cap = cap

class Driver(ctypes.Structure):
    _fields_ = [
        ("lock", ctypes.py_object),
        ("logs", ctypes.py_object)
    ]
    
    def __init__(self, logs: Dict[str, 'IoCloser']):
        super().__init__()
        self.lock = threading.Lock()
        self.logs = logs

class StopLoggingRequest(ctypes.Structure):
    _fields_ = [
        ("file", ctypes.c_char_p)
    ]
    
    def __init__(self, file: str):
        super().__init__()
        self.file = file.encode()

class Response(ctypes.Structure):
    _fields_ = [
        ("err", ctypes.c_char_p)
    ]
    
    def __init__(self, err: str):
        super().__init__()
        self.err = err.encode()

class Start(ctypes.Structure):
    _fields_ = [
        ("file", ctypes.c_char_p)
    ]
    
    def __init__(self, file: str):
        super().__init__()
        self.file = file.encode()

class GraphEventsCounter(ctypes.Structure):
    _fields_ = [
        ("activations", ctypes.c_int),
        ("creations", ctypes.c_int),
        ("removals", ctypes.c_int),
        ("gets", ctypes.c_int),
        ("puts", ctypes.c_int),
        ("stats", ctypes.c_int),
        ("cleanups", ctypes.c_int),
        ("exists", ctypes.c_int),
        ("init", ctypes.c_int),
        ("metadata", ctypes.c_int),
        ("diff", ctypes.c_int),
        ("apply_diff", ctypes.c_int),
        ("changes", ctypes.c_int),
        ("diff_size", ctypes.c_int)
    ]
    
    def __init__(self, activations: int, creations: int, removals: int, gets: int, puts: int, stats: int, cleanups: int, exists: int, init: int, metadata: int, diff: int, apply_diff: int, changes: int, diff_size: int):
        super().__init__()
        self.activations = activations
        self.creations = creations
        self.removals = removals
        self.gets = gets
        self.puts = puts
        self.stats = stats
        self.cleanups = cleanups
        self.exists = exists
        self.init = init
        self.metadata = metadata
        self.diff = diff
        self.apply_diff = apply_diff
        self.changes = changes
        self.diff_size = diff_size

class GraphDriverRequest(ctypes.Structure):
    _fields_ = [
        ("id", ctypes.c_char_p),
        ("parent", ctypes.c_char_p),
        ("mount_label", ctypes.c_char_p),
        ("read_only", ctypes.c_bool)
    ]
    
    def __init__(self, id: Optional[str] = None, parent: Optional[str] = None, mount_label: Optional[str] = None, read_only: Optional[bool] = None):
        super().__init__()
        self.id = id.encode() if id else None
        self.parent = parent.encode() if parent else None
        self.mount_label = mount_label.encode() if mount_label else None
        self.read_only = read_only

class GraphDriverResponse(ctypes.Structure):
    _fields_ = [
        ("err", ctypes.c_char_p),
        ("dir", ctypes.c_char_p),
        ("exists", ctypes.c_bool),
        ("status", ctypes.POINTER(ctypes.py_object)),
        ("metadata", ctypes.py_object),
        ("changes", ctypes.POINTER(ctypes.py_object)),
        ("size", ctypes.c_long)
    ]
    
    def __init__(self, err: Optional[str] = None, dir: Optional[str] = None, exists: Optional[bool] = None, status: Optional[List[Dict[str, str]]] = None, metadata: Optional[Dict[str, str]] = None, changes: Optional[List] = None, size: Optional[int] = None):
        super().__init__()
        self.err = err.encode() if err else None
        self.dir = dir.encode() if dir else None
        self.exists = exists
        self.status = (ctypes.py_object * len(status))(*status) if status else None
        self.metadata = metadata
        self.changes = (ctypes.py_object * len(changes))(*changes) if changes else None
        self.size = size

class OSContext(ctypes.Structure):
    _fields_ = [
        ("orig_ns", ctypes.py_object),
        ("new_ns", ctypes.py_object),
        ("tid", ctypes.c_int),
        ("caller", ctypes.c_char_p)
    ]
    
    def __init__(self, orig_ns, new_ns, tid: int, caller: str):
        super().__init__()
        self.orig_ns = orig_ns
        self.new_ns = new_ns
        self.tid = tid
        self.caller = caller.encode()

class SingleFileLayer(ctypes.Structure):
    _fields_ = [
        ("name", ctypes.c_char_p),
        ("content", ctypes.POINTER(ctypes.c_ubyte))
    ]
    
    def __init__(self, name: str, content: bytes):
        super().__init__()
        self.name = name.encode()
        self.content = (ctypes.c_ubyte * len(content))(*content)

class FileInLayer(ctypes.Structure):
    _fields_ = [
        ("path", ctypes.c_char_p),
        ("content", ctypes.POINTER(ctypes.c_ubyte))
    ]
    
    def __init__(self, path: str, content: bytes):
        super().__init__()
        self.path = path.encode()
        self.content = (ctypes.c_ubyte * len(content))(*content)

class JoinError(ctypes.Structure):
    _fields_ = [
        ("errs", ctypes.POINTER(ctypes.py_object))
    ]
    
    def __init__(self, errs: List[Exception]):
        super().__init__()
        self.errs = (ctypes.py_object * len(errs))(*errs)

class ErrNotAccessible(ctypes.Structure):
    _fields_ = [
        ("path", ctypes.c_char_p),
        ("cause", ctypes.py_object)
    ]
    
    def __init__(self, path: str, cause: Exception):
        super().__init__()
        self.path = path.encode()
        self.cause = cause

class ErrEscapesBase(ctypes.Structure):
    _fields_ = [
        ("base", ctypes.c_char_p),
        ("subpath", ctypes.c_char_p)
    ]
    
    def __init__(self, base: str, subpath: str):
        super().__init__()
        self.base = base.encode()
        self.subpath = subpath.encode()

class SafePath(ctypes.Structure):
    _fields_ = [
        ("path", ctypes.c_char_p),
        ("cleanup", ctypes.py_object),
        ("lock", ctypes.py_object),
        ("source_base", ctypes.c_char_p),
        ("source_subpath", ctypes.c_char_p)
    ]
    
    def __init__(self, path: str, cleanup: callable, source_base: str, source_subpath: str):
        super().__init__()
        self.path = path.encode()
        self.cleanup = cleanup
        self.lock = threading.Lock()
        self.source_base = source_base.encode()
        self.source_subpath = source_subpath.encode()

class EnvironCarrier(ctypes.Structure):
    _fields_ = [
        ("trace_parent", ctypes.c_char_p),
        ("trace_state", ctypes.c_char_p)
    ]
    
    def __init__(self, trace_parent: str, trace_state: str):
        super().__init__()
        self.trace_parent = trace_parent.encode()
        self.trace_state = trace_state.encode()

class ReaderCtx(ctypes.Structure):
    _fields_ = [
        ("ctx", ctypes.py_object),
        ("r", ctypes.py_object)
    ]
    
    def __init__(self, ctx, r):
        super().__init__()
        self.ctx = ctx
        self.r = r

class Composite(ctypes.Structure):
    _fields_ = [
        ("cleanups", ctypes.POINTER(ctypes.py_object))
    ]
    
    def __init__(self, cleanups: List[callable]):
        super().__init__()
        self.cleanups = (ctypes.py_object * len(cleanups))(*cleanups)

class ServerResponse(ctypes.Structure):
    _fields_ = [
        ("body", ctypes.py_object),
        ("header", ctypes.py_object),
        ("status_code", ctypes.c_int),
        ("req_url", ctypes.py_object)
    ]
    
    def __init__(self, body, header: Dict[str, str], status_code: int, req_url):
        super().__init__()
        self.body = body
        self.header = header
        self.status_code = status_code
        self.req_url = req_url

class HijackedConn(ctypes.Structure):
    _fields_ = [
        ("conn", ctypes.py_object),
        ("r", ctypes.py_object)
    ]
    
    def __init__(self, conn, r):
        super().__init__()
        self.conn = conn
        self.r = r

class HijackedConnCloseWriter(ctypes.Structure):
    _fields_ = [
        ("hijacked_conn", ctypes.py_object)
    ]
    
    def __init__(self, hijacked_conn: HijackedConn):
        super().__init__()
        self.hijacked_conn = hijacked_conn

class BytesBufferClose(ctypes.Structure):
    _fields_ = [
        ("buffer", ctypes.py_object)
    ]
    
    def __init__(self, buffer):
        super().__init__()
        self.buffer = buffer

class ConfigWrapper(ctypes.Structure):
    _fields_ = [
        ("config", ctypes.py_object),
        ("host_config", ctypes.py_object),
        ("networking_config", ctypes.py_object)
    ]
    
    def __init__(self, config, host_config, networking_config):
        super().__init__()
        self.config = config
        self.host_config = host_config
        self.networking_config = networking_config

class ErrConnectionFailed(ctypes.Structure):
    _fields_ = [
        ("error", ctypes.py_object)
    ]
    
    def __init__(self, error: Exception):
        super().__init__()
        self.error = error

class ObjectNotFoundError(ctypes.Structure):
    _fields_ = [
        ("object", ctypes.c_char_p),
        ("id", ctypes.c_char_p)
    ]
    
    def __init__(self, object: str, id: str):
        super().__init__()
        self.object = object.encode()
        self.id = id.encode()

class Client(ctypes.Structure):
    _fields_ = [
        ("scheme", ctypes.c_char_p),
        ("host", ctypes.c_char_p),
        ("proto", ctypes.c_char_p),
        ("addr", ctypes.c_char_p),
        ("base_path", ctypes.c_char_p),
        ("client", ctypes.py_object),
        ("version", ctypes.c_char_p),
        ("user_agent", ctypes.c_char_p),
        ("custom_http_headers", ctypes.py_object),
        ("manual_override", ctypes.c_bool),
        ("negotiate_version", ctypes.c_bool),
        ("negotiated", ctypes.c_bool),
        ("tp", ctypes.py_object),
        ("base_transport", ctypes.py_object)
    ]
    
    def __init__(self, scheme: str, host: str, proto: str, addr: str, base_path: str, client, version: str, user_agent: Optional[str], custom_http_headers: Dict[str, str], manual_override: bool, negotiate_version: bool, negotiated: bool, tp, base_transport):
        super().__init__()
        self.scheme = scheme.encode()
        self.host = host.encode()
        self.proto = proto.encode()
        self.addr = addr.encode()
        self.base_path = base_path.encode()
        self.client = client
        self.version = version.encode()
        self.user_agent = user_agent.encode() if user_agent else None
        self.custom_http_headers = custom_http_headers
        self.manual_override = manual_override
        self.negotiate_version = negotiate_version
        self.negotiated = negotiated
        self.tp = tp
        self.base_transport = base_transport
