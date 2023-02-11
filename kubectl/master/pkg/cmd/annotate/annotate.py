import ctypes
import json
import jsonpatch
import resource

class AnnotateFlags(ctypes.Structure):
    _fields_ = [
        ('all', ctypes.c_bool),
        ('allNamespaces', ctypes.c_bool),
        ('DryRunStrategy', ctypes.c_int),
        ('FieldManager', ctypes.c_char_p),
        ('FieldSelector', ctypes.c_char_p),
        ('FilenameOptions', ctypes.c_void_p),
        ('List', ctypes.c_bool),
        ('Local', ctypes.c_bool),
        ('OutputFormat', ctypes.c_char_p),
        ('overwrite', ctypes.c_bool),
        ('PrintFlags', ctypes.c_void_p),
        ('RecordFlags', ctypes.c_void_p),
        ('resourceVersion', ctypes.c_char_p),
        ('Selector', ctypes.c_char_p),
        ('IOStreams', ctypes.c_void_p)
    ]
    def __init__(self, streams):
        self.all = False
        self.allNamespaces = False
        self.dryRunStrategy = None # TODO: determine the equivalent in Python
        self.fieldManager = ""
        self.fieldSelector = ""
        self.filenameOptions = resource.FilenameOptions()
        self.list = False
        self.local = False
        self.outputFormat = ""
        self.overwrite = False
        self.printFlags = None # TODO: determine the equivalent in Python
        self.recordFlags = None # TODO: determine the equivalent in Python
        self.resourceVersion = ""
        self.selector = ""
        self.iostreams = streams

        
#
class AnnotateOptions(ctypes.Structure):
    _fields_ = [
        ('all', ctypes.c_bool),
        ('allNamespaces', ctypes.c_bool),
        ('builder', ctypes.c_void_p),
        ('dryRunStrategy', ctypes.c_int),
        ('enforceNamespace', ctypes.c_bool),
        ('fieldSelector', ctypes.c_char_p),
        ('fieldManager', ctypes.c_char_p),
        ('FilenameOptions', ctypes.c_void_p),
        ('IOStreams', ctypes.c_void_p),
        ('list', ctypes.c_bool),
        ('local', ctypes.c_bool),
        ('namespace', ctypes.c_char_p),
        ('newAnnotations', ctypes.c_void_p),
        ('overwrite', ctypes.c_bool),
        ('PrintObj', ctypes.c_void_p),
        ('Recorder', ctypes.c_void_p),
        ('resources', ctypes.c_void_p),
        ('resourceVersion', ctypes.c_char_p),
        ('removeAnnotations', ctypes.c_void_p),
        ('selector', ctypes.c_char_p),
        ('unstructuredClientForMapping', ctypes.c_void_p)
    ]
    def __init__(self):
        self.all = False
        self.allNamespaces = False
        self.builder = None # TODO: determine the equivalent in Python
        self.dryRunStrategy = None # TODO: determine the equivalent in Python
        self.enforceNamespace = False
        self.fieldSelector = ""
        self.fieldManager = ""
        self.filenameOptions = resource.FilenameOptions()
        self.iostreams = None # TODO: determine the equivalent in Python
        self.list = False
        self.local = False
        self.namespace = ""
        self.newAnnotations = {}
        self.overwrite = False
        self.printObj = None # TODO: determine the equivalent in Python
        self.recorder = None # TODO: determine the equivalent in Python
        self.resources = []
        self.resourceVersion = ""
        self.removeAnnotations = []
        self.selector = ""
        self.unstructuredClientForMapping = None # TODO: determine the equivalent in Python
