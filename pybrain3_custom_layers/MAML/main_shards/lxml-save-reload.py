# Dylan Kenneth Eliot

"""
using lxml, create a reloadable saved map.

Then use it after reload. 

"""

from lxml import etree
from pybrain3.tools.shortcuts import buildNetwork
from pybrain3.structure import FeedForwardNetwork
import pybrain3.structure.modules as modules_pkg
import pybrain3.structure.connections as connections_pkg

# 1. Setup a network with mixed layer types for testing
net = buildNetwork(2, 3, 1) # hiddenclass=None lets us manually add types if needed

def save_network_strict_lxml(net, filename):
    # Root element defining the overall network architecture type
    root = etree.Element("Network", 
                         name="StrictPyBrainNet", 
                         py_class=type(net).__name__)

    # 2. Serialize Modules with Class Types
    modules_xml = etree.SubElement(root, "Modules")
    for mod in net.modules:
        m_elem = etree.SubElement(modules_xml, "Module")
        m_elem.set("name", mod.name)
        m_elem.set("size", str(mod.indim))
        # Capture the specific class (e.g., SigmoidLayer, LinearLayer)
        m_elem.set("class", type(mod).__name__)
        
        if hasattr(mod, 'params') and mod.params is not None:
            # Join weights with high precision
            m_elem.text = ",".join(map(str, mod.params))

    # 3. Serialize Connections with Class Types
    connections_xml = etree.SubElement(root, "Connections")
    for conn_key in net.connections:
        for c in net.connections[conn_key]:
            c_elem = etree.SubElement(connections_xml, "Connection")
            c_elem.set("class", type(c).__name__)
            c_elem.set("from", c.inmod.name)
            c_elem.set("to", c.outmod.name)
            
            if hasattr(c, 'params') and c.params is not None:
                c_elem.text = ",".join(map(str, c.params))

    # 4. Write to file with formatting
    tree = etree.ElementTree(root)
    tree.write(filename, pretty_print=True, xml_declaration=True, encoding="utf-8")

def load_network_strict_lxml(filename):
    tree = etree.parse(filename)
    root = tree.getroot()

    net = eval(root.get('py_class'))()

    module_objs = {}
    modules_pending_params = []
    connections_pending_params = []

    # --- First pass: collect connection topology ---
    from_names = set()
    to_names = set()

    for c_elem in root.find("Connections"):
        from_names.add(c_elem.get("from"))
        to_names.add(c_elem.get("to"))

    # --- Recreate modules ---
    for m_elem in root.find("Modules"):
        class_name = m_elem.get("class")
        name = m_elem.get("name")
        size = int(m_elem.get("size"))

        mod_class = getattr(modules_pkg, class_name)
        module = mod_class(size)
        module.name = name

        module_objs[name] = module
        modules_pending_params.append((module, m_elem.text))

    # --- Register modules correctly based on topology ---
    for name, module in module_objs.items():
        if name not in to_names:
            net.addInputModule(module)
        elif name not in from_names:
            net.addOutputModule(module)
        else:
            net.addModule(module)

    # --- Recreate connections ---
    for c_elem in root.find("Connections"):
        class_name = c_elem.get("class")
        from_name = c_elem.get("from")
        to_name = c_elem.get("to")

        inmod = module_objs[from_name]
        outmod = module_objs[to_name]

        conn_class = getattr(connections_pkg, class_name)
        connection = conn_class(inmod, outmod)

        net.addConnection(connection)
        connections_pending_params.append((connection, c_elem.text))

    # Finalize structure
    net.sortModules()

    # Restore parameters
    for module, text in modules_pending_params:
        if text:
            module.params[:] = list(map(float, text.split(",")))

    for connection, text in connections_pending_params:
        if text:
            connection.params[:] = list(map(float, text.split(",")))

    return net

# Execute
save_network_strict_lxml(net, "detailed_network.xml")
print("Network with class types saved successfully.")

loaded_net = load_network_strict_lxml("detailed_network.xml")
print("Loaded.", loaded_net.activate([1, 2]))
