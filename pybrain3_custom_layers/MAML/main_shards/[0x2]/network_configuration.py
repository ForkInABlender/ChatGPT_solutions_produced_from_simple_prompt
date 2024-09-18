# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition )

"""

Basic configuration and setup :: To be revised later via config & training sets.


Product status: Keter ( Blue stage; push for Green before production ready )


"""

from pybrain3.structure import TanhLayer, LSTMLayer, FullConnection
from pybrain3.tools.shortcuts import buildNetwork
from pybrain3.structure.networks import RecurrentNetwork
from configparser import ConfigParser

# Load the configuration file
config = ConfigParser()
config.read('network_config.ini')

# Function to manually define a network
def manually_define_network(section):
    net = RecurrentNetwork()
    # split data by network layer type, layer node count, layer name
    in_layer=tuple(config[section]['input_module'].split(', '))
    hidden_layer=tuple(config[section]['hidden_module0'].split(', '))
    out_layer=tuple(config[section]['output_module'].split(', '))
    # Add layers
    net.addInputModule(eval(in_layer[0])(eval(in_layer[1]), name=in_layer[2]))
    net.addModule(eval(hidden_layer[0])(eval(hidden_layer[1]), name=hidden_layer[2]))
    net.addOutputModule(eval(out_layer[0])(eval(out_layer[1]), name=out_layer[2]))
    #split data by connection type, layer in, layer out
    in2h0=tuple(config[section]['connection_in_to_hidden0'].split(', '))
    h02out=tuple(config[section]['connection_hidden0_to_out'].split(', '))
    recurrent_hidden=tuple(config[section]['recurrent_connection_hidden0'].split(', '))
    # Add connections
    net.addConnection(eval(in2h0[0])(eval('net[\''+in2h0[1]+'\']'), eval('net[\''+in2h0[2]+'\']')))
    net.addConnection(eval(h02out[0])(eval('net[\''+h02out[1]+'\']'), eval('net[\''+h02out[2]+'\']')))
    net.addRecurrentConnection(eval(recurrent_hidden[0])(eval('net[\''+recurrent_hidden[1]+'\']'), eval('net[\''+recurrent_hidden[2]+'\']')))
    
    net.sortModules()
    return net

# Function to define a network using buildNetwork
def build_network(section):
    params = list(map(int, config[section]['structure'].split(', ')))
    hiddenclass = eval(config[section].get('hidden_class', 'TanhLayer'))
    return buildNetwork(*params, hiddenclass=hiddenclass)

# Dictionary to store networks
networks = {}

# Iterate through sections to setup each network
for section in config.sections():
    method = config[section]['defined_method']
    if method == 'manually':
        networks[section] = manually_define_network(section)
    elif method == 'buildNetwork':
        networks[section] = build_network(section)

# Now networks dictionary contains all the neural networks
for name, net in networks.items():
    print(f"Initialized {name} using {config[name]['defined_method']}")
