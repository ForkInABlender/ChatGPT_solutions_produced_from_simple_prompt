# Dylan Kenneth Eliot 


from pybrain3.tools.xml import NetworkReader

net = NetworkReader.readFrom("CCO_spike.xml")

print(net.activate([.1, .1]))
