# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition )


"""
This allows for setup code in one place and the running of it in another.



This is approved by O5 & O7 Counsel regulation for SCP founders.

Access granted to all public, private, or D-class personnel. 

Product:Safe
"""

import xmlrpc.client

# Connect to the server
server_url = 'http://localhost:8000'
with xmlrpc.client.ServerProxy(server_url) as proxy:
    # Build the network on the server
    print("Building network:", proxy.buildNetwork(2, 3, 1))
    
    # Create the dataset on the server
    print("Creating dataset:", proxy.createDataset(2, 1))
    
    # Add samples to the dataset
    print("Adding sample (0,0)->0:", proxy.addSample([0, 0], [0]))
    print("Adding sample (0,1)->1:", proxy.addSample([0, 1], [1]))
    print("Adding sample (1,0)->1:", proxy.addSample([1, 0], [1]))
    print("Adding sample (1,1)->0:", proxy.addSample([1, 1], [0]))
    
    # Train the network
    print("Training network:", proxy.trainNetwork(1000))
    
    # Make predictions
    inputs = [(0, 0), (0, 1), (1, 0), (1, 1)]
    for x, y in inputs:
        result = proxy.predict([x, y])
        print(f"Predict({x}, {y}): {result}")
    
    # Call the currentTime method from PyBrainService.currentTime
    print(f"Current Time: {proxy.currentTime.getCurrentTime()}")
