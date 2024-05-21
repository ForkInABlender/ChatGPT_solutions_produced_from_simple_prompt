# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition )

"""
It can handle up to 3d calculations. 

4D calculations may require a different approach against 3d tolopogical order of information. You may consider reshape of data to account for such.


"""

from javascript import require
import numpy as np

# Load the JavaScript file
gpu_browser = require('./js/gpu.js/gpu-browser.js')

# Create a GPU instance
gpu_instance = gpu_browser.GPU()

# Define the kernel function source as a single-line string
kernel_function_source = "function(a, b) { return a[this.thread.z][this.thread.y][this.thread.x] + b[this.thread.z][this.thread.y][this.thread.x]; }"

# Create the kernel function
# Adjust the setOutput dimensions to match the shape of your 3D array
output_shape = [3, 3, 3]  # Example shape, adjust as needed
kernel_function = gpu_instance.createKernel(kernel_function_source).setOutput(output_shape)

# Define the input 3D arrays
array_a = np.random.rand(3, 3, 3)
array_b = np.random.rand(3, 3, 3)

array_a_list = array_a.tolist()
array_b_list = array_b.tolist()

# Call the kernel function with the 3D arrays
result = kernel_function(array_a_list, array_b_list)

# Convert the result back to a NumPy array
#result_array = np.array(result)

print("Resulting 3D array:", result.valueOf())
