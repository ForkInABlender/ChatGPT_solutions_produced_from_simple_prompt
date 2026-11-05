# Dylan Kenneth Eliot

"""

Code has been fixed.

make sure you have javascript module installed and are using nodejs 18+ to nodejs 22;

This also works on android OS via python3.10

You will need to make sure libnode-dec is removed/purged correctly vefore proper install takes place. 

https://pypi.org/project/javascript/







"""


from javascript import require

# Load the JavaScript file
gpu_browser = require('./js/gpu.js/gpu-browser.js')

# Create a GPU instance
gpu_instance = gpu_browser.GPU()

# Define the kernel function source as a single-line string
kernel_function_source = "function(a, b) { return a[this.thread.x] + b[this.thread.x]; }"

# Create the kernel function
kernel_function = gpu_instance.createKernel(kernel_function_source).setOutput([10])

# Define the input arrays
array_a = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
array_b = [10, 9, 8, 7, 6, 5, 4, 3, 2, 1]

# Run the kernel function
result = kernel_function(array_a, array_b)

# Print the result
print(result)
