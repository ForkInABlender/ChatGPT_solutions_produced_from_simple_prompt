# Dylan Kenneth Eliot

"""
compute all 3d chunks and possibilities against a filter, aka "kernel set for 'embedding'..."....


This will later be refactored for pybrain3 network usage. For now, though, as tested code, it exists here as a testiment to what is needed part-wise for
 building someone an AI of their own mind. 

"""

import numpy as np


def f(b, d_start, d_end, h_start, h_end, w_start, w_end, input_tensor):
 return input_tensor[b, d_start:d_end, h_start:h_end, w_start:w_end, :]

#g-gemini 2024/12/09
def conv3d_with_list_map(input_tensor, kernel, stride=1, padding=0):
    """
    Performs 3D convolution on the input tensor using list, map, and lambda.

    Args:
        input_tensor: Input tensor with shape (batch_size, depth, height, width, channels).
        kernel: Convolution kernel with shape (kernel_depth, kernel_height, kernel_width, in_channels, out_channels).
        stride: Stride for convolution.
        padding: Amount of zero padding to apply to the input tensor.

    Returns:
        Output tensor after convolution.
    """

    batch_size, depth, height, width, channels = input_tensor.shape
    kernel_depth, kernel_height, kernel_width, in_channels, out_channels = kernel.shape

    # Apply padding
    if padding > 0:
        input_tensor = np.pad(input_tensor, ((0,), (padding,), (padding,), (padding,), (0,)))

    # Calculate output dimensions
    output_depth = (depth - kernel_depth) // stride + 1
    output_height = (height - kernel_height) // stride + 1
    output_width = (width - kernel_width) // stride + 1

    # Initialize the output tensor
    output = np.zeros((batch_size, output_depth, output_height, output_width, out_channels))

    # Perform convolution
    for b in range(batch_size):
        for d in range(output_depth):
            for h in range(output_height):
                for w in range(output_width):
                    output[b, d, h, w, :] = list(map(
                        lambda c: np.sum(f(b, d*stride, d*stride+kernel_depth, h*stride, h*stride+kernel_height, w*stride, w*stride+kernel_width, input_tensor) * kernel[:, :, :, :, c]),
                        range(out_channels)
                    ))

    return output

# Example usage
input_tensor = np.random.rand(1, 121, 145, 121, 1)  # 1 sample, 121x145x121 volume, 1 channel(/s)
kernel = np.random.rand(3, 3, 3, 1, 1)  # 3x3x3 kernel, 1 input channel(/s), 1 output channel(/s)
output_tensor = conv3d_with_list_map(input_tensor, kernel)

print(kernel.shape)
print(input_tensor.shape)
print(output_tensor.shape)
Ax=conv3d_with_list_map(output_tensor, kernel)
print(Ax.shape)
print("testing:")
for AA in range(10):
    Ax=conv3d_with_list_map(Ax, kernel)
    print(Ax.shape)
