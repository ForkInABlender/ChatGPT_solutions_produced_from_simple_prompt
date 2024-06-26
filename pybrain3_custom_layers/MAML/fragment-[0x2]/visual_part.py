# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition )

"""
This is the part for visualization. This now needs to be trained on image data.

This works well with MESA & the following when using this script:


LIBGL_ALWAYS_SOFTWARE=1 python3.[6, 10] visual_part.py

The command itself will allow for using mesa instead of your GPU hardware. Now, using such will make your code CPU focused, but it will render.

Currently, this script is setup so a display is not dependent on how it functions.

This is best used with the docker hub container (https://hub.docker.com/r/de3343/ai_mods_py3.10).

For pybrain3 documentation, use (https://web.archive.org/web/20230115194351/http://pybrain.org/pages/home)

For saving and reloading, use the python library "dill", use (https://pypi.org/project/dill/)

"""


import moderngl
import numpy as np
from pybrain3.datasets import SupervisedDataSet
from pybrain3.tools.shortcuts import buildNetwork
from pybrain3.supervised.trainers import BackpropTrainer
from PIL import Image

# Create a ModernGL context in a headless mode
ctx = moderngl.create_standalone_context()

# Define the size of the framebuffer
width, height = 401, 321

# Create a numpy array to hold pixel data
pixels = np.zeros((height, width, 4), dtype=np.uint8)

# Use a for-loop to set the color of each pixel (e.g., a gradient)
for y in range(height):
    for x in range(width):
        pixels[y, x] = [0, 0, 0, 255]  # RGBA format

# Create a texture from the pixel data
texture = ctx.texture((width, height), 4, data=pixels.tobytes())
texture.use()

# Define vertex data for a fullscreen quad
vertices = np.array([
    -1.0, -1.0,
     1.0, -1.0,
    -1.0,  1.0,
     1.0,  1.0,
], dtype='f4')

# Create a vertex buffer object (VBO)
vbo = ctx.buffer(vertices.tobytes())

# Create a vertex array object (VAO)
program = ctx.program(
    vertex_shader='''
    #version 330 core
    in vec2 in_vert;
    out vec2 v_vert;
    void main() {
        v_vert = in_vert;
        gl_Position = vec4(in_vert, 0.0, 1.0);
    }
    ''',
    fragment_shader='''
    #version 330 core
    in vec2 v_vert;
    out vec4 f_color;
    uniform sampler2D Texture;
    void main() {
        f_color = texture(Texture, v_vert);
    }
    '''
)

vao = ctx.simple_vertex_array(program, vbo, 'in_vert')

# Create a framebuffer object (FBO) and attach a texture
fbo = ctx.framebuffer(
    color_attachments=[ctx.texture((width, height), 4)]
)

# Bind the framebuffer
fbo.use()

# Clear the framebuffer with white color
ctx.clear(1.0, 1.0, 1.0, 1.0)

# Bind the texture to the fragment shader
program['Texture'].value = 0

# Render the fullscreen quad
vao.render(moderngl.TRIANGLE_STRIP)

# Read the pixels from the framebuffer
data = fbo.read(components=4)
image = np.frombuffer(data, dtype=np.uint8).reshape((height, width, 4))

# Convert the image array to a writable array
image = image.copy()

# Print the output (for verification)
print("Pixel data at (400, 300):", image[300, 400])  # Sample pixel

# Optionally save the image using an external library (e.g., Pillow)
Image.fromarray(image).save('output.png')

# Prepare the entire image data for pybrain
flattened_image = image[:, :, :3].flatten() / 255.0  # Normalize the RGB values

# Create a dataset
ds = SupervisedDataSet(width * height * 3, width * height * 3)
ds.addSample(flattened_image, flattened_image)  # Use the entire image as input and target

net = buildNetwork(width * height * 3, 10, width * height * 3)  # Input, hidden, output

trainer = BackpropTrainer(net, ds)

for a in range(30):
    print(f"{a+1}: {trainer.train()}")


updated_flattened_image = net.activate(flattened_image) * 255.0
updated_image = updated_flattened_image.reshape((height, width, 3)).astype(np.uint8)

# Replace the RGB channels in the original image with the updated image
image[:, :, :3] = updated_image

# Print the updated output (for verification)
print("Updated pixel data at (400, 300):", image[300, 400])  # Sample pixel

# Save the updated image using Pillow
Image.fromarray(image).save('updated_output.png')

# Release resources
vbo.release()
vao.release()
fbo.release()
ctx.release()
