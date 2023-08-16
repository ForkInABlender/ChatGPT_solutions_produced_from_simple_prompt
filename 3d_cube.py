# Dylan Kenneth Eliot & gpt-4-codeinterpreter

"""
Creates a basic 3d_cube and adds a texture image, then rotates the cube.

This a okay starting point for those looking to learn OpenGL programming.
"""

import pygame
from pygame.locals import *
from OpenGL.GL import *
from OpenGL.GLU import *
from PIL import ImageSequence, Image
import numpy as np

# Load the gif image
image_path = 'image_texture.png'
image = Image.open(image_path)

# Get all the frames
frames = [frame.copy() for frame in ImageSequence.Iterator(image)]

# Load one of the frames as the texture
texture_image = frames[0].convert("RGB")
texture_data = np.array(texture_image)

# Initialize Pygame
pygame.init()
display = (800, 600)
pygame.display.set_mode(display, DOUBLEBUF | OPENGL)

# Perspective projection

def gluPerspective(fov_y, aspect_ratio, near, far):
	height = near * np.tan(np.radians(fov_y) / 2)
	width = height * aspect_ratio
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
	glMatrixMode(GL_PROJECTION)
	glLoadIdentity()
	glFrustum(-width, width, -height, height, near, far)
	glMatrixMode(GL_MODELVIEW)

gluPerspective(45, float(display[0]/display[1]), 0.1, 50.0)
glTranslatef(0.0, 0.0, -5)

# Create texture from image
texture = glGenTextures(1)
glBindTexture(GL_TEXTURE_2D, texture)
glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, texture_image.width, texture_image.height, 0, GL_RGB, GL_UNSIGNED_BYTE, texture_data)
glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR)
glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR)
# Enable texture and depth test
glEnable(GL_TEXTURE_2D)
glEnable(GL_DEPTH_TEST)
# Main loop
while True:
		for event in pygame.event.get():
			if event.type == pygame.QUIT:
				pygame.quit()
				quit()
		# Clear buffers
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
		# Define the cube vertices and texture coordinates
		vertices = [(-1, -1,  1), ( 1, -1,  1), ( 1,  1,  1), (-1,  1,  1),(-1, -1, -1), ( 1, -1, -1), ( 1,  1, -1), (-1,  1, -1)]
		tex_coords = [(0, 0), (1, 0), (1, 1), (0, 1)]
		faces = [(0, 1, 2, 3), (3, 2, 6, 7), (7, 6, 5, 4), (4, 5, 1, 0), (0, 3, 7, 4), (1, 5, 6, 2)]
		# Draw the textured cube
		glBindTexture(GL_TEXTURE_2D, texture)
		glBegin(GL_QUADS)
		for face in faces:
			for i, vertex in enumerate(face):
				glTexCoord2fv(tex_coords[i])
				glVertex3fv(vertices[vertex])
		glEnd()
		# Rotate the cube
		glRotatef(1, 3, 1, 1)
		pygame.display.flip()
		pygame.time.wait(10)
