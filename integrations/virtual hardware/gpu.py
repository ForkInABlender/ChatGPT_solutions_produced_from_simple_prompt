# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition )

"""

This is a Virtual GPU rendering a singular triangle, and moves it around the screen.

This can be done with matplotlib if you know what you're doing.


It's practical benefit: "conform all opengl calls to work within its complexity", making it useful for rendering GPU calculations just as efficiently as on the CPU.

This can be done with 0 GLSL involved, compiled, or embedded.

"""


import sys
import numpy as np
import asyncio
from PyQt5.QtWidgets import QApplication, QWidget
from PyQt5.QtGui import QPainter, QImage
from PyQt5.QtCore import Qt, QTimer
from qasync import QEventLoop  # qasync integrates asyncio with PyQt5

class VirtualGPU:
    def __init__(self, width, height):
        # Initialize the frame buffer as an array of zeros (black)
        self.width = width
        self.height = height
        self.frame_buffer = np.zeros((height, width, 3), dtype=np.uint8)

    def clear(self, color=(0, 0, 0)):
        # Efficiently clear the frame buffer using numpy
        self.frame_buffer[:, :] = color

    def draw_triangle(self, vertices, color):
        # Using bounding box for optimization
        min_x = max(0, min(v[0] for v in vertices))
        max_x = min(self.width, max(v[0] for v in vertices))
        min_y = max(0, min(v[1] for v in vertices))
        max_y = min(self.height, max(v[1] for v in vertices))

        # Optimized rasterization loop using numpy vectorization
        X, Y = np.meshgrid(np.arange(min_x, max_x), np.arange(min_y, max_y))
        X_flat, Y_flat = X.ravel(), Y.ravel()

        # Vectorized point-in-triangle calculation
        v0 = np.array(vertices[2]) - np.array(vertices[0])
        v1 = np.array(vertices[1]) - np.array(vertices[0])
        v2 = np.vstack((X_flat, Y_flat)).T - np.array(vertices[0])

        dot00 = np.dot(v0, v0)
        dot01 = np.dot(v0, v1)
        dot02 = np.einsum('ij,j->i', v2, v0)
        dot11 = np.dot(v1, v1)
        dot12 = np.einsum('ij,j->i', v2, v1)

        inv_denom = 1 / (dot00 * dot11 - dot01 * dot01)
        u = (dot11 * dot02 - dot01 * dot12) * inv_denom
        v = (dot00 * dot12 - dot01 * dot02) * inv_denom

        inside_triangle = (u >= 0) & (v >= 0) & (u + v < 1)
        points_inside = np.vstack((X_flat[inside_triangle], Y_flat[inside_triangle])).T

        # Update only those points inside the triangle
        self.frame_buffer[points_inside[:, 1], points_inside[:, 0]] = color

class RenderWidget(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle('Async Optimized Virtual GPU with PyQt5')
        self.width = 640
        self.height = 480
        self.setFixedSize(self.width, self.height)

        # Set Qt optimization for opaque painting to reduce flickering
        self.setAttribute(Qt.WA_OpaquePaintEvent)

        # Create Virtual GPU
        self.gpu = VirtualGPU(self.width, self.height)

        # Triangle initial position
        self.triangle_vertices = [(320, 100), (100, 380), (540, 380)]
        self.dx = 1  # Reduced step size for smoother motion
        self.dy = 1  # Reduced step size for smoother motion

        # Preallocate the QImage that will be reused for rendering
        self.image = QImage(self.gpu.frame_buffer.data, self.width, self.height, QImage.Format_RGB888)

        # Set up QTimer for consistent updates (around 60 FPS)
        self.timer = QTimer()
        self.timer.timeout.connect(self.update_frame)
        self.timer.start(1)

    def update_frame(self):
        # Update the frame buffer
        self.gpu.clear(color=(0, 0, 0))

        # Draw a red triangle
        self.gpu.draw_triangle(self.triangle_vertices, color=(255, 0, 0))

        # Update triangle position for the next frame
        self.triangle_vertices = [(x + self.dx, y + self.dy) for (x, y) in self.triangle_vertices]

        # Bounce back if it hits the window boundary
        for i, (x, y) in enumerate(self.triangle_vertices):
            if x < 0 or x >= self.width:
                self.dx = -self.dx
            if y < 0 or y >= self.height:
                self.dy = -self.dy

        # Trigger a repaint of the widget
        self.repaint()

    def paintEvent(self, event):
        # Directly paint the QImage without creating new instances
        painter = QPainter(self)
        painter.drawImage(0, 0, self.image)  # Draw the preallocated QImage

async def main():
    app = QApplication(sys.argv)
    loop = QEventLoop(app)
    asyncio.set_event_loop(loop)

    window = RenderWidget()
    window.show()

    with loop:
        loop.run_forever()

if __name__ == '__main__':
    asyncio.run(main())
