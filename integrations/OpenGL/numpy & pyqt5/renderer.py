# Dylan Kenneth Eliot


import sys
import numpy as np
from PyQt5.QtWidgets import QApplication, QMainWindow, QFrame
from PyQt5.QtGui import QPainter, QColor
from PyQt5.QtCore import Qt, QTimer

"""

"""


class PyOpenGLEmulator:
    def __init__(self):
        self.model_matrix = np.identity(4, dtype=np.float32)
        self.view_matrix = np.identity(4, dtype=np.float32)
        self.projection_matrix = np.identity(4, dtype=np.float32)
        self.vertices = []
        self.transformed_vertices = []
        self.matrix_stack = []
        self.current_matrix_mode = 'MODEL'

    def glLoadIdentity(self):
        """Resets the current matrix to the identity matrix."""
        if self.current_matrix_mode == 'MODEL':
            self.model_matrix = np.identity(4, dtype=np.float32)
        elif self.current_matrix_mode == 'VIEW':
            self.view_matrix = np.identity(4, dtype=np.float32)
        elif self.current_matrix_mode == 'PROJECTION':
            self.projection_matrix = np.identity(4, dtype=np.float32)

    def glMatrixMode(self, mode):
        """Sets the current matrix mode (GL_MODELVIEW, GL_VIEW, or GL_PROJECTION)."""
        if mode in ['GL_MODELVIEW', 'GL_VIEW', 'GL_PROJECTION']:
            if mode == 'GL_MODELVIEW':
                self.current_matrix_mode = 'MODEL'
            elif mode == 'GL_VIEW':
                self.current_matrix_mode = 'VIEW'
            elif mode == 'GL_PROJECTION':
                self.current_matrix_mode = 'PROJECTION'
        else:
            raise ValueError("Invalid matrix mode. Use 'GL_MODELVIEW', 'GL_VIEW', or 'GL_PROJECTION'.")

    def glPushMatrix(self):
        """Pushes the current matrix onto the stack."""
        if self.current_matrix_mode == 'MODEL':
            self.matrix_stack.append(np.copy(self.model_matrix))
        elif self.current_matrix_mode == 'VIEW':
            self.matrix_stack.append(np.copy(self.view_matrix))
        elif self.current_matrix_mode == 'PROJECTION':
            self.matrix_stack.append(np.copy(self.projection_matrix))

    def glPopMatrix(self):
        """Pops the top matrix from the stack and sets it as the current matrix."""
        if not self.matrix_stack:
            raise ValueError("Matrix stack is empty.")
        if self.current_matrix_mode == 'MODEL':
            self.model_matrix = self.matrix_stack.pop()
        elif self.current_matrix_mode == 'VIEW':
            self.view_matrix = self.matrix_stack.pop()
        elif self.current_matrix_mode == 'PROJECTION':
            self.projection_matrix = self.matrix_stack.pop()

    def glTranslatef(self, x, y, z):
        """Translates the current matrix by x, y, z."""
        translation_matrix = np.array([
            [1, 0, 0, x],
            [0, 1, 0, y],
            [0, 0, 1, z],
            [0, 0, 0, 1]
        ], dtype=np.float32)
        if self.current_matrix_mode == 'MODEL':
            self.model_matrix = np.dot(self.model_matrix, translation_matrix)
        elif self.current_matrix_mode == 'VIEW':
            self.view_matrix = np.dot(self.view_matrix, translation_matrix)
        elif self.current_matrix_mode == 'PROJECTION':
            self.projection_matrix = np.dot(self.projection_matrix, translation_matrix)

    def glRotatef(self, angle, x, y, z):
        """Rotates the current matrix around the given axis by an angle (in degrees)."""
        rad = np.deg2rad(angle)
        cos_theta, sin_theta = np.cos(rad), np.sin(rad)
        if x == 0 and y == 0 and z == 1:  # Z-axis rotation
            rotation_matrix = np.array([
                [cos_theta, -sin_theta, 0, 0],
                [sin_theta, cos_theta,  0, 0],
                [0,         0,          1, 0],
                [0,         0,          0, 1]
            ], dtype=np.float32)
        elif x == 1 and y == 0 and z == 0:  # X-axis rotation
            rotation_matrix = np.array([
                [1, 0,         0,          0],
                [0, cos_theta, -sin_theta, 0],
                [0, sin_theta, cos_theta,  0],
                [0, 0,         0,          1]
            ], dtype=np.float32)
        elif x == 0 and y == 1 and z == 0:  # Y-axis rotation
            rotation_matrix = np.array([
                [cos_theta,  0, sin_theta, 0],
                [0,          1, 0,         0],
                [-sin_theta, 0, cos_theta, 0],
                [0,          0, 0,         1]
            ], dtype=np.float32)
        else:
            raise ValueError("Rotation axis must be one of the unit axes (X, Y, Z).")
        
        if self.current_matrix_mode == 'MODEL':
            self.model_matrix = np.dot(self.model_matrix, rotation_matrix)
        elif self.current_matrix_mode == 'VIEW':
            self.view_matrix = np.dot(self.view_matrix, rotation_matrix)
        elif self.current_matrix_mode == 'PROJECTION':
            self.projection_matrix = np.dot(self.projection_matrix, rotation_matrix)

    def glScalef(self, x, y, z):
        """Scales the current matrix by x, y, z."""
        scale_matrix = np.array([
            [x, 0, 0, 0],
            [0, y, 0, 0],
            [0, 0, z, 0],
            [0, 0, 0, 1]
        ], dtype=np.float32)
        if self.current_matrix_mode == 'MODEL':
            self.model_matrix = np.dot(self.model_matrix, scale_matrix)
        elif self.current_matrix_mode == 'VIEW':
            self.view_matrix = np.dot(self.view_matrix, scale_matrix)
        elif self.current_matrix_mode == 'PROJECTION':
            self.projection_matrix = np.dot(self.projection_matrix, scale_matrix)

    def glVertex3f(self, x, y, z):
        """Adds a vertex to the list of vertices."""
        vertex = np.array([x, y, z, 1.0], dtype=np.float32)
        self.vertices.append(vertex)

    def glBegin(self, mode):
        """Begins a new shape (not implemented in detail, only placeholder)."""
        self.vertices = []

    def glEnd(self):
        """Ends the current shape (not implemented in detail, only placeholder)."""
        self.transform_vertices()

    def transform_vertices(self):
        """Applies the model, view, and projection matrices to all vertices."""
        self.transformed_vertices = []
        for vertex in self.vertices:
            transformed_vertex = np.dot(self.model_matrix, vertex)
            transformed_vertex = np.dot(self.view_matrix, transformed_vertex)
            transformed_vertex = np.dot(self.projection_matrix, transformed_vertex)
            self.transformed_vertices.append(transformed_vertex)

    def glDraw(self):
        """Returns the transformed vertices, simulating a draw call."""
        if not self.transformed_vertices:
            self.transform_vertices()
        return self.transformed_vertices

class OpenGLFrame(QFrame):
    def __init__(self, parent=None):
        super(OpenGLFrame, self).__init__(parent)
        self.emulator = PyOpenGLEmulator()
        self.setMinimumSize(400, 400)
        self.translation_x = 0

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.Antialiasing)
        painter.fillRect(event.rect(), Qt.black)
        self.draw_scene(painter)

    def draw_scene(self, painter):
        """Draws the transformed vertices on the frame."""
        self.emulator.glLoadIdentity()
        self.emulator.glTranslatef(self.translation_x, 0, 0)
        self.emulator.glRotatef(45, 0, 0, 1)
        self.emulator.glScalef(100, 100, 1)

        # Adding vertices of a simple triangle
        self.emulator.glBegin('GL_TRIANGLES')
        self.emulator.glVertex3f(0, 1, 0)
        self.emulator.glVertex3f(-1, -1, 0)
        self.emulator.glVertex3f(1, -1, 0)
        self.emulator.glEnd()

        painter.setPen(QColor(0, 255, 0))
        transformed_vertices = self.emulator.glDraw()
        if len(transformed_vertices) == 3:
            v1, v2, v3 = transformed_vertices
            # Transform from model space to screen space (basic approach)
            width, height = self.width(), self.height()
            v1_screen = (int(v1[0] + width / 2), int(height / 2 - v1[1]))
            v2_screen = (int(v2[0] + width / 2), int(height / 2 - v2[1]))
            v3_screen = (int(v3[0] + width / 2), int(height / 2 - v3[1]))
            painter.drawLine(v1_screen[0], v1_screen[1], v2_screen[0], v2_screen[1])
            painter.drawLine(v2_screen[0], v2_screen[1], v3_screen[0], v3_screen[1])
            painter.drawLine(v3_screen[0], v3_screen[1], v1_screen[0], v1_screen[1])

class MainWindow(QMainWindow):
    def __init__(self):
        super(MainWindow, self).__init__()
        self.setWindowTitle("PyQt5 OpenGL Emulator")
        self.setGeometry(100, 100, 600, 600)
        self.opengl_frame = OpenGLFrame(self)
        self.setCentralWidget(self.opengl_frame)

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()

    # Initialize the timer to make the triangle move around the screen and update the position
    def update_position():
        window.opengl_frame.translation_x += 1
        if window.opengl_frame.translation_x > window.opengl_frame.width() / 2:
            window.opengl_frame.translation_x = -window.opengl_frame.width() / 2
        window.opengl_frame.update()

    timer = QTimer()
    timer.timeout.connect(update_position)
    timer.start(16)  # Approximately 60 frames per second

    sys.exit(app.exec_())
