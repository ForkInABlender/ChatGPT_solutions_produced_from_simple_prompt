# Dylan Kenneth Eliot

"""

This is the physics for a dynamic shift of spatial displacement of mass and height given quantum tunnelling given a controlled environmental temperature 
 & quantum state of flux. 

When tuned correctly, psig will give the initial point alpha for pressure given medium saturation given ambient temperature,
 with height and mass fluctuation in reference to shift in vacuum pressure or temperature shift even at a quantile level given a finite initial size.
When the point space for alpha is 0.0 and mass is 0.0, you've got a singularity or black-hole, where it hits fission, quantumly shifts given the limit to 
 things moving at the speed of light, or blips from the universe given vacuum pressure or acrossed it.

https://www.igasusa.com/files/R22-PT-Chart.pdf


O5 & O7 Counsel approved for public use and general purpose use.

All rights reserved. 
"""



import sys
import numpy as np
from PyQt5.QtWidgets import QApplication, QMainWindow, QGraphicsView, QGraphicsScene, QWidget, QVBoxLayout, QGraphicsEllipseItem
from PyQt5.QtCore import QTimer, Qt
from PyQt5.QtGui import QColor, QBrush, QPen
from scipy.spatial.transform import Rotation as R

# Constants
hbar = 1.0545718e-34  # Reduced Planck's constant (J·s)
c = 3.0e8  # Speed of light (m/s)
G = 6.67430e-11  # Gravitational constant (m^3/kg/s^2)
v_higgs = 246e9  # Higgs VEV in eV
alpha_initial = 1.2 # Starting suppression factor
beta = 5e6  # Increased suppression strength
mass_initial = 70  # Example mass in kg (human scale)
r_initial = 1.7  # Initial height (m)

# Visual Parameters
BASE_DIAMETER = 100  # Diameter in pixels corresponding to mass_initial
focal_length = 400  # Perspective projection focal length for 2D view
center = np.array([400, 300])  # 2D center of the display area

# 3D Sphere Data (using spherical coordinates)
num_points = 100  # Number of points to sample on the sphere

# Generate 3D sphere points
theta = np.linspace(0, np.pi, num_points)  # polar angle (0 to pi)
phi = np.linspace(0, 2 * np.pi, num_points)  # azimuthal angle (0 to 2pi)

# Generate the 3D coordinates of the sphere
def generate_sphere_points(radius):
    x = radius * np.outer(np.sin(theta), np.cos(phi))
    y = radius * np.outer(np.sin(theta), np.sin(phi))
    z = radius * np.outer(np.cos(theta), np.ones_like(phi))
    return np.vstack((x.flatten(), y.flatten(), z.flatten())).T

# Function for rotating the 3D object
def rotate_points(points_3d, angle_x, angle_y, angle_z):
    rotation = R.from_euler('xyz', [angle_x, angle_y, angle_z], degrees=True)
    rotated_points = rotation.apply(points_3d)
    return rotated_points

# Perspective projection from 3D to 2D
def project_to_2d(points_3d, focal_length=400):
    z = points_3d[:, 2]
    scale = focal_length / (focal_length + z)  # Perspective formula
    x_2d = points_3d[:, 0] * scale
    y_2d = points_3d[:, 1] * scale
    return np.vstack((x_2d, y_2d)).T + center  # Translate to center

# Mass Suppression & Vacuum Control Functions
def mass_suppression(alpha, m):
    return alpha * m

def casimir_effect(alpha):
    return (hbar * c / (240 * (alpha + 1e-20))) ** 0.25

def higgs_suppression(m):
    return 1 / (1 + np.log(1 + m / (v_higgs * 1e-9)))

def spacetime_stability(alpha, m, r):
    epsilon = 1e-20  # To prevent division by zero
    return alpha < (r * c**2) / (2 * G * (m + epsilon))

def dynamic_vacuum_control(iteration, base=1.0):
    # Oscillatory vacuum control factor (can be tuned as needed)
    return base + 0.1 * np.sin(iteration / 10.0)

# Main Window with QGraphicsView
class SimulationWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Dynamic 3D Mass Suppression Simulation")
        self.resize(800, 600)

        # Create central widget with vertical layout
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        layout = QVBoxLayout(central_widget)
        
        # Create a scene and view
        self.scene = QGraphicsScene(0, 0, 800, 600)
        self.view = QGraphicsView(self.scene)
        layout.addWidget(self.view)

        # Simulation state variables
        self.alpha = alpha_initial
        self.mass_current = mass_initial
        self.iteration = 0
        self.radius = r_initial  # Starting size of the object
        
        # Generate initial sphere points
        self.sphere_points = generate_sphere_points(self.radius)
        self.sphere_graphics = []
        
        # Create the visual object (points on sphere surface)
        color = QColor("blue")
        for _ in range(len(self.sphere_points)):
            ellipse = self.scene.addEllipse(0, 0, 5, 5, QPen(Qt.black), QBrush(color))
            self.sphere_graphics.append(ellipse)
        
        # Initial rotation angles
        self.angle_x = 0
        self.angle_y = 0
        self.angle_z = 0

        # Timer for dynamic simulation (updates every 100ms)
        self.timer = QTimer()
        self.timer.timeout.connect(self.update_simulation)
        self.timer.start(100)

    def update_simulation(self):
        # Compute base Casimir effect and modulate with dynamic vacuum control.
        base_a = casimir_effect(self.alpha)
        vac_control_factor = dynamic_vacuum_control(self.iteration)
        a = base_a * vac_control_factor

        # Update alpha using exponential decay and Higgs suppression
        self.alpha = np.exp(-beta * a) * higgs_suppression(self.mass_current)
        self.mass_current = mass_suppression(self.alpha, mass_initial)

        # Check for spacetime stability (if fails, stop simulation)
        if not spacetime_stability(self.alpha, self.mass_current, r_initial):
            self.timer.stop()
            print(f"WARNING: Spacetime collapse risk at step {self.iteration}!")
            return

        # Update radius dynamically (simulating the shrinking and enlarging)
        self.radius = r_initial * (self.mass_current / mass_initial)
        self.radius *= (1 + 0.01 * np.sin(self.iteration / 10.0))  # Additional dynamic factor

        # Generate new 3D points based on the updated radius
        self.sphere_points = generate_sphere_points(self.radius)
        
        # Rotate the sphere
        rotated_points = rotate_points(self.sphere_points, self.angle_x, self.angle_y, self.angle_z)
        
        # Project 3D points to 2D
        points_2d = project_to_2d(rotated_points)

        # Update the visual representation
        self.update_visual(points_2d)

        # Log the current state to the terminal
        self.log_state()
        
        # Increment rotation angles for next frame
        self.angle_x += 1
        self.angle_y += 1
        self.angle_z += 1

        self.iteration += 1

    def update_visual(self, points_2d):
        # Optionally, change color based on alpha (lower alpha -> redder)
        red_component = int(255 * (1 - self.alpha))
        blue_component = int(255 * self.alpha)
        color = QColor(red_component, 0, blue_component)
        
        # Update the position and color of each point
        for i, point in enumerate(points_2d):
            if i < len(self.sphere_graphics):
                self.sphere_graphics[i].setRect(point[0] - 2, point[1] - 2, 4, 4)
                self.sphere_graphics[i].setBrush(QBrush(color))

    def log_state(self):
        # Print out the dynamic state of the simulation to the terminal
        print(f"Iteration: {self.iteration}")
        print(f"Alpha: {self.alpha}")
        print(f"Mass (kg): {self.mass_current}")
        print(f"Radius (m): {self.radius}")
        print(f"Vacuum Control Factor: {dynamic_vacuum_control(self.iteration)}")
        print("-" * 40)

# Main execution
if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = SimulationWindow()
    window.show()
    sys.exit(app.exec_())
