# Dylan Kenneth Eliot

"""

Usable with pybrain3, rdkit, brian2, biopython for purposes of mimicking tissue function of the brain.


"""

# Re-import necessary libraries and reattempt 3D scatter plot generation
import numpy as np
import struct
import matplotlib.pyplot as plt

# Reload the raw NIfTI data as done previously
file_path = '/mnt/data/spmT_0001.nii'

# Open the NIfTI file and read its content as raw binary data
with open(file_path, 'rb') as f:
    raw_data = f.read()

# Function to extract the NIfTI header and image data from raw binary content
def parse_nifti(raw_data):
    # NIfTI header is 348 bytes
    header = raw_data[:348]
    # Extract the dimensionality and datatype information from the header
    dim = struct.unpack('<8H', header[40:56])
    
    # The datatype starts at byte 70, 2 bytes
    datatype = struct.unpack('<H', header[70:72])[0]

    # Extract the size of the image data (after the header, starting from byte 352 onwards)
    img_data_offset = 352
    # Calculate the number of elements in the image based on the dimensions
    total_elements = 1
    for d in dim[1:]:  # Skip the first dimension as it's the number of dimensions
        total_elements *= d
    
    # The image data starts after the header
    img_data = raw_data[img_data_offset:img_data_offset + total_elements * 4]  # Assuming 32-bit float
    
    # Convert to numpy array based on datatype (assuming float32)
    image = np.frombuffer(img_data, dtype=np.float32)
    
    # Reshape the image data according to the dimensions
    image = image.reshape(dim[1:])
    
    return dim, datatype, image

# Parse the raw data again
dim, datatype, image = parse_nifti(raw_data)


# The code below is a example for displaying such outside of grayscale as a 3d graphical image.
"""
# Reshape the data to remove singleton dimensions
image_3d = image.squeeze()

# Extract non-zero coordinates from the 3D data
x, y, z = np.nonzero(image_3d)  # Get coordinates of non-zero values
values = image_3d[x, y, z]  # Corresponding values at those coordinates

# Create a 3D scatter plot for the full data
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# Scatter plot the non-zero points
scatter = ax.scatter(x, y, z, c=values, cmap='viridis', marker='o', s=0.5)

# Label the axes and show the color bar
ax.set_xlabel('X Dimension')
ax.set_ylabel('Y Dimension')
ax.set_zlabel('Z Dimension')
fig.colorbar(scatter, ax=ax, label="Voxel Intensity")

plt.title("3D Layout of the Full NIfTI Image Data")
plt.show()
"""
# Another way to use it is like this:
"""
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# Slice along the first dimension (dim 0)
axes[0].imshow(image[60, :, :, 0, 0, 0, 0], cmap='gray')
axes[0].set_title("Slice along Dim 1 (60th Slice)")

# Slice along the second dimension (dim 1)
axes[1].imshow(image[:, 60, :, 0, 0, 0, 0], cmap='gray')
axes[1].set_title("Slice along Dim 2 (60th Slice)")

# Slice along the third dimension (dim 2)
axes[2].imshow(image[:, :, 60, 0, 0, 0, 0], cmap='gray')
axes[2].set_title("Slice along Dim 3 (60th Slice)")

plt.tight_layout()
plt.show()
"""
