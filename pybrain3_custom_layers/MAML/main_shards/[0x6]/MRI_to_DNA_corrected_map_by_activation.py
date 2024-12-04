#Dylan Kenneth Eliot

"""

This is configured for mapping by activation.

All 0's are excluded as they're empty, unused, or not active, or outside the brain itself and filler space.


"""

import numpy as np
import struct

# Load the uploaded NIfTI file to inspect its structure and data
file_path = '/mnt/data/20170831/spmT_0001.nii'

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
        if d == 0:  # Handle cases where the dimension is 0
            break
        total_elements *= d
    
    # The image data starts after the header
    img_data = raw_data[img_data_offset:img_data_offset + total_elements * 4]  # Assuming 32-bit float
    
    # Convert to numpy array based on datatype (assuming float32)
    image = np.frombuffer(img_data, dtype=np.float32)
    
    # Reshape the image data according to the dimensions
    valid_dims = [d for d in dim[1:] if d != 0]  # Filter out 0 dimensions
    image = image.reshape(valid_dims)
    
    return dim, datatype, image

# Parse the raw data
dim, datatype, image = parse_nifti(raw_data)

# Focus only on the 3D spatial data (removing singleton dimensions)
spatial_data = image[:, :, :, 0, 0, 0, 0]
subregion_size = (3, 3, 3)

# Function to group spatial data into (3, 3, 3) subregions
def group_into_subregions(spatial_data, subregion_size):
    """
    Groups the entire 3D spatial data into (3, 3, 3) subregions.
    """
    grouped_data = []

    for x in range(0, spatial_data.shape[0], subregion_size[0]):
        for y in range(0, spatial_data.shape[1], subregion_size[1]):
            for z in range(0, spatial_data.shape[2], subregion_size[2]):
                # Extract (3, 3, 3) subregion from the spatial data
                subregion = spatial_data[x:x + subregion_size[0],
                                         y:y + subregion_size[1],
                                         z:z + subregion_size[2]]
                
                # Add to grouped data with coordinates and subregion values
                grouped_data.append({
                    "coordinates": (x, y, z),
                    "values": subregion
                })

    return grouped_data

# Group the spatial data into subregions
all_subregions = group_into_subregions(spatial_data, subregion_size)

# Function to extract non-zero values and their global coordinates
def extract_non_zero_voxels(all_subregions):
    """
    Extract non-zero voxel values and their global coordinates.
    """
    non_zero_data = []

    for subregion in all_subregions:
        subregion_root_coord = subregion["coordinates"]
        subregion_values = subregion["values"]

        # Find non-zero indices
        non_zero_indices = np.argwhere(subregion_values != -0.0)

        # Map indices to global coordinates in the 3D space
        non_zero_global_coords = [
            (
                subregion_root_coord[0] + idx[0],  # x coordinate
                subregion_root_coord[1] + idx[1],  # y coordinate
                subregion_root_coord[2] + idx[2]   # z coordinate
            )
            for idx in non_zero_indices
        ]

        # Pair global coordinates with their respective values
        non_zero_values_with_coords = [
            {"coordinates": coord, "value": subregion_values[idx[0], idx[1], idx[2]]}
            for coord, idx in zip(non_zero_global_coords, non_zero_indices)
        ]

        non_zero_data.extend(non_zero_values_with_coords)
    
    return non_zero_data

# Extract non-zero voxel data
non_zero_voxel_data = extract_non_zero_voxels(all_subregions)

dna_bases = ['A', 'T', 'C', 'G']

# Normalize voxel values and map to DNA
def voxel_to_dna(voxel_value):
    # Normalize value to range 0-3 (assuming voxel values are between 0 and 1)
    normalized_value = int((voxel_value * 4) % 4)
    return dna_bases[normalized_value]

# Convert all voxel values to DNA sequences
voxel_dna_translation = [
    {
        "coordinates": voxel['coordinates'],
        "dna_sequence": voxel_to_dna(voxel['value'])
    }
    for voxel in non_zero_voxel_data
]

# Display results
for voxel in voxel_dna_translation:
    print(f"Coordinates: {voxel['coordinates']}, DNA Sequence: {voxel['dna_sequence']}")
