# Dylan Kenneth Eliot

"""
Parts needed for making a brain tape from fMRI/MRI images.

This will be used with pybrain3 for playback after snapshot.

"""

import nibabel as nib
import numpy as np
from scipy.ndimage import zoom

def load_nifti_to_array(filepath, target_shape=(40, 32, 32)):
    # Load NIfTI file and convert to numpy array
    img = nib.load(filepath)
    data = img.get_fdata()
    # Normalize data to range 0-1
    data = (data - np.min(data)) / (np.max(data) - np.min(data))
    # Resize to target shape
    zoom_factors = [t / s for t, s in zip(target_shape, data.shape)]
    resized_data = zoom(data, zoom_factors, order=1)  # Use order=1 for linear interpolation
    return resized_data

# Load and preprocess data for a specific brain region/whole brain then segment to region
# Example: loading the hippocampus region segmented file
hippocampus_data = load_nifti_to_array("path/to/hippocampus.nii.gz")
