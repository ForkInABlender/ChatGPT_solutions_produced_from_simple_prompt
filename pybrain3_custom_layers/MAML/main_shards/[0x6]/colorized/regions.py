# Dylan Kenneth Eliot & Julius.ai

"""
This was tested with julius. It aided in shaping how to see it in full color per region.
Each region in the CSV file is mapped to a highlighted region. This lets one better isolate for the subnetwork regions that need to trained on
 via pybrain3 inside of a docker image.


"""

import nibabel as nib
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.colors import ListedColormap
import numpy as np
from matplotlib.colors import rgb2hex


# Load the NIfTI file
img = nib.load('temp.nii')

# Fetch the MNI152 template with brain regions labeled
atlas = datasets.fetch_atlas_harvard_oxford('cort-maxprob-thr25-1mm')
atlas_filename = atlas.maps

# Resample the atlas to match the MRI image dimensions
resampled_atlas = image.resample_to_img(atlas_filename, img, interpolation='nearest')
display = plotting.plot_roi(resampled_atlas,bg_img=img,display_mode='ortho',title='Brain Regions -- gray 2 full color',colorbar=True,cmap='Spectral')


# Extract the colormap used in the visualization
cmap = plt.cm.Spectral

# Get the labels from the atlas
region_labels = atlas.labels

# Generate a color for each region (excluding background)
num_regions = len(region_labels) - 1  # Exclude the background
colors = cmap(np.linspace(0, 1, num_regions))

# Create a DataFrame mapping regions to colors
region_color_mapping = pd.DataFrame({
    'Region': region_labels[1:],  # Exclude the background label
    'Color': [rgb2hex(color[:3]) for color in colors]  # Convert RGB to HEX
})

# Create a visual representation of the colors
fig, ax = plt.subplots(figsize=(12, 8))
ax.axis('tight')
ax.axis('off')

# Create table
table = ax.table(cellText=region_color_mapping.values,
                colLabels=region_color_mapping.columns,
                cellLoc='left',
                loc='center',
                cellColours=[[color, color] for color in colors])

# Adjust table properties
table.auto_set_font_size(False)
table.set_fontsize(9)
table.scale(1.2, 1.5)

region_color_mapping.to_csv('brain_region_colors.csv', index=False)

plt.show()
