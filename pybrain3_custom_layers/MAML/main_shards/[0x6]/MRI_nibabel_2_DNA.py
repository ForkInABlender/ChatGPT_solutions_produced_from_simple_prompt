# Dylan Kenneth Eliot

"""
from [x:x', y:y', z:z'], look at the active and inactive regions, and note the genetic markers by location given state of recorded activation captured as MRI NiFTi ".nii" numpy based datasets.



"""

import nibabel as nib
import numpy as np
from sklearn.decomposition import PCA

# Load the NIfTI (.nii) file
nii_file = "20170831/spmT_0001_1.nii"
img = nib.load(nii_file)
data = img.get_fdata()

# 
features = data[76:80, 76:80, 76:80]

pca = PCA(n_components=1)

print("Extracted Features Shape: ", list(map(lambda a: a.tolist()[-1], pca.fit_transform(features.reshape(-1, 1)))))
