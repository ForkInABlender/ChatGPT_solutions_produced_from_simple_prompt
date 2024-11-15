# Dylan Kenneth Eliot

"""

This is how to read the MRI/fMRI data raw in the terminal. All points that are zero tell you "nothing physically present" or "no activation of region/function(/s) if any".

From this point forward, the development path has been accelerated. Meaning we'll know how to apply the values the snapshots via MRI/fMRI data as well as eeg.

The MRI/fMRI data will only increase the precision it gets during training.

This will be used with the brain emulator and neuro receptors & transmitters that are also emulateable.


"""

import nibabel as nib
import numpy as np

def load_and_visualize_nifti(file_path, slice_index=60):
    nii = nib.load(file_path)
    data = nii.get_fdata()
    
    print("Data dimensions:", data.shape)
    print("Data type:", data.dtype)
    
    # Display a middle slice
    print(list(data[: :, slice_index]))

def main():
    path_1 = './20170831/spmT_0001.nii'
    path_2 = './20170831/spmT_0001_1.nii'
    
    load_and_visualize_nifti(path_1)
    load_and_visualize_nifti(path_2)

if __name__ == '__main__':
    main()
