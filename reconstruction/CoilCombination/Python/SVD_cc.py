import numpy as np
from scipy.linalg import svd
from skimage.restoration import denoise_wavelet

def SVD_cc(image_mc):
    
    #  SVD based Coil Combination for Hyperpolarized MRI data
    #  Reference:
    #    Hall, Emma L., et al. "Methodology for improved detection of low 
    #    concentration metabolites in MRS: optimised combination of signals from 
    #    multi-element coil arrays." Neuroimage 86 (2014): 35-42.

    #  INPUTS:
    #    Image_mc - multi-channel image data. dim:[Nx, Ny, Nz, Ncoil, N_freq, Ndyn]
        
    #  OUTPUTS:
    #    Image_cc - coil combination Image.

    # (c) 2017-2018 The Regents of the University of California
    # All Rights Reserved.

    # Written by :Parimal Joshi 

    Nx, Ny, Nz, Ncoil, N_freq, Ndyn = image_mc.shape
    image_cc = np.zeros((Nx, Ny, Nz, N_freq, Ndyn), dtype=image_mc.dtype)

    for i in range(Nx):
        for j in range(Ny):
            for k in range(Nz):
                for t in range(Ndyn):
                    tmp_spectrum = image_mc[i, j, k, :, :, t]  # [Ncoil, N_freq]
                    U, s, Vh = svd(tmp_spectrum, full_matrices=False)
                    thresh = s * (np.abs(s) > 1)
                    S = np.diag(thresh)
                    combined = np.sum(U @ S @ Vh, axis=0)  # [N_freq]
                    image_cc[i, j, k, :, t] = combined

    return image_cc
