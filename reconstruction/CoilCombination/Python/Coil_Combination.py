import numpy as np
import SVD_cc
import prewhitening_cc
import Sum_of_Square_cc
import RefPeak_cc

def coil_combination(Image_mc, cc_method_flag=0, prewhitening_flag=0):
# Coil Combination for Hyperpolarized MRI data
# Reference:
#    Hall, Emma L., et al. "Methodology for improved detection of low 
#    concentration metabolites in MRS: optimised combination of signals from 
#    multi-element coil arrays." Neuroimage 86 (2014): 35-42.

# INPUTS:
#    Image_mc - multi-channel image data. dim:[Nx, Ny, Nz, Ncoil, N_freq, Ndyn]
#    cc_method_flag -
#        0 - sum of square (SOS).
#        1 - Singular value decomposition(SVD). Only works for spetroscopy
#        data
#        2 - Reference Peak method.
#    prewhitening_flag - 
#        0 - No prewhitening
#        1 - Prewhitening based on edge of image       
#  OUTPUTS:
#    Image_cc - coil combination Image.
#    Smap - sensitivity maps, empty for 0,1,2 results.

#  (c) 2017-2018 The Regents of the University of California
#  All Rights Reserved.

#  Written by :Parimal Joshi 

    if prewhitening_flag:
        # Assuming the last time point is noise
        Image_mc, correlation_matrix = prewhitening_cc(Image_mc)

    Nx, Ny, Nz, Ncoil, N_freq, Ndyn = Image_mc.shape

    # Coil combination method
    if cc_method_flag == 0:
        Image_cc = Sum_of_Square_cc.sum_of_square_cc(Image_mc)
        Smap = None
    elif cc_method_flag == 1:
        Image_cc = SVD_cc.svd_cc(Image_mc)
        Smap = None
    elif cc_method_flag == 2:
        Image_cc, Smap = RefPeak_cc.ref_peak_cc(Image_mc)
    else:
        raise ValueError("Invalid cc_method_flag. Must be 0, 1, or 2.")

    return Image_cc, Smap