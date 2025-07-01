import numpy as np

def ref_peak_cc(Image_mc, peak_pos=None):
    """
    Reference Peak method for image combination.

    Reference:
        Hall, Emma L., et al. "Methodology for improved detection of low 
        concentration metabolites in MRS: optimised combination of signals from 
        multi-element coil arrays." Neuroimage 86 (2014): 35-42.

    Parameters:
        Image_mc (ndarray): Multi-channel image data. Shape: [Nx, Ny, Nz, Ncoil, N_freq, Ndyn]
        peak_pos (int, optional): Index of metabolite peak. If None, selects peak automatically.

    Returns:
        Image_cc (ndarray): Coil-combined image.
        Smap (ndarray): Sensitivity map.

    Written by: Parimal Joshi
    """
    Nx, Ny, Nz, Ncoil, N_freq, Ndyn = Image_mc.shape

    # Sum over dynamics dimension (axis=5)
    Image_t = np.sum(Image_mc, axis=5, keepdims=True)  # shape: [Nx, Ny, Nz, Ncoil, N_freq]
    Spectrum_t = np.squeeze(np.sum(Image_t, axis=(0, 1, 2, 3)))  # shape: [N_freq]

    if peak_pos is None:
        peak_pos = np.argmax(np.abs(Spectrum_t))
        print(f'Auto select metabolite # {peak_pos} for RefPeak combination.')
    else:
        print(f'Manual select metabolite # {peak_pos} for RefPeak combination.')

    Image_ref = Image_t[:, :, :, :, peak_pos]

    Image_sos = np.sqrt(np.sum(np.abs(Image_ref)**2, axis=3, keepdims=True))  
    eps = np.finfo(Image_sos.dtype).eps
    Image_sos[Image_sos == 0] = eps
    # print(f'Image_ref shape: {Image_ref.shape}')
    # print(f'Image_sos shape: {Image_sos.shape}')

    Smap = Image_ref / Image_sos 
    # print(f'Smap shape: {Smap.shape}')
    Smap_exp = np.expand_dims(Smap, axis=(5)) 
    Smap_exp = np.tile(Smap_exp, (1, 1, 1, 1, 1, Ndyn))
    
    Image_cc = np.sum(Image_mc * np.conj(Smap_exp), axis=3)

    return Image_cc, Smap
