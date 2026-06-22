import numpy as np

def prewhitening_cc(image_mc, noise_data=None):
    """
    Prewhitening step for multi-channel image data.
    
    Reference:
    Hall, Emma L., et al. "Methodology for improved detection of low 
    concentration metabolites in MRS: optimised combination of signals from 
    multi-element coil arrays." Neuroimage 86 (2014): 35-42.
    
    Parameters:
    image_mc : numpy.ndarray
        Multi-channel image data with dimensions [Nx, Ny, Nz, Ncoil, N_freq, Ndyn].
    noise_data : numpy.ndarray, optional
        Noise-only data. If not supplied, defaults to using the final timepoint from all
        frequencies and all slices from image_mc.
    
    Returns:
    tuple:
        - image_w : numpy.ndarray
            Prewhitened image data.
        - correlation_matrix : numpy.ndarray
            Correlation matrix.

    Written by: Parimal Joshi
    """
    Nx, Ny, Nz, Ncoil, N_freq, Ndyn = image_mc.shape
    
    if noise_data is None:
        noise = image_mc[0:1, :, :, :, :, -1:] 
    else:
        noise = noise_data
    
    noise = np.transpose(noise, (3, 0, 1, 2, 4, 5))  # [Ncoil, Nx, Ny, Nz, N_freq, Ndyn]
    nsample = noise.size // Ncoil
    noise = noise.reshape(Ncoil, nsample)
    
    eps = np.finfo(noise.dtype).eps
    correlation_matrix = np.abs(noise @ noise.T.conj()) / (nsample + eps)
    
    # Eigen decomposition
    eigvals, eigvecs = np.linalg.eigh(correlation_matrix)
    # Inverse square root of the covariance matrix
    inv_sqrt_eigvals = np.diag(1.0 / np.sqrt(np.maximum(eigvals, 1e-12)))
    Whitening_Matrix = eigvecs @ inv_sqrt_eigvals @ eigvecs.T.conj()

    sigma = np.sqrt(np.diag(correlation_matrix)).reshape(1, 1, 1, Ncoil, 1, 1)
    sigma = np.tile(sigma, (Nx, Ny, Nz, 1, N_freq, Ndyn))
    sigma = np.clip(sigma, 1e-6, None)
    
    image_w = image_mc / (sigma + eps)
    image_whitened = np.tensordot(image_w, Whitening_Matrix, axes=([3], [1]))
    image_whitened = np.moveaxis(image_whitened, -1, 3)  # Ensure [Nx, Ny, Nz, Ncoil, N_freq, Ndyn]


    return image_whitened, correlation_matrix