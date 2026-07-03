import numpy as np

def compute_mean_time(S, dT):
    """
    Computes the mean time for a signal, defined as the center of mass of the time curve.
    
    Parameters:
    S : numpy array
        Signal dynamics, with time in the last dimension [voxels and/or # of metabolites, # of time points].
    dT : float
        Spacing of time points.
    
    Returns:
    numpy array
        Mean time computed for all voxels and/or metabolites.
    """
    S_shape = S.shape
    Nt = S_shape[-1]  # Number of time points
    
    t_values = np.arange(Nt) * dT
    
    if len(S_shape) > 1:
        tmat = np.reshape(t_values, [*(1,) * (len(S_shape) - 1), Nt])
        tmat = np.broadcast_to(tmat, S_shape)  # Match dimensions of S
    else:
        tmat = t_values
    
    mean_time = np.sum(S * tmat, axis=-1) / np.sum(S, axis=-1)
    
    return mean_time
