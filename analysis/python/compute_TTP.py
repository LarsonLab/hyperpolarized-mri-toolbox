import numpy as np

def compute_TTP(S, dT):
    """
    Computes the time-to-peak (TTP) for a signal.
    
    Parameters:
    S : np.ndarray
        Signal dynamics, time must be in the last dimension.
    dT : float
        Spacing of time points.
    
    Returns:
    TTP : np.ndarray
        Time-to-peak computed for all voxels and/or metabolites.
    """
    S = np.asarray(S)
    
    # Find the index of the peak value along the last axis
    Imax = np.argmax(S, axis=-1)
    Nt = S.shape[-1]
    
    # Initialize Ipeak with zeros
    Ipeak = np.zeros_like(Imax, dtype=np.float64)
    
    # Quadratic peak fit
    it = np.nditer(Imax, flags=['multi_index'])
    for idx in it:
        n = it.multi_index
        if 0 < Imax[n] < Nt - 1:
            Ipeak[n] = 0.5 * (S[n + (Imax[n] - 1,)] - S[n + (Imax[n] + 1,)]) / \
                        (S[n + (Imax[n] - 1,)] - 2 * S[n + (Imax[n],)] + S[n + (Imax[n] + 1,)])
        else:
            Ipeak[n] = 0
    
    # Compute final TTP
    TTP = (Imax + Ipeak - 1) * dT
    
    return TTP