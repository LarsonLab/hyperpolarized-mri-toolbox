import numpy as np

def compute_AUCratio(S):
    """
    Computes the Area-under-curve (AUC) ratio for metabolites.
    
    Parameters:
    S : numpy array
        Signal dynamics, must have substrate then product in the second-to-last dimension,
        and time in the last dimension [voxels(optional), metabolites, time points].
    
    Returns:
    numpy array
        AUC ratio calculated as the sum of the product signal divided by the sum of the substrate signal.
    """
    S_shape = S.shape
    
    if len(S_shape) > 2:
        Stemp = S.reshape((-1, 2, S_shape[-1]))
        AUCratio = np.sum(Stemp[:, 1, :], axis=1) / np.sum(Stemp[:, 0, :], axis=1)
        if len(S_shape) > 3:
            AUCratio = AUCratio.reshape(S_shape[:-2])
    else:
        AUCratio = np.sum(S[1, :]) / np.sum(S[0, :])
    
    return AUCratio
