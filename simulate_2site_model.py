import numpy as np
from scipy.linalg import expm
import matplotlib.pyplot as plt

def simulate_2site_model(Tin, R1, k, flips, TR, input_function=None):
    """
    Simulates the magnetization evolution in a two-site exchange model 
    with varying flip angles, following the original 2x2 matrix exponential 
    approach in the MATLAB code.

    Parameters
    ----------
    Tin : float or array_like
        - If a single float, the model starts with magnetization [1, 0]
          (all in site 1), then evolves for 'Tin' seconds before the first flip.
        - If a 2-element vector, it is directly interpreted as [Mz1, Mz2] 
          at the start (before the first flip).
    R1 : float or array_like
        Longitudinal relaxation rates (1/s) for the two sites.
        - If a single float, it will be replicated for both sites.
        - Otherwise, must be length 2, e.g. [R1_site1, R1_site2].
    k : array_like
        Exchange rates, shape = (2,) or (2,). 
        Interpreted as:
            k(0) = forward rate (site1 -> site2)
            k(1) = reverse rate (site2 -> site1)
        Example: k = [0.05, 0.02]
    flips : ndarray
        A 2 x N array (two rows, N columns) specifying the flip angles 
        in radians for each site across N time points:
            flips[0, :] for site1
            flips[1, :] for site2
    TR : float
        Repetition time (seconds) between flips.
    input_function : array_like, optional
        Additional input magnetization to be added to site 1 at each TR. 
        Must have length N if used. If None or all zeros, no extra input.

    Returns
    -------
    Mxy : ndarray
        A (2 x N) array with the transverse magnetization immediately after
        each flip, for each site.
    Mz : ndarray
        A (2 x N) array with the longitudinal magnetization immediately after
        each flip, for each site.

    Notes
    -----
    - This reproduces the original 2-site logic found at the end of the
      MATLAB function simulate_2site_model.m (i.e., the part after the 
      `return` in Peder Larson's code).
    - If you prefer a more general approach (up to 4 sites), see simulate_Nsite_model().
    """
    # Convert everything to NumPy arrays
    flips = np.asarray(flips, dtype=float)
    if flips.shape[0] != 2:
        raise ValueError("For a 2-site model, flips must have 2 rows (one for each site).")

    # Handle R1
    if R1 is None:
        R1 = [0.0, 0.0]
    elif np.isscalar(R1):
        R1 = [R1, R1]
    R1 = np.asarray(R1, dtype=float)
    if R1.size != 2:
        raise ValueError("R1 must be a scalar or length-2 array.")

    # Handle k
    k = np.asarray(k, dtype=float)
    if k.size != 2:
        raise ValueError("k must be length 2 for a two-site model (forward & reverse).")

    # Number of time points (flips)
    N = flips.shape[1]

    # Build system matrix A and exponentiate
    # A = [ -R1(1) - k(1),      +k(2)
    #       +k(1),        -R1(2) - k(2) ]
    A = np.array([
        [-R1[0] - k[0],  +k[1]     ],
        [ +k[0],         -R1[1] - k[1]]
    ], dtype=float)
    Ad = expm(A * TR)  # evolution over TR without extra input

    # Check input_function
    if input_function is None:
        input_function = np.zeros(N, dtype=float)
    else:
        input_function = np.asarray(input_function, dtype=float)
        if input_function.size != N:
            raise ValueError("input_function must have length = N or be None/zeros.")
    use_input_function = not np.allclose(input_function, 0.0)

    # If using input, precompute Bd = A^(-1)*(Ad - I)
    # from the original code: Bd = A\(Ad-eye(2)) 
    # translates to Bd = np.linalg.inv(A) @ (Ad - np.eye(2)).
    if use_input_function:
        Bd = np.linalg.inv(A) @ (Ad - np.eye(2))

    # Determine initial magnetization Mz0
    # If Tin is a scalar => evolve from [1, 0] for 'Tin' seconds
    # If 2-element vector => directly use [Mz1, Mz2]
    if np.isscalar(Tin):
        Tin = float(Tin)
        Mz_init = np.array([1.0, 0.0], dtype=float)
        Mz0 = expm(A * Tin) @ Mz_init
    else:
        Tin = np.asarray(Tin, dtype=float)
        if Tin.size != 2:
            raise ValueError("Tin must be scalar or a 2-element vector.")
        Mz0 = Tin

    # Allocate output arrays
    Mxy = np.zeros((2, N), dtype=float)
    Mz  = np.zeros((2, N), dtype=float)

    # Apply first flip
    # Mxy(:,1) = Mz0 .* sin(flips(:,1));
    # Mz(:,1)  = Mz0 .* cos(flips(:,1));
    Mxy[:, 0] = Mz0 * np.sin(flips[:, 0])
    Mz[:, 0]  = Mz0 * np.cos(flips[:, 0])

    # Loop over flips
    for n in range(1, N):
        # Evolve from previous step
        # if use_input_function:
        #     Mz_m = Ad*Mz(:,n-1) + Bd(:,1)*input_function(n-1);
        # else:
        #     Mz_m = Ad*Mz(:,n-1);
        Mz_prev = Mz[:, n - 1]
        if use_input_function:
            # Only site 1 gets the input; "Bd[:,0]" picks the first column if needed
            Mz_m = Ad @ Mz_prev + Bd[:, 0] * input_function[n - 1]
        else:
            Mz_m = Ad @ Mz_prev

        # Apply flip n
        Mxy[:, n] = Mz_m * np.sin(flips[:, n])
        Mz[:, n]  = Mz_m * np.cos(flips[:, n])

    return Mxy, Mz




if __name__ == "__main__":
    # Example usage


    # 6 flip points
    N = 6
    TR = 1.0

    # 2 sites, 6 flips each
    flips = np.array([
        [np.deg2rad(10), np.deg2rad(20), np.deg2rad(30), np.deg2rad(40), np.deg2rad(50), np.deg2rad(60)],
        [np.deg2rad(5),  np.deg2rad(10), np.deg2rad(15), np.deg2rad(15), np.deg2rad(20), np.deg2rad(25)]
    ])

    # Relaxation rates
    R1 = [1.0, 1.2]

    # Exchange rates [forward, reverse]
    k = [0.05, 0.02]

    # Start from 3 seconds evolution from [1,0]
    Tin = 3.0

    # No input function
    input_function = np.zeros(N)

    Mxy, Mz = simulate_2site_model(Tin, R1, k, flips, TR, input_function)

    print("Mxy =\n", Mxy)
    print("Mz =\n", Mz)

    # Simple plots
    time_points = np.arange(N) * TR

    plt.figure()
    for site_idx in range(2):
        plt.plot(time_points, Mz[site_idx, :], marker='o', label=f"Site {site_idx+1}")
    plt.xlabel("Time (s)")
    plt.ylabel("Mz After Flip")
    plt.title("Mz for Each Site (2-Site Model)")
    plt.legend()
    plt.savefig("Mz2")

    plt.figure()
    for site_idx in range(2):
        plt.plot(time_points, Mxy[site_idx, :], marker='o', label=f"Site {site_idx+1}")
    plt.xlabel("Time (s)")
    plt.ylabel("Mxy After Flip")
    plt.title("Mxy for Each Site (2-Site Model)")
    plt.legend()
    plt.savefig("Mxy2")
