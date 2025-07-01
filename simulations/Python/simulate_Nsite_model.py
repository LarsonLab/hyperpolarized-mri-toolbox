import numpy as np
from scipy.linalg import expm
import matplotlib.pyplot as plt


def simulate_Nsite_model(Tin, R1, k, flips, TR, input_function=None):
    """
    Simulates the magnetization evolution in an N-site exchange model 
    with a single substrate and multiple products for varying flip angles.
    
    Parameters
    ----------
    Tin : float or array_like
        If Tin is a single scalar, it is interpreted as a time (in seconds) 
        during which the system evolves from all magnetization in the first site.
        If Tin is an array of length Nmet, it is the initial magnetization vector.
    R1 : float or array_like
        Longitudinal relaxation rates (1/s) for each site. 
        If a single float is provided, it is broadcasted to all sites.
    k : array_like
        A matrix of shape (Nmet, 1) or (Nmet, 2). The first column 
        contains forward exchange rates (1/s). The optional second column 
        contains reverse exchange rates (1/s).
    flips : ndarray
        A 2D array of shape (Nmet, N) specifying the flip angles (in radians) 
        for each metabolite across N time points.
    TR : float
        Repetition time in seconds (time between flips).
    input_function : array_like or None, optional
        External input added to the first site at each time point. 
        If None or all zeros, no input is added. 
        Otherwise, we assume `input_function[n]` is added over the nth interval.
    
    Returns
    -------
    Mxy : ndarray
        A 2D array of shape (Nmet, N) containing the transverse magnetization 
        immediately after each flip.
    Mz : ndarray
        A 2D array of shape (Nmet, N) containing the longitudinal magnetization 
        immediately after each flip.
    """
    
    flips = np.asarray(flips, dtype=float)
    # Number of metabolites
    Nmet = flips.shape[0]
    # Number of time points
    N = flips.shape[1]
    
    # Ensure R1 is a numpy array of length Nmet
    if np.isscalar(R1):
        R1 = np.full(Nmet, R1, dtype=float)
    else:
        R1 = np.asarray(R1, dtype=float)
        if R1.size != Nmet:
            raise ValueError("R1 must be either a scalar or an array of length Nmet.")
    
    # k should be of shape (Nmet, 2) at the end
    k = np.asarray(k, dtype=float)
    if k.ndim == 1:
        # If only forward rates provided, assume no reverse rates
        if k.size == Nmet:
            k = np.column_stack((k, np.zeros(Nmet)))
        else:
            raise ValueError("k must have shape (Nmet,) or (Nmet,2).")
    elif k.shape[1] == 1:
        # If only one column, add a column of zeros for reverse
        k = np.column_stack((k, np.zeros(Nmet)))
    
    # Check input_function
    if input_function is None:
        input_function = np.zeros(N, dtype=float)
    else:
        input_function = np.asarray(input_function, dtype=float)
        if input_function.size != N and not np.allclose(input_function, 0):
            raise ValueError("input_function must be None, all zeros, or length N.")
    
    use_input_function = not np.allclose(input_function, 0)
    Nsim = 100  # number of subdivisions if using input
    
    # Build A matrix depending on Nmet
    # The structure is:
    #   A[i,i]   = - R1[i] - sum_of_forward_rates_from_i
    #   A[i,j]   = + some_exchange_rate
    # For 2, 3, 4 sites, we follow the pattern in the original code.
    def build_A(Nmet, R1, k):
        if Nmet == 2:
            A = np.array([
                [-R1[0] - k[0,0],       +k[0,1]     ],
                [      +k[0,0],  -R1[1] - k[0,1]    ]
            ], dtype=float)
        
        elif Nmet == 3:
            A = np.array([
                [-R1[0] - k[0,0] - k[1,0], +k[0,1],           +k[1,1]          ],
                [         +k[0,0],        -R1[1] - k[0,1],    0                ],
                [         +k[1,0],        0,                 -R1[2] - k[1,1]  ]
            ], dtype=float)
        
        elif Nmet == 4:
            A = np.array([
                [-R1[0] - k[0,0] - k[1,0], +k[0,1],           +k[1,1],          +k[2,1]],
                [           +k[0,0],      -R1[1] - k[0,1],    0,                0      ],
                [           +k[1,0],      0,                 -R1[2] - k[1,1],  0      ],
                [           +k[2,0],      0,                 0,                -R1[3] - k[2,1]]
            ], dtype=float)
        else:
            raise ValueError("This code handles up to 4-site models (Nmet <= 4).")
        
        return A
    
    A = build_A(Nmet, R1, k)
    
    # Matrix exponential for full TR
    Ad_TR = expm(A * TR)
    if use_input_function:
        # Matrix exponential for TR subdivided into smaller segments
        Ad_Nsim = expm(A * (TR / Nsim))
    else:
        Ad_Nsim = None
    
    # Compute initial Mz0
    # If Tin is a single value, interpret it as an evolution time from
    # an initial state of [1, 0, 0, ...].
    # Otherwise, interpret Tin as the initial magnetization vector.
    if np.isscalar(Tin):
        Tin = float(Tin)
        # Start with Mz = [1, 0, 0, ...]
        Mz_init = np.zeros(Nmet, dtype=float)
        Mz_init[0] = 1.0
        Mz0 = expm(A * Tin).dot(Mz_init)
    else:
        Tin = np.asarray(Tin, dtype=float)
        if Tin.size != Nmet:
            raise ValueError("Tin must be a scalar or an array of length Nmet.")
        Mz0 = Tin
    
    # Prepare output arrays
    Mxy = np.zeros((Nmet, N), dtype=float)
    Mz  = np.zeros((Nmet, N), dtype=float)
    
    # Immediately after the first flip
    Mxy[:, 0] = Mz0 * np.sin(flips[:, 0])
    Mz[:, 0]  = Mz0 * np.cos(flips[:, 0])
    
    # Loop over subsequent flips
    for n in range(1, N):
        Mz_prev = Mz[:, n - 1]
        if use_input_function:
            # Subdivide TR into small steps, each step adding a fraction of input
            for _ in range(Nsim):
                Mz_prev = Ad_Nsim @ (Mz_prev + np.r_[[input_function[n - 1] / Nsim], 
                                                     np.zeros(Nmet - 1, dtype=float)])
            Mz_m = Mz_prev
        else:
            # Evolve magnetization over TR with no additional input
            Mz_m = Ad_TR @ Mz_prev
        
        # Apply flip
        Mxy[:, n] = Mz_m * np.sin(flips[:, n])
        Mz[:, n]  = Mz_m * np.cos(flips[:, n])
    
    return Mxy, Mz

print("hi")
N = 4  # number of time points (flips)
TR = 1.0  # repetition time, in seconds

# b) Flip angles for each metabolite across N time points
#    Here, we have 2 sites (rows = 2) and 6 time points (columns = 6).
flips = np.array([
    [np.deg2rad(10), np.deg2rad(20), np.deg2rad(30), np.deg2rad(40), np.deg2rad(50), np.deg2rad(60)],
    [np.deg2rad(5),  np.deg2rad(10), np.deg2rad(15), np.deg2rad(15), np.deg2rad(20), np.deg2rad(25)]
], dtype=float)

# c) Relaxation rates (R1) for each site
#    If you give a single scalar, it will be broadcast to both sites.
#    Otherwise, give a length-2 array like [1.0, 0.8].
R1 = [1.0, 0.8]  # 1/s

# d) Exchange rates, shape (Nmet,2) = (2,2) for forward & reverse
#    Suppose site1 -> site2 has k(1,1)=0.05, site2 -> site1 has k(1,2)=0.02
k = [0.05, 0.02]  
# site1 forward & reverse # site2 doesn't convert to a "third site" in 2-site model
# Actually, in a pure 2-site scenario, you'd only need the first row for k(1,*).
# We place zeros in the second row to match the shape from the original code logic.

# e) Initial magnetization ("Tin")
#    If scalar => evolve from [1, 0] for 'Tin' seconds before the first flip
#    If array => direct initial Mz
Tin = 3.0  # Evolve 3 seconds from [1, 0] before first flip

# f) Input function
#    Let's set it to all zeros for no external infusion
input_function = np.zeros(N)  # shape = (N,)

# 3) Run the simulation
Mxy, Mz = simulate_Nsite_model(
    Tin=Tin,
    R1=R1,
    k=k,
    flips=flips,
    TR=TR,
    input_function=input_function
)

# 4) Inspect or print results
print("Mxy =\n", Mxy)
print("Mz =\n", Mz)

# 5) (Optional) Plot the results
time_points = np.arange(N) * TR  # times at which flips occur

#Plot Mz
plt.figure()
for site_idx in range(Mz.shape[0]):
    plt.plot(time_points, Mz[site_idx, :], marker='o', label=f"Site {site_idx+1}")
plt.xlabel("Time (s)")
plt.ylabel("Mz After Flip")
plt.title("Mz for Each Site After Each Flip")
plt.legend()
plt.savefig("Mz.png")
# Plot Mxy
plt.figure()
for site_idx in range(Mxy.shape[0]):
    plt.plot(time_points, Mxy[site_idx, :], marker='o', label=f"Site {site_idx+1}")
plt.xlabel("Time (s)")
plt.ylabel("Mxy After Flip")
plt.title("Mxy for Each Site After Each Flip")
plt.legend()
plt.savefig("Mxy.png")