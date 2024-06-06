import numpy as np
from scipy.optimize import least_squares
import matplotlib.pyplot as plt
import time
from scipy.linalg import expm
import tqdm


def py_fit_pyr_kinetics(S, TR, flips, params_fixed, params_est, plot_flag):
    """ 
    Python implementation of fit_pyr_kinetics.m 
    
    Fits precursor-product kinetic model, assuming origination from a single substrate/precursor to multiple products
    An "input-less" method is used, eliminating need to make any assumptions about the input function.

    INPUTS: 
    S - multi-voxel or one voxel data over time [voxels, # of metabolites, # of time points]
        substrate signal followed by metabolites (assumes order: pyr, lac, bic, ala)
    TR - temporal resolution/ time of repetition in sec [scalar]
    flips - flip angle in degrees for each metabolite per time pt [# of mets, # of timepts]
    params_fixed - dict of fixed params and their values
    params_fit: dict of fit params and their values
    plot_flag  - bool, 1 will plot the fits and print results
    
    OUTPUTS: 
    Sfit - fit metabolite time courses [voxel, # of metabolites - 1, # of time points]
           Note: pyruvate fit is not retuned as pyruvate signal is used for input estimation
    ufit - fit input function time course [voxel, # of time points]
    params_fit - vector of fit parameters [voxel, # of params to fit]
    Rsq - Rsq error for each voxel per metabolite [voxel, # of metabolites-1]
    CHIsq - chi sq error for each voxel per metabolite [voxel, # of metabolites-1]
    NRMSE - normalized root mean square error [voxel, # of metabolites-1]

    Translated to Python by Sule Sahin
    Copyright 2024 UCSF
    """

    # extract and define data dimensions
    size_S = S.shape
    ndimsx = len(size_S) - 2
    Nt = size_S[-1]
    t = np.arange(Nt) * TR
    Nx = size_S[0:ndimsx]
    prodNx = int(np.prod(Nx))
    Nmets = size_S[-2]
    if not Nx:
        Nx = 1

    # set parameter defaults
    params_all = np.array(['kPL', 'kPB', 'kPA', 'R1P', 'R1L', 'R1A', 'R1B', 'S0_P', 'S0_L', 'S0_B', 'S0_A'],
                          dtype=object)
    params_default_est = {'kPL': 0.01, 'kPB': 0.01, 'kPA': 0.01, 'R1P': 1 / 30, 'R1L': 1 / 25, 'R1A': 1 / 25,
                          'R1B': 1 / 15, 'S0_P': 0, 'S0_L': 0, 'S0_B': 0, 'S0_A': 0}
    params_default_lb = {'kPL': -np.inf, 'kPB': -np.inf, 'kPA': -np.inf, 'R1P': 1 / 50, 'R1L': 1 / 50, 'R1A': 1 / 50,
                         'R1B': 1 / 50, 'S0_P': -np.inf, 'S0_L': -np.inf, 'S0_B': -np.inf, 'S0_A': -np.inf}
    params_default_ub = {'kPL': np.inf, 'kPB': np.inf, 'kPA': np.inf, 'R1P': 1 / 10, 'R1L': 1 / 10, 'R1A': 1 / 10,
                         'R1B': 1 / 5, 'S0_P': np.inf, 'S0_L': np.inf, 'S0_B': np.inf, 'S0_A': np.inf}

    # setup fitting based on number of metabolites
    products_string = ['lactate', 'bicarbonate', 'alanine']
    if Nmets == 2:  # assume pyruvate & lactate
        params_fixed['kPA'] = 0
        params_fixed['S0_A'] = 0
        params_fixed['R1A'] = 1
        params_fixed['kPB'] = 0
        params_fixed['S0_B'] = 0
        params_fixed['R1B'] = 1
        products_string = ['lactate']
    elif Nmets == 3:  # assume pyruvate & lactate & bicarbonate
        params_fixed['kPA'] = 0
        params_fixed['S0_A'] = 0
        params_fixed['R1A'] = 1
        products_string = ['lactate', 'bicarbonate']

    # determine which parameters need to be fit
    l_params_est = np.zeros([len(params_all) - len(params_fixed)], dtype=object)
    nt = 0
    for param in params_all:
        if param not in params_fixed:
            l_params_est[nt] = param
            nt = nt + 1
    Nparams_to_fit = l_params_est.shape[0]

    # set estimated or fixed values for parameters
    params_est_vec = np.zeros([Nparams_to_fit])
    params_lb = np.zeros([Nparams_to_fit])
    params_ub = np.zeros([Nparams_to_fit])
    for i in range(Nparams_to_fit):
        if l_params_est[i] in params_est:
            params_est_vec[i] = params_est[l_params_est[i]]
        else:
            params_est_vec[i] = params_default_est[l_params_est[i]]
        if l_params_est[i] + '_lb' in params_est:
            params_lb[i] = params_est[l_params_est[i] + '_lb']
        else:
            params_lb[i] = params_default_lb[l_params_est[i]]
        if l_params_est[i] + '_ub' in params_est:
            params_ub[i] = params_est[l_params_est[i] + '_ub']
        else:
            params_ub[i] = params_default_ub[l_params_est[i]]

    # reshape data
    Sreshape = np.reshape(S, [prodNx, Nmets, Nt])
    if Nmets < 4:
        Sreshape = np.concatenate((Sreshape, np.zeros([prodNx, 4 - Nmets, Nt])), axis=1)
        flips = np.concatenate((flips, np.ones([4 - Nmets, flips.shape[1]])), axis=0)

    [Sscale, Mzscale] = flips_scaling_factors(flips, Nt)  # extract scaling for longitudinal/transverse component

    # define arrays to fill
    params_fit_vec = np.zeros([np.prod(Nx), Nparams_to_fit])
    objective_val = np.zeros(np.prod(Nx))
    CHIsq = np.zeros([np.prod(Nx), Nmets - 1])
    Rsq = np.zeros([np.prod(Nx), Nmets - 1])
    Sfit = np.zeros([np.prod(Nx), Nmets - 1, Nt])
    ufit = np.zeros([np.prod(Nx), Nt])
    NRMSE = np.zeros([np.prod(Nx), Nmets-1])
    for i in tqdm.tqdm(range(Sreshape.shape[0])):  # loop over each voxel
        Mxy = np.reshape(Sreshape[i, :, :], [4, Nt])

        if Mxy.any():
            Mz = Mxy / Sscale
            Istart = 0  # option to propogate model from various pts in time

            # nonlinear least square fitting using difference_inputless as obj fxn
            res = least_squares(difference_inputless, params_est_vec, bounds=(params_lb, params_ub), ftol=1e-06,
                                xtol=1e-06, gtol=1e-06, args=(params_fixed, TR, Mzscale, Sscale, Mz, Istart, Nmets))
            params_fit_vec[i, :] = res.x
            objective_val[i] = res.cost

            # get fit curves using fit parameters
            (Mzfit, ufit[i, :]) = trajectories_inputless(params_fit_vec[i, :], params_fixed, TR, Mzscale, Mz[0, :],
                                                         Istart)
            Sfit[i, :, :] = Mzfit[1:Nmets, :] * Sscale[1:Nmets, :]
            ufit[i, :] = ufit[i, :] * Sscale[0, :]

            # calculate errors
            Rsq_denom = np.sum((Sreshape[i, 1:Nmets, :] - np.repeat(np.reshape(np.mean(Sreshape[i, 1:Nmets, :], axis=1),
                                                                               [Nmets-1, 1]), Nt, axis=1)) ** 2, axis=1)
            Rsq[i, :] = 1 - np.sum((Sreshape[i, 1:Nmets, :] - Sfit[i, :, :]) ** 2, axis=1) / Rsq_denom
            CHIsq[i, :] = np.sum((Sreshape[i, 1:Nmets, :] - Sfit[i, :, :]) ** 2 / (Sreshape[i, 1:Nmets, :] + Sfit[i, :, :]), axis=1) / 2
            NRMSE[i, :] = nrmse(Sfit[i, :, :], Sreshape[i, 1:Nmets, :])

            if plot_flag:  # plot fits and print parameters if plot flag
                fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10,10))
                ax1.plot(t, Mz[0:Nmets, :].T, t, Mzfit[1:Nmets, :].T, '--', t, ufit[i, :] / Sscale[0, :], 'k')
                ax1.set_xlabel('time(s)')
                ax1.set_ylabel('state magnetization (au)')
                ax1.set_xlim(0, TR * Nt)

                ax2.plot(t, Mxy[0:Nmets, :].T, t, Sfit[i, 0:Nmets - 1, :].T, '--', t, ufit[i, :], 'k')
                ax2.set_xlabel('time(s)')
                ax2.set_ylabel('signal (au)')
                ax2.set_xlim(0, TR * Nt)

                fit_results_string = ''
                for n in range(Nparams_to_fit):
                    fit_results_string = fit_results_string + ' ' + l_params_est[n] + '=' + str(
                        np.round(params_fit_vec[i, n], decimals=4))
                plt.title(fit_results_string)
                print(fit_results_string)

                products_legend = ['pyruvate'] + products_string + [p + ' fit' for p in products_string] + ['input '
                                                                                                            'estimate']
                plt.legend(products_legend)
                plt.show()
                time.sleep(0.5)

    # reshape fitted curves back into original shape
    if ndimsx > 0:
        Sfit = np.reshape(Sfit, Nx + (Nmets-1, Nt))
        ufit = np.reshape(ufit, Nx + (Nt,))

    return Sfit, ufit, params_fit_vec, objective_val, Rsq, CHIsq, NRMSE


def flips_scaling_factors(flips, Nt):
    # calculates scale for transverse and longitudinal components using flip angle for each metabolite
    # allows for multiple excitations per time point

    Nflips = int(flips.shape[1] / Nt)
    Nmets = flips.shape[0]

    Mzscale = np.zeros([Nmets, Nt])
    Sscale = np.zeros([Nmets, Nt])
    for t in range(Nt):
        Iflips = np.arange(Nflips) + (t * Nflips)
        Mzscale[:, t] = np.prod(np.cos(flips[:, Iflips]), axis=1)
        for n in range(Nflips):
            Sscale[:, t] = Sscale[:, t] + (np.sin(flips[:, Iflips[n]]) * np.prod(np.cos(flips[:, Iflips[0:n]]), axis=1))
    Sscale = Sscale / Nflips
    return Sscale, Mzscale


def difference_inputless(params_fit, params_fixed, TR, Mzscale, Sscale, Mz, Istart, Nmets):
    # calculates difference between data and fit (objective function)
    [Mzfit, ufit] = trajectories_inputless(params_fit, params_fixed, TR, Mzscale, Mz[0, :], Istart)
    temp_diff = (Mz - Mzfit) * Sscale
    diff_products = temp_diff[1:Nmets, :]
    diff_products = diff_products.flatten()
    return diff_products


def trajectories_inputless(params_fit, params_fixed, TR, Mzscale, Mz_pyr, Istart):
    # Compute product magnetizations using a uni-directional two-site model
    # Uses substrate magnetization measurements, estimated relaxation and conversion rates
    Nmets = Mzscale.shape[0]
    N = Mzscale.shape[1]
    u = np.zeros(N)

    # extract fixed and fit params and assign them to a dict: params_all
    params_all = np.array(['kPL', 'kPB', 'kPA', 'R1P', 'R1L', 'R1A', 'R1B', 'S0_P', 'S0_L', 'S0_B', 'S0_A'],
                          dtype=object)
    params_full = {p: 0 for p in params_all}
    l_params_est = np.zeros([len(params_all) - len(params_fixed)], dtype=object)
    nt = 0
    for param in params_all:
        if param not in params_fixed:
            l_params_est[nt] = param
            nt = nt + 1

    # assign fixed/fit values ot each param
    for param in params_all:
        if param in params_fixed:
            params_full[param] = params_fixed[param]
        else:
            params_full[param] = params_fit[list(l_params_est).index(param)]

    # initialize Mz
    Mz_all = np.zeros([Nmets, N])
    Mz_all[0, :] = Mz_pyr
    Mz_all[1, Istart] = params_full['S0_L']
    Mz_all[2, Istart] = params_full['S0_B']
    Mz_all[3, Istart] = params_full['S0_A']

    A = np.array([[-params_full['R1P'] - params_full['kPL'] - params_full['kPB'] - params_full['kPA'], 0, 0, 0],
                  [params_full['kPL'], -params_full['R1L'], 0, 0],
                  [params_full['kPB'], 0, -params_full['R1B'], 0],
                  [params_full['kPA'], 0, 0, -params_full['R1A']]])  # define state matrix

    for It in range(Istart, N - 1):  # model forward in time
        # model equations based on solution of differential equations for the metabolite exchange
        Mz_init = Mz_all[:, [It]] * Mzscale[:, [It]]

        exp_term = np.exp((-params_full['R1P'] - params_full['kPL'] - params_full['kPB'] - params_full['kPA']) * TR)
        add_term = params_full['R1P'] + params_full['kPL'] + params_full['kPB'] + params_full['kPA']
        u[It] = (Mz_pyr[It + 1] - Mz_init[0] * exp_term) * add_term / (1 - exp_term)

        xstar = - np.linalg.inv(A) @ np.array([[u[It]], [0], [0], [0]])
        Mz_all[:, [It + 1]] = xstar + expm(A * TR) @ (Mz_init - xstar)

    for It in range(Istart, 0, -1):  # model backwards in time
        Mz_init = Mz_all[:, It]

        exp_term = np.exp((- params_full['R1P'] - params_full['kPL'] - params_full['kPB'] - params_full['kPA']) * -TR)
        add_term = params_full['R1P'] + params_full['kPL'] + params_full['kPB'] + params_full['kPA']
        u[It - 1] = (Mz_pyr[0, It - 1] * Mzscale[0, It - 1] - Mz_init[0, 0] * exp_term) * add_term / (1 - exp_term)

        xstar = - np.linalg.inv(A) @ np.array([[u[It - 1]], [0], [0], [0]])
        Mz_plus = xstar + expm(A * -TR) @ (Mz_init - xstar)
        Mz_all[:, It - 1] = Mz_plus / Mzscale[:, It - 1]

    return Mz_all, u


def nrmse(fit, real):
    rmse = np.sqrt(np.mean(((fit - real) ** 2), axis=1))
    mean = np.mean(fit, axis=1)
    return rmse/mean