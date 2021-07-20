#!/usr/bin/env python3
# python version 3.6.14 +


# Standard modules
import math
import pickle

# Third-party modules
import numpy as np

# Local modules
import diffusion3d_virus_town_large as dv


def load_init_terms(gConf, hosts):
    s = np.zeros((gConf['nx'], gConf['ny'], gConf['nz'], gConf['ng']))
    T = np.zeros((gConf['nx'], gConf['ny'], gConf['nz'], gConf['ng']))

    for h in hosts:
        s[h['x'], h['y'], h['z'], [1, 5]] = h['lam']
        T[h['x'], h['y'], h['z'], :] = (
            h['T1_0_U'], h['T2_0_U'], h['I_0_U'], h['V_0_U'],
            h['T1_0_L'], h['T2_0_L'], h['I_0_L'], h['V_0_L'],
            gConf['env_0']
        )

    return s, T


def update_A(A_diag, A_off_diag, hosts, T, t, gConf):
    def delta(delta_I, sigma, mu, t):
        return delta_I * np.e ** (sigma * (t - mu)) if (t >= mu) else delta_I

    for h in hosts:
        A_diag[h['x'], h['y'], h['z'], :] = (
            # T1_U, T2_U
            h['beta_U'] * T[h['x'], h['y'], h['z'], 3],
            h['beta_U'] * T[h['x'], h['y'], h['z'], 3],
            # I_U
            delta(h['delta_I'], h['sigma_U'], h['mu'], t)
            + h['w_U'] * T[h['x'], h['y'], h['z'], 1],
            # V_U
            h['c_U'],
            # T1_L, T2_L
            h['beta_L'] * T[h['x'], h['y'], h['z'], 7],
            h['beta_L'] * T[h['x'], h['y'], h['z'], 7],
            # I_L
            delta(h['delta_I'], h['sigma_L'], h['mu'], t)
            + h['w_L'] * T[h['x'], h['y'], h['z'], 5],
            # V_L
            h['c_L'],
            # V_env
            gConf['halflife'] + gConf['ventilation']
        )

        A_off_diag[h['x'], h['y'], h['z'],
                   [2, 2, 3, 3, 6, 6, 7, 7, 8, 8],
                   [0, 1, 2, 8, 4, 5, 3, 6, 3, 7]
                   ] = (
            -h['beta_U'] * T[h['x'], h['y'], h['z'], 3],
            -h['beta_U'] * T[h['x'], h['y'], h['z'], 3],
            -h['p_U'], -h['a_inhale'],
            -h['beta_L'] * T[h['x'], h['y'], h['z'], 7],
            -h['beta_L'] * T[h['x'], h['y'], h['z'], 7],
            -h['a_conduct'], -h['p_L'],
            -gConf['gamma_U'], -gConf['gamma_L']
        )


def inhost_viral_model_solver(conf, gConf, hosts):
    """ An inhost viral model solver using `sim_time_steping_diffusion_calc`
    to solve one or many within-host viral dynamic models simultaneously

    input
    =====
    gConf: dict, values about the grid
        nx, ny, nz, ng, ng2
        error_solv, error_solv_ng, i_upwind, i_harmonic,
        ndim_vel, nits, nits_solv_ng, relax
        dx, dy, dz, dt, ntime

    hostCells: dict, info about each host
        position - x, y, z
        parameters - T1_0, T2_0, I_0, V_0, beta, lam, delta, w, p, c
    """
    # initialise matrices
    #
    v_n = np.ones(gConf['ng'])
    # advection-diffusion in the system
    kdiff = np.zeros((gConf['nx'], gConf['ny'], gConf['nz'], gConf['ng']))
    kdiff[:, :, :, -1] = 1.e10
    # advection velocity
    u = np.zeros((conf['ndim_vel'], gConf['nx']+1, gConf['ny']+1,
                  gConf['nz']+1, gConf['ng']))
    #
    A_diag = np.zeros((gConf['nx'], gConf['ny'], gConf['nz'], gConf['ng']))
    A_off_diag = np.zeros((gConf['nx'], gConf['ny'], gConf['nz'],
                           gConf['ng'], gConf['ng2']))
    # source term and initial values
    s, T = load_init_terms(gConf, hosts)

    time_list = [0]
    T_list = [T]
    for _ in range(conf['nstep']):
        t = time_list[-1] + conf['dt'] * conf['ntime']
        update_A(A_diag, A_off_diag, hosts, T, t, gConf)

        T = dv.sim_time_steping_diffusion_calc(
            T, conf['ntime'], conf['nits'], conf['nits_solv_ng'],
            conf['relax'], gConf['error_solv'], gConf['error_solv_ng'],
            gConf['dx'], gConf['dy'], gConf['dz'], conf['dt'],
            v_n, A_diag, A_off_diag, kdiff, s, u,
            gConf['i_upwind'], gConf['i_harmonic'],
            conf['ndim_vel'], gConf['nx'], gConf['ny'], gConf['nz'],
            gConf['ng'], gConf['ng2']
        )

        time_list.append(t)
        T_list.append(T)

    return np.array(T_list), time_list


if __name__ == '__main__':
    # global constants (Wang, S. et al. 2020)
    T1_0 = 6.e4
    T2_0 = 0
    I_0 = 0
    lam = 1.e4
    delta_I = 2

    # hyper-param
    const_sec_per_day = 86400
    d_sec = 60  # amount of second for each time step
    t_sec = 30 * const_sec_per_day  # amount of second for running the model
    ntime = 1  # output every n step
    dt = d_sec / const_sec_per_day  # convert time step unit to day
    nstep = math.ceil(t_sec / d_sec / ntime)  # amount of steps

    conf = {
        # dt, internal and external iterations
        'dt': dt, 'ntime': ntime, 'nstep': nstep,
        # solve a 2D problem
        'ndim_vel': 2,
        # non-linear iterations
        'nits': 3, 'nits_solv_ng': 3,
        # relaxation value
        'relax': 1.0,
    }
    gConf = {
        # number of cells in the simulation
        'nx': 5, 'ny': 5, 'nz': 3,
        # number of equations (energy group)
        'ng': 9, 'ng2': 9,
        # error tolerances
        'error_solv': 1.e-6, 'error_solv_ng': 1.e-6,
        # upwind differencing, and harmonic averaging for the diffusion coeff
        'i_upwind': 1, 'i_harmonic': 0,
        # spatial discretization
        'dx': 0.1, 'dy': 0.1, 'dz': 0.1,
        #
        'env_0': 0, 'halflife': 0.03125, 'ventilation': 0.,
        #
        'gamma_U': 1.e-2, 'gamma_L': 1.e-2
    }
    hosts = [
        # patient #7 URT & LRT (Wang, S. et al. 2020)
        {
            'x': 1, 'y': 1, 'z': 1,
            'T1_0_U': T1_0, 'T2_0_U': T2_0, 'I_0_U': I_0, 'V_0_U': 1.e-4,
            'T1_0_L': T1_0, 'T2_0_L': T2_0, 'I_0_L': I_0, 'V_0_L': 0,
            'lam': lam, 'delta_I': delta_I,
            'beta_U': 9.9e-7, 'p_U': 1.08e4, 'c_U': 48, 'w_U': 2.1e-4,
            'sigma_U': 1.e-4,
            'beta_L': 1.e-6, 'p_L': 1.1e5, 'c_L': 209, 'w_L': 4.5e-4,
            'sigma_L': 0.11, 'mu': 9,
            'a_conduct': 2 ** -6, 'a_inhale': 0.
        },
        # patient #7 URT & LRT (Wang, S. et al. 2020)
        {
            'x': 1, 'y': 3, 'z': 1,
            'T1_0_U': T1_0, 'T2_0_U': T2_0, 'I_0_U': I_0, 'V_0_U': 0,
            'T1_0_L': T1_0, 'T2_0_L': T2_0, 'I_0_L': I_0, 'V_0_L': 0,
            'lam': lam, 'delta_I': delta_I,
            'beta_U': 9.9e-7, 'p_U': 1.08e4, 'c_U': 48, 'w_U': 2.1e-4,
            'sigma_U': 1.e-4,
            'beta_L': 1.e-6, 'p_L': 1.1e5, 'c_L': 209, 'w_L': 4.5e-4,
            'sigma_L': 0.11, 'mu': 9,
            'a_conduct': 2 ** -6, 'a_inhale': 1.
        },
        # patient #7 URT & LRT (Wang, S. et al. 2020)
        {
            'x': 3, 'y': 1, 'z': 1,
            'T1_0_U': T1_0, 'T2_0_U': T2_0, 'I_0_U': I_0, 'V_0_U': 0,
            'T1_0_L': T1_0, 'T2_0_L': T2_0, 'I_0_L': I_0, 'V_0_L': 0,
            'lam': lam, 'delta_I': delta_I,
            'beta_U': 9.9e-7, 'p_U': 1.08e4, 'c_U': 48, 'w_U': 2.1e-4,
            'sigma_U': 1.e-4,
            'beta_L': 1.e-6, 'p_L': 1.1e5, 'c_L': 209, 'w_L': 4.5e-4,
            'sigma_L': 0.11, 'mu': 9,
            'a_conduct': 2 ** -6, 'a_inhale': .5
        },
        # patient #7 URT & LRT (Wang, S. et al. 2020)
        {
            'x': 3, 'y': 3, 'z': 1,
            'T1_0_U': T1_0, 'T2_0_U': T2_0, 'I_0_U': I_0, 'V_0_U': 0,
            'T1_0_L': T1_0, 'T2_0_L': T2_0, 'I_0_L': I_0, 'V_0_L': 0,
            'lam': lam, 'delta_I': delta_I,
            'beta_U': 9.9e-7, 'p_U': 1.08e4, 'c_U': 48, 'w_U': 2.1e-4,
            'sigma_U': 1.e-4,
            'beta_L': 1.e-6, 'p_L': 1.1e5, 'c_L': 209, 'w_L': 4.5e-4,
            'sigma_L': 0.11, 'mu': 9,
            'a_conduct': 2 ** -6, 'a_inhale': 1.e-3
        },
    ]

    URTs = []
    LRTs = []
    T_list, t_list = inhost_viral_model_solver(conf, gConf, hosts)
    URTs.append(T_list[:, 1, 1, 1, 3])
    URTs.append(T_list[:, 1, 3, 1, 3])
    URTs.append(T_list[:, 3, 1, 1, 3])
    URTs.append(T_list[:, 3, 3, 1, 3])
    LRTs.append(T_list[:, 1, 1, 1, 7])
    LRTs.append(T_list[:, 1, 3, 1, 7])
    LRTs.append(T_list[:, 3, 1, 1, 7])
    LRTs.append(T_list[:, 3, 3, 1, 7])

    saveDict = {
        'URTs': URTs,
        'LRTs': LRTs,
        't_list': t_list,
    }
    with open('tmp.pkl', 'wb') as f:
        pickle.dump(saveDict, f, pickle.HIGHEST_PROTOCOL)
