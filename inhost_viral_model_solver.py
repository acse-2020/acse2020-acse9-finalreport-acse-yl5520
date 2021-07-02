# Standard modules
import math

# Third-party modules
import matplotlib.pyplot as plt
import numpy as np

# Local modules
import diffusion3d_virus_town_large as dv


def load_init_terms(gConf, hosts):
    s = np.zeros((gConf['nx'], gConf['ny'], gConf['nz'], gConf['ng']))
    T = np.zeros((gConf['nx'], gConf['ny'], gConf['nz'], gConf['ng']))

    for h in hosts:
        s[h['x'], h['y'], h['z'], 1] = h['lam']
        T[h['x'], h['y'], h['z'], :] = (h['T1_0'], h['T2_0'],
                                        h['I_0'], h['V_0'])

    return s, T


def update_A(A_diag, A_off_diag, hosts, T, t):
    def delta(delta_I, sigma, mu, t):
        return delta_I * np.e ** (sigma * (t - mu)) if (t >= mu) else delta_I

    for h in hosts:
        A_diag[h['x'], h['y'], h['z'], :] = (
            # T1, T2
            h['beta'] * T[h['x'], h['y'], h['z'], 3],
            h['beta'] * T[h['x'], h['y'], h['z'], 3],
            # I
            delta(h['delta_I'], h['sigma'], h['mu'], t)
            + h['w'] * T[h['x'], h['y'], h['z'], 1],
            # V
            h['c']
        )

        A_off_diag[h['x'], h['y'], h['z'], [2, 2, 3], [0, 1, 2]] = (
            -h['beta'] * T[h['x'], h['y'], h['z'], 3],
            -h['beta'] * T[h['x'], h['y'], h['z'], 3],
            -h['p']
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
    # advection velocity
    v_n = np.ones(gConf['ng'])
    # advection-diffusion in the system
    kdiff = np.zeros((gConf['nx'], gConf['ny'], gConf['nz'], gConf['ng']))
    #
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
        update_A(A_diag, A_off_diag, hosts, T, t)

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
        'nx': 5, 'ny': 3, 'nz': 3,
        # number of equations (energy group)
        'ng': 4, 'ng2': 4,
        # error tolerances
        'error_solv': 1.e-6, 'error_solv_ng': 1.e-6,
        # upwind differencing, and harmonic averaging for the diffusion coeff
        'i_upwind': 1, 'i_harmonic': 0,
        # spatial discretization
        'dx': 0.1, 'dy': 0.1, 'dz': 0.1,
    }
    hosts = [
        # patient #7 LRT (Wang, S. et al. 2020)
        {
            'x': 1, 'y': 1, 'z': 1,
            'T1_0': T1_0, 'T2_0': T2_0, 'I_0': I_0,
            'lam': lam, 'delta_I': delta_I,
            'beta': 1.e-6, 'p': 1.1e5, 'c': 209, 'w': 4.5e-4,
            'sigma': 0.11, 'V_0': 1.e-3, 'mu': 9,
        },
        # patient #14 LRT (Wang, S. et al. 2020)
        {
            'x': 3, 'y': 1, 'z': 1,
            'T1_0': T1_0, 'T2_0': T2_0, 'I_0': I_0,
            'lam': lam, 'delta_I': delta_I,
            'beta': 4.e-5, 'p': 6.5e3, 'c': 20, 'w': 1.4e-2,
            'sigma': 1.8, 'V_0': 1.e-5, 'mu': 6,
        },
    ]

    T_list, t_list = inhost_viral_model_solver(conf, gConf, hosts)
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    ax.plot(t_list, np.log10(T_list[:, 1, 1, 1, 3]), 'r')
    ax.plot(t_list, np.log10(T_list[:, 3, 1, 1, 3]), 'b')
    ax.plot([-5, 35], [2, 2], ':')
    ax.set_xlim(-1, 31)
    ax.set_ylim(0, 10)
    ax.set_xlabel('time (days)')
    ax.set_ylabel('log10 copies/ml')
    plt.show()
