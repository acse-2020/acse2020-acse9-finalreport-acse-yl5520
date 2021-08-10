#!/usr/bin/env python3
# encoding: utf-8
# version: 3.8.2

# standard modules
from typing import Dict, List, Tuple

# third-party modeuls
import numpy as np

# local modules
import diffusion3d_virus_town_large as dv
from .Hosts import Hosts


def load_init_terms(gConf: Dict, hosts: Hosts) -> Tuple[List, List]:
    s = np.zeros((gConf['nx'], gConf['ny'], gConf['nz'], gConf['ng']))
    T = np.zeros((gConf['nx'], gConf['ny'], gConf['nz'], gConf['ng']))

    for h in hosts.hosts:
        s[h['x'], h['y'], h['z'], [1, 5]] = h['lam']
        T[h['x'], h['y'], h['z'], :] = (
            h['T1_0_U'], h['T2_0_U'], h['I_0_U'], h['V_0_U'],
            h['T1_0_L'], h['T2_0_L'], h['I_0_L'], h['V_0_L'],
            gConf['env_0']
        )

    return s, T


def update_A(
    xi: int, gConf: List, hosts: List,
    A_diag: List, A_off_diag: List, T: List, t: float
) -> None:
    def delta(delta_I, sigma, mu, t):
        return delta_I * np.e ** (sigma * (t - mu)) if (t >= mu) else delta_I

    for h in hosts.hosts:
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
            -xi * gConf['gamma_U'], -xi * gConf['gamma_L']
        )


def extract(T: List, hosts: Hosts) -> List:
    values = []
    for h in hosts.hosts:
        values.append(T[h['x'], h['y'], h['z'], :])
    return np.array(values)


def solver(
    mConf: Dict, gConf: Dict, hosts: Hosts, xis: List
) -> Tuple[List, List]:
    """ An inhost viral model solver using `sim_time_steping_diffusion_calc`
    to solve one or many within-host viral dynamic models simultaneously

    :param mConf: model configuration
    :param gConf: grid configuration
    :param hosts: hosts configuration
    :param xis: `xi` status at each timestep

    :return 1: T1, T2, I, V values of each host at each out-step
    :return 2: time of each out-step
    """

    """ === initialise matrices === """
    v_n = np.ones(gConf['ng'])

    # advection-diffusion in the system
    kdiff = np.zeros((gConf['nx'], gConf['ny'], gConf['nz'], gConf['ng']))
    kdiff[:, :, :, -1] = 1.e8

    # advection velocity
    u = np.zeros((
        mConf['ndim_vel'], gConf['nx'] + 1, gConf['ny'] + 1,
        gConf['nz'] + 1, gConf['ng']
    ))

    # diag and off-diag matrix in the matrix multiplication
    A_diag = np.zeros((gConf['nx'], gConf['ny'], gConf['nz'], gConf['ng']))
    A_off_diag = np.zeros((
        gConf['nx'], gConf['ny'], gConf['nz'], gConf['ng'], gConf['ng2']
    ))

    # source term and initial values
    s, T = load_init_terms(gConf, hosts)

    """ === start iteration === """
    # pre processing
    time_list = [0]
    T_list = [extract(T, hosts)]

    # iteration body
    for i, xi in enumerate(xis):
        print(f"{i + 1}/{mConf['nstep']}", end='\r', flush=True)
        t = time_list[-1] + mConf['dt'] * mConf['ntime']
        update_A(xi, gConf, hosts, A_diag, A_off_diag, T, t)

        T = dv.sim_time_steping_diffusion_calc(
            T, 1, mConf['nits'], mConf['nits_solv_ng'],
            mConf['relax'], mConf['error_solv'], mConf['error_solv_ng'],
            gConf['dx'], gConf['dy'], gConf['dz'], mConf['dt'], v_n, A_diag,
            A_off_diag, kdiff, s, u, gConf['i_upwind'], gConf['i_harmonic'],
            mConf['ndim_vel'], gConf['nx'], gConf['ny'], gConf['nz'],
            gConf['ng'], gConf['ng2']
        )

        if (i + 1) % mConf['ntime'] == 0:
            time_list.append(t)
            T_list.append(extract(T, hosts))

    # post processing
    print()
    if (i + 1) % mConf['ntime']:
        time_list.append(t)
        T_list.append(extract(T, hosts))

    return np.array(T_list), np.array(time_list)


if __name__ == '__main__':
    pass
