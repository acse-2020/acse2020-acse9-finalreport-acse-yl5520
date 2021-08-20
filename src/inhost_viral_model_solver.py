#!/usr/bin/env python3
# encoding: utf-8
# version: 3.8.2

# standard modules
from random import randint, uniform
from typing import Dict, List, Tuple

# third-party modeuls
import numpy as np
from tqdm import tqdm

# local modules
import diffusion3d_virus_town_large as dv
from .Hosts import Hosts


def load_init_values(gConf: Dict, hosts: Hosts) -> List:
    T = np.zeros((gConf['nx'], gConf['ny'], gConf['nz'], gConf['ng']))

    for h in hosts.host:
        T[h['x'], h['y'], h['z'], :] = (
            h['T1-U-0'], h['T2-U-0'], h['I-U-0'], h['V-U-0'],
            h['T1-L-0'], h['T2-L-0'], h['I-L-0'], h['V-L-0'],
            gConf['env_0']
        )

    return T


def update_matrix(
    xi: int, gConf: List, hosts: List,
    A_diag: List, A_off_diag: List, s: List,
    T: List, t: float, dt: float
) -> None:
    def delta(delta_I, sigma, mu, t):
        return delta_I * np.e ** (sigma * (t - mu)) if (t >= mu) else delta_I

    for h in hosts.host:
        # lymphocyte recruitment only during infection
        # prolong emergence day before infection
        if T[h['x'], h['y'], h['z'], 3] > 0:  # URT
            s[h['x'], h['y'], h['z'], 1] = h['lambda-U']
            h['mu-U'] += dt
        else:
            s[h['x'], h['y'], h['z'], 1] = 0
        if T[h['x'], h['y'], h['z'], 7] > 0:  # LRT
            s[h['x'], h['y'], h['z'], 5] = h['lambda-L']
            h['mu-L'] += dt
        else:
            s[h['x'], h['y'], h['z'], 5] = 0
            
        A_diag[h['x'], h['y'], h['z'], :] = (
            # T1_U, T2_U
            h['beta-U'] * T[h['x'], h['y'], h['z'], 3],
            h['beta-U'] * T[h['x'], h['y'], h['z'], 3],
            # I_U
            delta(h['delta-I-U'], h['sigma-U'], h['mu-U'], t)
            + h['w-U'] * T[h['x'], h['y'], h['z'], 1],
            # V_U
            h['c-U'],
            # T1_L, T2_L
            h['beta-L'] * T[h['x'], h['y'], h['z'], 7],
            h['beta-L'] * T[h['x'], h['y'], h['z'], 7],
            # I_L
            delta(h['delta-I-L'], h['sigma-L'], h['mu-L'], t)
            + h['w-L'] * T[h['x'], h['y'], h['z'], 5],
            # V_L
            h['c-L'],
            # V_env
            gConf['halflife'] + gConf['ventilation']
        )

        A_off_diag[h['x'], h['y'], h['z'],
                   [2, 2, 3, 3, 6, 6, 7, 7, 8],
                   [0, 1, 2, 8, 4, 5, 3, 6, 3]
                   ] = (
            h['beta-U'] * T[h['x'], h['y'], h['z'], 3],
            h['beta-U'] * T[h['x'], h['y'], h['z'], 3],
            h['p-U'], xi * h['gamma-inhale'],
            h['beta-L'] * T[h['x'], h['y'], h['z'], 7],
            h['beta-L'] * T[h['x'], h['y'], h['z'], 7],
            h['gamma-conduct'], h['p-L'],
            xi * h['gamma-shed']
        )


def extract(T: List, hosts: Hosts) -> List:
    values = []
    for h in hosts.host:
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
    kdiff = -np.ones((gConf['nx'], gConf['ny'], gConf['nz'], gConf['ng']))
    kdiff[1:gConf['nx'] - 1, 1:gConf['ny'] - 1, 1:gConf['nz'] - 1, -1] *= -1

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
    s = np.zeros((gConf['nx'], gConf['ny'], gConf['nz'], gConf['ng']))
    T = load_init_values(gConf, hosts)

    """ === start iteration === """
    # pre processing
    t = 0
    time_list = [t]
    T_list = [extract(T, hosts)]

    # iteration body
    for i, xi in enumerate(tqdm(xis, dynamic_ncols=True, disable=False)):
        update_matrix(xi, gConf, hosts, A_diag, A_off_diag, s, T, t, mConf['dt'])

        T = dv.sim_time_stepping_diffusion_calc(
            T, 1, mConf['nits'], mConf['nits_solv_ng'],
            mConf['relax'], mConf['error_solv'], mConf['error_solv_ng'],
            gConf['dx'], gConf['dy'], gConf['dz'], mConf['dt'], v_n, A_diag,
            A_off_diag, kdiff, s, u, gConf['i_upwind'], gConf['i_harmonic'],
            mConf['ndim_vel'], gConf['nx'], gConf['ny'], gConf['nz'],
            gConf['ng'], gConf['ng2']
        )

        t += mConf['dt']
        if (i + 1) % mConf['ntime'] == 0:
            time_list.append(t)
            T_list.append(extract(T, hosts))

    # post processing
    if (i + 1) % mConf['ntime']:
        time_list.append(t)
        T_list.append(extract(T, hosts))

    return np.array(T_list), np.array(time_list)


if __name__ == '__main__':
    pass
