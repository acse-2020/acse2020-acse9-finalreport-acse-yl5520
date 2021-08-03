#!/usr/bin/env python3
# encoding: utf-8
# version: 3.8.2

# standard modules
import argparse
import json
import math
import pathlib
import pickle
import sys
from typing import Tuple

# third-party modules
import numpy as np

# local modules
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
    kdiff[:, :, :, -1] = 1.e8
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


def parseInput() -> Tuple[argparse.Namespace, argparse.ArgumentParser]:
    """ parse input and return arguments and parser
    """
    parser = argparse.ArgumentParser(
        description='inhost viral dynamic model'
    )
    parser.add_argument(
        'modelConf',
        type=pathlib.Path,
        help='path/to/model-conf/json/file'
    )
    parser.add_argument(
        'gridConf',
        type=pathlib.Path,
        help='path/to/grid-conf/json/file'
    )
    parser.add_argument(
        'hostConf',
        type=pathlib.Path,
        help='path/to/host-conf/json/file'
    )
    args = parser.parse_args()

    return args, parser


const_sec_per_day = 86400
if __name__ == '__main__':
    args, parser = parseInput()

    with open(args.modelConf, 'r') as f:
        args.modelConf = json.load(f)

    # global constants (Wang, S. et al. 2020)
    T1_0 = args.modelConf['T1_0']
    T2_0 = args.modelConf['T2_0']
    I_0 = args.modelConf['I_0']
    lam = args.modelConf['lam']
    delta_I = args.modelConf['delta_I']

    # hyper-param
    d_sec = args.modelConf['d_sec']  # amount of second for each time step
    t_sec = args.modelConf['t_day'] * const_sec_per_day  # amount of second
    ntime = args.modelConf['ntime']  # output every n step

    # prepare model configuration
    del args.modelConf['T1_0']
    del args.modelConf['T2_0']
    del args.modelConf['I_0']
    del args.modelConf['lam']
    del args.modelConf['delta_I']
    del args.modelConf['d_sec']
    del args.modelConf['t_day']
    args.modelConf['dt'] = d_sec / const_sec_per_day  # convert time step unit
    args.modelConf['nstep'] = math.ceil(t_sec / d_sec / ntime)  # total steps

    # prepare grid configuration
    with open(args.gridConf, 'r') as f:
        args.gridConf = json.load(f)

    # prepare hosts configuration
    with open(args.hostConf, 'r') as f:
        args.hostConf = json.load(f)
    for i in range(len(args.hostConf)):
        args.hostConf[i]['T1_0_U'] = T1_0
        args.hostConf[i]['T1_0_L'] = T1_0
        args.hostConf[i]['T2_0_U'] = T2_0
        args.hostConf[i]['T2_0_L'] = T2_0
        args.hostConf[i]['I_0_U'] = I_0
        args.hostConf[i]['I_0_L'] = I_0
        args.hostConf[i]['lam'] = lam
        args.hostConf[i]['delta_i'] = delta_I

    """
    a_inhale_list = [2 ** -i for i in range(3)] + [0]
    URTs = []
    LRTs = []
    for a_inhale in a_inhale_list:
        hosts[0]['a_inhale'] = a_inhale
        T_list, t_list = inhost_viral_model_solver(
            args.modelConf, args.gridConf, hosts
        )
        URTs.append(T_list[:, 1, 1, 1, 3])
        LRTs.append(T_list[:, 1, 1, 1, 7])

    saveDict = {
        'URTs': URTs,
        'LRTs': LRTs,
        't_list': t_list,
        'a_inhale_list': a_inhale_list
    }
    with open('tmp.pkl', 'wb') as f:
        pickle.dump(saveDict, f, pickle.HIGHEST_PROTOCOL)
    """
