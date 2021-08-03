#!/usr/bin/env python3
# encoding: utf-8
# version: 3.8.2

# standard modules
import argparse
import json
import logging
import math
import pathlib
import pickle
import sys
from typing import Dict, List, Tuple

# third-party modules
import numpy as np

# local modules
import diffusion3d_virus_town_large as dv




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
        'hostsConf',
        type=pathlib.Path,
        help='path/to/host-conf/json/file'
    )
    parser.add_argument(
        'gridConf',
        type=pathlib.Path,
        help='path/to/grid-conf/json/file'
    )
    args = parser.parse_args()

    return args, parser


class Model:
    constPath = 'conf/constant.json'

    def __init__(
        self,
        modelConf: pathlib.Path,
        hostsConf: pathlib.Path,
        gridConf: pathlib.Path
    ):
        Model.loadConst()
        Model.loadModel(modelConf)
        Model.loadHosts(hostsConf)
        Model.loadGrid(gridConf)

    def load_init_terms(self) -> Tuple[np.ndarray, np.ndarray]:
        s = np.zeros((self.gConf['nx'], self.gConf['ny'],
                      self.gConf['nz'], self.gConf['ng']))
        T = np.zeros((self.gConf['nx'], self.gConf['ny'],
                      self.gConf['nz'], self.gConf['ng']))

        for h in self.hosts:
            s[h['x'], h['y'], h['z'], [1, 5]] = h['lam']
            T[h['x'], h['y'], h['z'], :] = (
                h['T1_0_U'], h['T2_0_U'], h['I_0_U'], h['V_0_U'],
                h['T1_0_L'], h['T2_0_L'], h['I_0_L'], h['V_0_L'],
                self.gConf['env_0']
            )

        return s, T

    def update_A(self, A_diag: np.ndarray, A_off_diag: np.ndarray,
                 T: np.ndarray, t: float) -> None:
        def delta(delta_I, sigma, mu, t):
            return delta_I * np.e ** (sigma * (t - mu)) if (t >= mu) else delta_I

        for h in self.hosts:
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
                self.gConf['halflife'] + self.gConf['ventilation']
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
                -self.gConf['gamma_U'], -self.gConf['gamma_L']
            )

    def inhost_viral_model_solver(self):
        """ An inhost viral model solver using `sim_time_steping_diffusion_calc`
        to solve one or many within-host viral dynamic models simultaneously
        """
        # initialise matrices
        #
        v_n = np.ones(self.gConf['ng'])
        # advection-diffusion in the system
        kdiff = np.zeros((self.gConf['nx'], self.gConf['ny'],
                          self.gConf['nz'], self.gConf['ng']))
        kdiff[:, :, :, -1] = 1.e8
        # advection velocity
        u = np.zeros((
            self.conf['ndim_vel'], self.gConf['nx'] + 1, self.gConf['ny'] + 1,
            self.gConf['nz'] + 1, self.gConf['ng']
        ))
        #
        A_diag = np.zeros((self.gConf['nx'], self.gConf['ny'],
                           self.gConf['nz'], self.gConf['ng']))
        A_off_diag = np.zeros((
            self.gConf['nx'], self.gConf['ny'], self.gConf['nz'],
            self.gConf['ng'], self.gConf['ng2']
        ))
        # source term and initial values
        s, T = self.load_init_terms()

        time_list = [0]
        T_list = [T]
        for i in range(self.conf['nstep']):
            print(f"{i + 1}/{self.conf['nstep']}", end='          \r', flush=True)
            t = time_list[-1] + self.conf['dt'] * self.conf['ntime']
            self.update_A(A_diag, A_off_diag, T, t)

            T = dv.sim_time_steping_diffusion_calc(
                T, self.conf['ntime'], self.conf['nits'], self.conf['nits_solv_ng'],
                self.conf['relax'], self.gConf['error_solv'],
                self.gConf['error_solv_ng'],
                self.gConf['dx'], self.gConf['dy'], self.gConf['dz'], self.conf['dt'],
                v_n, A_diag, A_off_diag, kdiff, s, u,
                self.gConf['i_upwind'], self.gConf['i_harmonic'],
                self.conf['ndim_vel'], self.gConf['nx'], self.gConf['ny'], self.gConf['nz'],
                self.gConf['ng'], self.gConf['ng2']
            )

            time_list.append(t)
            T_list.append(T)
        print()

        return np.array(T_list), time_list

    @classmethod
    def loadConst(cls) -> None:
        try:
            with open(cls.constPath, 'r') as f:
                const = json.load(f)
            cls.T1_0 = const['T1_0']
            cls.T2_0 = const['T2_0']
            cls.I_0 = const['I_0']
            cls.lam = const['lam']
            cls.delta_I = const['delta_I']
        except KeyError as err:
            logging.error(err)
            sys.exit(1)
        except FileNotFoundError as err:
            logging.error(err)
            sys.exit(1)

    @classmethod
    def loadModel(cls, conf: pathlib.Path) -> None:
        try:
            with open(conf, 'r') as f:
                conf = json.load(f)
            cls.d_sec = conf['d_sec']  # amount of second for each time step
            cls.t_sec = conf['t_day'] * 86400  # amount of second
            cls.ntime = conf['ntime']  # output every n step
            cls.conf = {
                'dt': cls.d_sec / 86400,  # convert time step unit to day
                'ntime': conf['ntime'],
                'nstep': math.ceil(cls.t_sec / cls.d_sec / cls.ntime),  # steps
                'ndim_vel': conf['ndim_vel'],
                'nits': conf['nits'],
                'nits_solv_ng': conf['nits_solv_ng'],
                'relax': conf['relax']
            }
        except KeyError as err:
            logging.error(err)
            sys.exit(1)
        except FileNotFoundError as err:
            logging.error(err)
            sys.exit(1)

    @classmethod
    def loadHosts(cls, conf: pathlib.Path) -> None:
        try:
            with open(conf, 'r') as f:
                conf = json.load(f)
            cls.hostAmount = conf['hints']['amount']
            cls.generateHost(conf['hints']['sick'], conf['hosts'])
        except KeyError as err:
            logging.error(err)
            sys.exit(1)
        except ValueError as err:
            logging.error(err)
            sys.exit(1)
        except FileNotFoundError as err:
            logging.error(err)
            sys.exit(1)

    @classmethod
    def loadGrid(cls, conf: pathlib.Path) -> None:
        try:
            with open(conf, 'r') as f:
                conf = json.load(f)
            cls.gConf = {
                'nx': cls.size * 2 + 1, 'ny': cls.size * 2 + 1, 'nz': 3,
                'ng': 9, 'ng2': 9,
                "error_solv": conf['error_solv'],
                "error_solv_ng": conf['error_solv_ng'],
                "i_upwind": 1, "i_harmonic": 0,
                "dx": 0.1, "dy": 0.1, "dz": 0.1,
                "env_0": conf['env_0'], "halflife": conf['halflife'],
                "ventilation": conf['ventilation'],
                "gamma_U": conf['gamma_U'], "gamma_L": conf['gamma_L']
            }
        except KeyError as err:
            logging.error(err)
            sys.exit(1)
        except FileNotFoundError as err:
            logging.error(err)
            sys.exit(1)

    @classmethod
    def generateHost(cls, sick: int, hosts: List[Dict]) -> None:
        for h in hosts:
            h['T1_0_U'] = cls.T1_0
            h['T1_0_L'] = cls.T1_0
            h['T2_0_U'] = cls.T2_0
            h['T2_0_L'] = cls.T2_0
            h['I_0_U'] = cls.I_0
            h['I_0_L'] = cls.I_0
            h['lam'] = cls.lam
            h['delta_I']= cls.delta_I
            if h['V_0_U'] > 0 or h['V_0_L'] > 0:
                sick -= 1
        if sick < 0:
            raise ValueError('number of sick hosts exceed hints.sick')

        for _ in range(sick):
            hosts.append(cls.newHost(sick=True))
        while len(hosts) < cls.hostAmount:
            hosts.append(cls.newHost())
        cls.placeHost(hosts)

    @classmethod
    def newHost(cls, sick: bool = False) -> Dict:
        h = {
            'T1_0_U': cls.T1_0, 'T1_0_L': cls.T1_0,
            'T2_0_U': cls.T2_0, 'T2_0_L': cls.T2_0,
            'I_0_U': cls.I_0, 'I_0_L': cls.I_0,
            'lam': cls.lam, 'delta_I': cls.delta_I,
            'V_0_U': 0, 'V_0_L': 0,
            'beta_U': 0, 'p_U': 0, 'c_U': 0,
            'w_U': 0, 'sigma_U': 0,
            'beta_L': 0, 'p_L': 0, 'c_L': 0,
            'w_L': 0, 'sigma_L': 0,
            'mu': 0,
            'a_conduct': 0, 'a_inhale': 0
        }

        if sick:
            h['V_0_U']: 1.e-4
            h['V_0_L']: 1.e-3
        return h

    @classmethod
    def placeHost(cls, hosts: List[Dict]) -> None:
        cls.size = math.ceil(math.sqrt(cls.hostAmount))
        xs = [i // cls.size * 2 + 1 for i in range(len(hosts))]
        ys = [i % cls.size * 2 + 1 for i in range(len(hosts))]
        for h, x, y in zip(hosts, xs, ys):
            h['x'] = x
            h['y'] = y
            h['z'] = 1
        cls.hosts = hosts


if __name__ == '__main__':
    args, parser = parseInput()

    model = Model(args.modelConf, args.hostsConf, args.gridConf)

    URTs = []
    LRTs = []
    T_list, t_list = model.inhost_viral_model_solver()
    URTs.append(T_list[:, 1, 1, 1, 3])
    LRTs.append(T_list[:, 1, 1, 1, 7])

    saveDict = {
        'URTs': URTs,
        'LRTs': LRTs,
        't_list': t_list
    }
    with open('tmp.pkl', 'wb') as f:
        pickle.dump(saveDict, f, pickle.HIGHEST_PROTOCOL)
