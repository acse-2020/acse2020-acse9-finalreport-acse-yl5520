#!/usr/bin/env python3
# encoding: utf-8
# version: 3.8.2

# standard modules
import argparse
import json
import logging
import sys
from typing import Tuple

# third-party modules

# local modules


def parseInput() -> Tuple[argparse.Namespace, argparse.ArgumentParser]:
    """ parse input and return arguments and parser

    :return 1: Namespace object from argparse
    :return 2: ArgumentParser object from argparse
    """
    parser = argparse.ArgumentParser(
        description='inhost viral dynamic solver driver'
    )
    parser.add_argument(
        'mConf',
        type=pathlib.Path,
        help='path/to/model-conf/json/file'
    )
    parser.add_argument(
        'hosts',
        type=pathlib.Path,
        help='path/to/host-conf/json/file'
    )
    parser.add_argument(
        'gConf',
        type=pathlib.Path,
        help='path/to/grid-conf/json/file'
    )
    args = parser.parse_args()

    return args, parser


def loadConst() -> None:
    """ update global constant from file
    """
    global T1_0, T2_0, I_0, lam, delta_I

    try:
        with open('conf/constant.json', 'r') as fin:
            const = json.load(fin)
        T1_0 = const['T1_0']
        T2_0 = const['T2_0']
        I_0 = const['I_0']
        lam = const['lam']
        delta_I = const['delta_I']
    except KeyError as err:
        logging.error(err)
        sys.exit(1)
    except FileNotFoundError as err:
        logging.error(err)
        sys.exit(1)


def loadModelConf(conf: pathlib.Path) -> None:
    """ update model configuration from file
    """
    global mConf

    try:
        with open(conf, 'r') as fin:
            conf = json.load(fin)
        mConf = {
            'dt': conf['d_sec'] / 86400,  # convert time step unit to day
            'ntime': conf['ntime'],  # output every n step
            'nstep': math.ceil(  # steps
                conf['t_day'] * 86400 / conf['d_sec'] / conf['ntime']
            ),
            'ndim_vel': conf['ndim_vel'],
            'nits': conf['nits'],
            'nits_solv_ng': conf['nits_solv_ng'],
            'error_solv': conf['error_solv'],
            'error_solv_ng': conf['error_solv_ng'],
            'relax': conf['relax']
        }
    except KeyError as err:
        logging.error(err)
        sys.exit(1)
    except FileNotFoundError as err:
        logging.error(err)
        sys.exit(1)


def loadHosts(conf: pathlib.Path) -> Hosts:
    """ update hosts instance from file
    """
    try:
        with open(conf, 'r') as fin:
            hosts = Hosts(json.load(fin))
    except FileNotFoundError as err:
        logging.error(err)
        sys.exit(1)


class Model:
    constPath = 'conf/constant.json'

    def __init__(
        self,
        modelConf: pathlib.Path,
        hostsConf: pathlib.Path,
        gridConf: pathlib.Path
    ):
        Model.loadHosts(hostsConf)
        Model.loadGrid(gridConf)


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
            h['delta_I'] = cls.delta_I
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
    URTs.append(T_list[:, :, :, :, 3])
    LRTs.append(T_list[:, :, :, :, 7])

    saveDict = {
        'URTs': URTs,
        'LRTs': LRTs,
        't_list': t_list
    }
    with open('tmp.pkl', 'wb') as f:
        pickle.dump(saveDict, f, pickle.HIGHEST_PROTOCOL)
