#!/usr/bin/env python3
# encoding: utf-8
# version: 3.8.2

# standard modules
import argparse
import json
import logging
import math
import pickle
import sys
from pathlib import Path

# local modules
from src.Hosts import Hosts
from src.inhost_viral_model_solver import solver
from src.plotting import plot_solution


def parseInput() -> argparse.Namespace:
    """ parse input and return arguments and parser

    :return 1: Namespace object from argparse
    """
    parser = argparse.ArgumentParser(
        description='inhost viral dynamic solver driver'
    )
    parser.add_argument(
        'mConf', nargs='?',
        default='conf/model-conf.json',
        help='path/to/model-conf/json-file'
    )
    parser.add_argument(
        'hosts', nargs='?',
        default='conf/hosts-conf.json',
        help='path/to/hosts-conf/json-file'
    )
    parser.add_argument(
        'gConf', nargs='?',
        default='conf/grid-conf.json',
        help='path/to/grid-conf/json-file'
    )
    parser.add_argument(
        '--output',
        help='path/to/store/pkl-file'
    )
    parser.add_argument(
        '--plot', nargs='?',
        action='append',
        help='plot the input `pkl` or the newly calculate values'
    )
    parser.add_argument(
        '--mode', default='000100010',
        help='9-bits binary number indicate which to plot, ' +
             'T1, T2, I, V for URT and LRT, then virus in Environment.'
    )
    parser.add_argument(
        '--figName',
        help='path/to/store/plotting'
    )
    args = parser.parse_args()

    return args


def loadConst() -> None:
    """ update global constant from file
    """
    global const

    try:
        with open('conf/constant.json', 'r') as fin:
            const = json.load(fin)
        const = {
            'T1_0': const['T1_0'],
            'T2_0': const['T2_0'],
            'I_0': const['I_0'],
            'lam': const['lam'],
            'delta_I': const['delta_I']
        }
    except KeyError as err:
        logging.error(err)
        sys.exit(1)
    except FileNotFoundError as err:
        logging.error(err)
        sys.exit(1)


def loadModel(conf: Path) -> None:
    """ update model configuration from file

    :param conf: path/to/model-conf/json-file
    """
    global mConf

    try:
        with open(conf, 'r') as fin:
            conf = json.load(fin)
        mConf = {
            'dt': conf['d_sec'] / 86400,  # convert time step unit to day
            'ntime': conf['ntime'],  # output every n step
            'nstep': math.ceil(  # steps
                conf['t_day'] * 86400 / conf['d_sec']
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


def loadHosts(conf: Path) -> None:
    """ update hosts instance from file

    :param conf: path/to/hosts-conf/json-file
    """
    global hosts

    try:
        with open(conf, 'r') as fin:
            hosts = Hosts(json.load(fin), const)
    except FileNotFoundError as err:
        logging.error(err)
        sys.exit(1)


def loadGrid(conf: Path, hosts: Hosts) -> None:
    """ update grid configuration from file

    :param conf: path/to/grid-conf/json-file
    :param hosts: Hosts instance
    """
    global gConf

    try:
        with open(conf, 'r') as fin:
            conf = json.load(fin)
        gConf = {
            'nx': hosts.size * 2 + 1, 'ny': hosts.size * 2 + 1, 'nz': 3,
            'ng': 9, 'ng2': 9,
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


if __name__ == '__main__':
    args = parseInput()

    if not args.plot or not args.plot[0]:
        loadConst()
        loadModel(args.mConf)
        loadHosts(args.hosts)
        loadGrid(args.gConf, hosts)

        host_values, t_values = solver(mConf, gConf, hosts)
        attr = {
            'size': hosts.size,
            'patients': hosts.patients,
            'data': host_values,
            't': t_values
        }

        with open(f'{args.output}.pkl', 'wb') as fout:
            pickle.dump(attr, fout, pickle.HIGHEST_PROTOCOL)

    if args.plot:
        if args.plot[0]:
            with open(args.plot[0], 'rb') as fin:
                attr = pickle.load(fin)
        plot_solution(attr, args.mode, args.figName)
