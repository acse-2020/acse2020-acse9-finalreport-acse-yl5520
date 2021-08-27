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
from datetime import datetime, timedelta
from pathlib import Path
from typing import Dict, Tuple

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
        '--constant',
        default='conf/constant-default.json',
        help='path/to/model-conf/json-file'
    )
    parser.add_argument(
        '--mConf',
        default='conf/model-conf-default.json',
        help='path/to/model-conf/json-file'
    )
    parser.add_argument(
        '--hosts',
        default='conf/hosts-conf-default.json',
        help='path/to/hosts-conf/json-file'
    )
    parser.add_argument(
        '--dist',
        default='conf/distribution.json',
        help='path/to/parameter-distribution/json-file'
    )
    parser.add_argument(
        '--gConf',
        default='conf/grid-conf-default.json',
        help='path/to/grid-conf/json-file'
    )
    parser.add_argument(
        '--schedule',
        default='conf/schedule-default.json',
        help='path/to/schedule/json-file'
    )
    parser.add_argument(
        '--output',
        type=Path,
        help='path/to/store/pkl-file'
    )
    parser.add_argument(
        '--plot',
        type=Path,
        help='path/to/load/pkl-file'
    )
    parser.add_argument(
        '--mode', default='000100010',
        help='9-bits binary number indicate which to plot, ' +
             'T1, T2, I, V for URT and LRT, then virus in Environment.'
    )
    parser.add_argument(
        '--figname',
        type=Path,
        help='path/to/store/plotting'
    )
    args = parser.parse_args()

    return args


def loadConst() -> None:
    """ update global constant from file
    """
    global const

    try:
        with open(args.constant, 'r') as fin:
            const = json.load(fin)
        const = {
            'T1-U-0': const['default-URT-initial-primary-target-cells-concentration'],
            'T1-L-0': const['default-LRT-initial-primary-target-cells-concentration'],
            'T2-U-0': const['default-URT-initial-secondary-target-cells-concentration'],
            'T2-L-0': const['default-LRT-initial-secondary-target-cells-concentration'],
            'I-U-0': const['default-URT-initial-infected-cells-concentration'],
            'I-L-0': const['default-LRT-initial-infected-cells-concentration'],
            'lambda-U': const['default-URT-secondary-target-cells-recruitment-rate'],
            'lambda-L': const['default-LRT-secondary-target-cells-recruitment-rate'],
            'delta-I-U': const['default-URT-infected-cells-base-killing-rate'],
            'delta-I-L': const['default-LRT-infected-cells-base-killing-rate']
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
            'dt': conf['timestep-size-sec'] / 86400,
            'ntime': conf['output-step'],
            'nstep': math.ceil(
                conf['simulation-duration-day'] * 86400 / conf['timestep-size-sec']
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


def loadHosts(host_conf: Path, dist_conf: Path) -> None:
    """ update hosts instance from file

    :param conf: path/to/hosts-conf/json-file
    """
    global hosts

    try:
        with open(host_conf, 'r') as fin_host, open(dist_conf, 'r') as fin_dist:
            hosts = Hosts(const, json.load(fin_host), json.load(fin_dist))
    except FileNotFoundError as err:
        logging.error(err)
        sys.exit(1)


def loadGrid(conf: Path) -> None:
    """ update grid configuration from file

    :param conf: path/to/grid-conf/json-file
    """
    global gConf

    try:
        with open(conf, 'r') as fin:
            conf = json.load(fin)
        length = hosts.length
        gConf = {
            'nx': length * 2 + 1, 'ny': length * 2 + 1, 'nz': 3,
            'ng': 9, 'ng2': 9,
            'i_upwind': 1, 'i_harmonic': 0,
            'dx': 0.1, 'dy': 0.1, 'dz': 0.1,
            'env_0': conf['initial-environment-virus-concentration'],
            'halflife': 1 - (1 / 2 ** (1 / (conf['in-vitro-virus-min-halflife'] / 1440))),
            'ventilation': conf['environment-day-ventilation-rate']
        }
#         gConf['halflife'] = 0
#         gConf['ventilation'] = 0
    except KeyError as err:
        logging.error(err)
        sys.exit(1)
    except FileNotFoundError as err:
        logging.error(err)
        sys.exit(1)


def getStartDatetime(start: str, days: Tuple) -> datetime:
    day, time = start.split(' ')
    day = days.index(day.lower())
    hour, minute = int(time[:2]), int(time[2:])
    return datetime(1, 1, day, hour, minute, 0, 0)


def getPeriod(period: Dict) -> Dict:
    rtn = {}
    for k, v in period.items():
        rtn[k.lower()] = v.split('-')
    return rtn


def loadSchedule(conf: Path) -> None:
    """ create function `xi` from schedule configuration

    :param conf: path/to/schedule/json-file

    :return 1: `xi` function
    """
    global mConf, xis
    days = ('', 'mon', 'tue', 'wed', 'thu', 'fri', 'sat', 'sun')

    try:
        with open(conf, 'r') as fin:
            conf = json.load(fin)
        xis = [conf['type']] * mConf['nstep']  # default type

        if conf['details'].items():  # exclude non period time
            datetime = getStartDatetime(conf['start'], days)
            period = getPeriod(conf['details'])

            delta = timedelta(days=mConf['dt'])
            for i in range(mConf['nstep']):
                datetime += delta
                time = f'{datetime.hour:02}{datetime.minute:02}'

                # not within period
                if (
                    days[datetime.weekday()] not in period
                    or (
                        period[days[datetime.weekday()]][0] > time
                        or period[days[datetime.weekday()]][1] < time
                    )
                ):
                    xis[i] = 1 - conf['type']
    except KeyError as err:
        logging.error(err)
        sys.exit(1)
    except FileNotFoundError as err:
        logging.error(err)
        sys.exit(1)


if __name__ == '__main__':
    args = parseInput()

    # run simulation and store result
    if args.output:
        loadConst()
        loadModel(args.mConf)
        loadHosts(args.hosts, args.dist)
        loadGrid(args.gConf)
        loadSchedule(args.schedule)

        host_values, t_values = solver(mConf, gConf, hosts, xis)
        attr = {
            'size': hosts.length,
            'title': [h['attribute'] for h in hosts.host],
            'data': host_values,
            't': t_values
        }

        if args.output.suffix == '.pkl':
            args.output = Path(args.output.stem)
        with open(f'data/{args.output.stem}.pkl', 'wb') as fout:
            pickle.dump(attr, fout, pickle.HIGHEST_PROTOCOL)

    # plot simulation result
    elif args.plot:
        with open(args.plot, 'rb') as fin:
            attr = pickle.load(fin)

        if args.figname and args.figname[-4:] == '.png':
            args.figname = f'images/{args.figname[:-4]}'

        plot_solution(attr, args.mode, args.figname)
