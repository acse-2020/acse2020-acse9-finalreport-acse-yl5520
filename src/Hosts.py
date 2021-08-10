#!/usr/bin/env python3
# encoding: utf-8
# version: 3.8.2

# standard modules
import math
from typing import Dict, Tuple


class Hosts:
    host_keys = [
        'T1_0_U', 'T2_0_U', 'I_0_U', 'V_0_U',
        'T1_0_L', 'T2_0_L', 'I_0_L', 'V_0_L',
        'beta_U', 'p_U', 'c_U', 'w_U', 'sigma_U',
        'beta_L', 'p_L', 'c_L', 'w_L', 'sigma_L',
        'mu', 'a_conduct', 'a_inhale', 'lam', 'delta_I'
    ]

    def __init__(self, conf: Dict, const: Dict):
        self.total = conf['hints']['amount']
        self.sick = conf['hints']['sick']
        self.const = const
        self.hosts = []
        self.patients = []
        self.getHosts(conf)

    def getHosts(self, conf: Dict) -> None:
        """ subset and augment the input host list to fit `total` and `sick`

        :param conf: input hosts configuration
        """
        sick = 0
        healthy = 0

        for h in conf['hosts']:  # filter input hosts
            if (  # append sick
                (h['V_0_U'] > 0 or h['V_0_L'] > 0)
                and sick < self.sick
            ):
                sick += 1
                for k in self.host_keys:
                    if k not in h:
                        self.generateKey(k, h)
                self.hosts.append(h)
                self.patients.append(1)
            elif (  # append healthy
                h['V_0_U'] == 0 and h['V_0_L'] == 0
                and healthy < self.total - self.sick
            ):
                healthy += 1
                for k in self.host_keys:
                    if k not in h:
                        self.generateKey(k, h)
                self.hosts.append(h)
                self.patients.append(0)

        while sick < self.sick:
            self.generateHost(sick=True)
            self.patients.append(1)
        while len(self.hosts) < self.total:
            self.generateHost()
            self.patients.append(0)

        self.placeHost()

    def generateKey(self, k: str, h: Dict) -> None:
        if 'T1_0' in k:
            h[k] = self.const['T1_0']
        elif 'T2_0' in k:
            h[k] = self.const['T2_0']
        elif 'I_0' in k:
            h[k] = self.const['I_0']
        elif k in ('lam', 'delta_I'):
            h[k] = self.const[k]
        elif k == 'a_conduct':
            h[k] = 0.00625
        elif k == 'a_inhale':
            h[k] = 0.00625
        else:
            raise NotImplementedError(k)

    def generateHost(self, sick: bool = False) -> None:
        h = {
            'T1_0_U': self.const['T1_0'], 'T1_0_L': self.const['T1_0'],
            'T2_0_U': self.const['T2_0'], 'T2_0_L': self.const['T2_0'],
            'I_0_U': self.const['I_0'], 'I_0_L': self.const['I_0'],
            'lam': self.const['lam'], 'delta_I': self.const['delta_I'],
            'V_0_U': 0, 'V_0_L': 0,
            'beta_U': 0, 'p_U': 0, 'c_U': 0, 'w_U': 0, 'sigma_U': 0,
            'beta_L': 0, 'p_L': 0, 'c_L': 0, 'w_L': 0, 'sigma_L': 0,
            'mu': 0, 'a_conduct': 0, 'a_inhale': 0
        }

        if sick:
            h['V_0_U']: 1.e-4
            h['V_0_L']: 1.e-3
        self.hosts.append(h)

    def placeHost(self) -> None:
        self.size = math.ceil(math.sqrt(self.total))
        xs = [i // self.size * 2 + 1 for i in range(self.total)]
        ys = [i % self.size * 2 + 1 for i in range(self.total)]
        for h, x, y in zip(self.hosts, xs, ys):
            h['x'] = x
            h['y'] = y
            h['z'] = 1

    def getIdFromPosition(self, x: int, y: int) -> int:
        raise NotImplementedError

    def getPositionFromId(self, Id: int) -> Tuple[int, int]:
        raise NotImplementedError


if __name__ == '__main__':
    pass
