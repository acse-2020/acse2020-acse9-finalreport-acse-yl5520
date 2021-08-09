#!/usr/bin/env python3
# encoding: utf-8
# version: 3.8.2

# standard modules
from typing import Dict


class Hosts:
    host_keys = [
        'beta_U', 'p_U', 'c_U', 'w_U', 'sigma_U',
        'beta_L', 'p_L', 'c_L', 'w_L', 'sigma_L',
        'mu', 'a_conduct', 'a_inhale'
    ]

    def __init__(self, conf: Dict):
        self.total = conf['hints']['amount']
        self.sick = conf['hints']['sick']
        self.hosts = []
        generateHosts()

    def generateHosts(self) -> None:
        """ subset and augment the input host list to fit `total` and `sick`
        """
        sick = 0

        # filter input hosts
        for h in conf['hosts']:
            if (
                (h['V_0_U'] > 0 or h['V_0_L'] > 0)
                and sick < self.sick
            ):
                sick += 1
                for k in self.host_keys:
                    if k not in h:
                        self.generateKey(k, h)
                self.hosts.append(h)

    def generateKey(self, k: str, h: Dict) -> None:
        raise NotImplementedError
