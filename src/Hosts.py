#!/usr/bin/env python3
# encoding: utf-8
# version: 3.8.2

# standard modules
import math
from random import gauss, uniform
from typing import Dict, Tuple


class Hosts:
    def __init__(self, const: Dict, conf: Dict, dist: Dict):
        self.const = const
        self.dist = dist
        self.host = []
        self.generateHosts(conf)

    def getTypeStringFromTypeId(self, typeId: int) -> str:
        attrs = list(f'{typeId:03b}')
        for i, name in zip(range(len(attrs)), ['mask', 'vaccine', 'infect']):
            attrs[i] = name if int(attrs[i]) else 'un' + name
        return '-'.join(attrs)

    def getTypeIdFromTypeString(self, typeStr: str) -> int:
        digits = typeStr.split('-')
        for i in range(len(digits)):
            digits[i] = '0' if digits[i][:2] == 'un' else '1'
        return int(''.join(digits), base=2)

    def getPosition(self) -> Tuple[int, int]:
        size = len(self.host)
        return (size // self.length * 2 + 1, size % self.length * 2 + 1)

    def getMusAndSigmas(self, genType: int) -> Dict:
        distribution = {
            'beta-U': self.dist['URT-infection-rate']
                               [self.getTypeStringFromTypeId(genType).split('-')[1]],
            'beta-L': self.dist['LRT-infection-rate']
                               [self.getTypeStringFromTypeId(genType).split('-')[1]],
            'p-U': self.dist['URT-virus-production-rate']['normal'],
            'p-L': self.dist['LRT-virus-production-rate']['normal'],
            'c-U': self.dist['URT-virus-clearance-rate']
                            [self.getTypeStringFromTypeId(genType).split('-')[1]],
            'c-L': self.dist['LRT-virus-clearance-rate']
                            [self.getTypeStringFromTypeId(genType).split('-')[1]],
            'w-U': self.dist['URT-lymphocyte-killing-rate']['normal'],
            'w-L': self.dist['LRT-lymphocyte-killing-rate']['normal'],
            'sigma-U': self.dist['URT-killing-rate-rate']['normal'],
            'sigma-L': self.dist['LRT-killing-rate-rate']['normal'],
            'mu': self.dist['adaptive-immune-response-emergence-day']
                            [self.getTypeStringFromTypeId(genType).split('-')[1]],
        }

        return distribution

    def getPositiveGaussRandom(self, mean: float, std: float) -> float:
        rtn = gauss(mean, std)
        while rtn <= 0:
            rtn = gauss(mean, std)
        return rtn

    def generateHost(self, genType: int, genAmount: int) -> None:
        x, y = self.getPosition()
        distribution = self.getMusAndSigmas(genType)
        for _ in range(genAmount):
            host = {
                'attribute': self.getTypeStringFromTypeId(genType),
                'x': x, 'y': y, 'z': 1,
                'gamma-conduct': .01,
                'lambda-U': self.const['lambda-U'], 'lambda-L': self.const['lambda-L'],
                'delta-I-U': self.const['delta-I-U'], 'delta-I-L': self.const['delta-I-L'],
                'T1-U-0': self.const['T1-U-0'], 'T2-U-0': self.const['T2-U-0'],
                'I-U-0': self.const['I-U-0'],
                'V-U-0': uniform(1e-6, 1e-3) if genType % 2 else 0,
                'beta-U': self.getPositiveGaussRandom(*distribution['beta-U']),
                'p-U': self.getPositiveGaussRandom(*distribution['p-U']),
                'c-U': self.getPositiveGaussRandom(*distribution['c-U']),
                'w-U': self.getPositiveGaussRandom(*distribution['w-U']),
                'sigma-U': self.getPositiveGaussRandom(*distribution['sigma-U']),
                'T1-L-0': self.const['T1-L-0'], 'T2-L-0': self.const['T2-L-0'],
                'I-L-0': self.const['I-L-0'], 'V-L-0': 0,
                'beta-L': self.getPositiveGaussRandom(*distribution['beta-L']),
                'p-L': self.getPositiveGaussRandom(*distribution['p-L']),
                'c-L': self.getPositiveGaussRandom(*distribution['c-L']),
                'w-L': self.getPositiveGaussRandom(*distribution['w-L']),
                'sigma-L': self.getPositiveGaussRandom(*distribution['sigma-L']),
                'gamma-inhale': uniform(.45, .55) if genType < 4 else 1e-20,
                'gamma-shed': uniform(.30, .40) * .01 if genType < 4 else 1e-20
            }
            host['mu-U'] = 0 if genType // 2 % 2 else self.getPositiveGaussRandom(*distribution['mu'])
            host['mu-L'] = host['mu-U']
            self.host.append(host)

    def generateHosts(self, conf: Dict) -> None:
        """ subset and augment the input host list to fit `total` and `sick`

        :param conf: input hosts configuration
        """
        self.total = sum(conf.values())
        self.length = math.ceil(math.sqrt(self.total))
        for k, v in conf.items():
            self.generateHost(self.getTypeIdFromTypeString(k), v)
        print(self.host)
        exit(1)


if __name__ == '__main__':
    pass
