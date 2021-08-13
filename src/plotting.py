#!/usr/bin/env python3
# encoding: utf-8
# version: 3.8.2

# standard modules
from typing import Dict

# third-party modules
import matplotlib.pyplot as plt
import numpy as np


def plot_solution(attr: Dict, mode: str, figName: str):
    ls = [
        'T1_URT', 'T2_URT', 'I_URT', 'V_URT',
        'T1_LRT', 'T2_LRT', 'I_LRT', 'V_LRT',
        'V_ENV'
    ]
    fig, axs = plt.subplots(
        attr['size'], attr['size'],
        figsize=(8 * attr['size'], 5 * attr['size'])
    )
    if not isinstance(axs, np.ndarray):
        axs = np.array([axs])

    for i, ax in enumerate(axs.flatten()):
        for j in range(len(mode)):
            if mode[j] == '1':
                ax.semilogy(attr['t'], attr['data'][:, i, j], label=f'{ls[j]}')
        ax.set_ylim(1, 1e15)
        ax.set_xlabel('time (day)')
        title = 'normal'
        if attr['patients'][i] == 1:
            title = 'patient'
        ax.set_title(title)
        ax.legend(loc='best')

    plt.tight_layout()
    if figName:
        plt.savefig(figName)
    plt.show()


if __name__ == '__main__':
    pass
