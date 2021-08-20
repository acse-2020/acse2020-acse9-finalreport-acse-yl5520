#!/usr/bin/env python3
# encoding: utf-8
# version: 3.8.2

# standard modules
from typing import Dict

# third-party modules
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm


def plot_solution(attr: Dict, mode: str, figName: str):
    ls = [
        'T1-URT', 'T2-URT', 'I-URT', 'V-URT',
        'T1-LRT', 'T2-LRT', 'I-LRT', 'V-LRT',
        'V-ENV'
    ]
    fig, axs = plt.subplots(
        attr['size'], attr['size'],
        figsize=(8 * attr['size'], 5 * attr['size'])
    )
    if not isinstance(axs, np.ndarray):
        axs = np.array([axs])

    for i, ax in tqdm(list(enumerate(axs.flatten()[:attr['data'].shape[1]]))):
        for j in range(len(mode)):
            if mode[j] == '1':
                ax.plot(attr['t'], attr['data'][:, i, j], label=f'{ls[j]}')
        ax.set_yscale('symlog')
        ax.set_ylim(bottom=1e-9)
        ax.set_xlabel('time (day)')
        ax.set_xticks([i for i in range(int(attr['t'][-1]) + 1)])
        ax.set_title(attr['title'][i])
        ax.legend(loc='best')

    plt.tight_layout()
    if figName:
        print(f'saving figure {figName}')
        plt.savefig(figName)
    else:
        plt.show()


if __name__ == '__main__':
    pass
