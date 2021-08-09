#!/usr/bin/env python3
# encoding: utf-8
# version: 3.8.2

# standard modules
import pickle

# third-party modules
import matplotlib.pyplot as plt


if __name__ == '__main__':
    with open('tmp.pkl', 'rb') as f:
        loadDict = pickle.load(f)
    URTs = loadDict['URTs'][0]
    LRTs = loadDict['LRTs'][0]
    t_list = loadDict['t_list']

    fig, ax = plt.subplots(1, 1, figsize=(12, 9))
    for (lrt, inhale) in zip(
        [LRTs[:, 1, 1, 1], LRTs[:, 1, 3, 1], LRTs[:, 3, 1, 1], LRTs[:, 3, 3, 1]],
        ['0.0625', '1e-2', '1e-3', '1e-4']
    ):
        ax.semilogy(t_list, lrt, label=f'LRT {inhale}')
    ax.plot([-5, 35], [1e2, 1e2], ':')
    ax.set_xlim(-1, 31)
    ax.set_ylim(1, 1e14)
    ax.set_xlabel('time (day)')
    ax.set_ylabel('viral load')
    ax.set_title('LRT comparison')
    ax.legend(loc='best')
    plt.savefig('tmp.png')
    exit(0)

    fig, axs = plt.subplots(2, 2, figsize=(24, 18))
    for (ax, urt, lrt, inhale) in zip (
        axs.flatten(),
        [URTs[:, 1, 1, 1], URTs[:, 1, 3, 1], URTs[:, 3, 1, 1], URTs[:, 3, 3, 1]],
        [LRTs[:, 1, 1, 1], LRTs[:, 1, 3, 1], LRTs[:, 3, 1, 1], LRTs[:, 3, 3, 1]],
        ['0.0625', '1e-2', '1e-3', '1e-4']
    ):
        ax.semilogy(t_list, urt, label='URT')
        ax.semilogy(t_list, lrt, label='LRT')
        ax.plot([-5, 35], [1e2, 1e2], ':')
        ax.set_xlim(-1, 31)
        ax.set_ylim(1, 1e14)
        ax.set_xlabel('time (day)')
        ax.set_ylabel('viral load')
        ax.set_title(f'inhale {inhale}')
        ax.legend(loc='best')
    plt.tight_layout()
    plt.savefig('tmp.png')
