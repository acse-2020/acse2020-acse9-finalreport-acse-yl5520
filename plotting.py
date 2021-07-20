#!/usr/bin/env python3
# python version 3.6.14

# standard modules
import pickle

# third-party modules
import matplotlib.pyplot as plt


if __name__ == '__main__':
    with open('tmp.pkl', 'rb') as f:
        loadDict = pickle.load(f)
    URTs = loadDict['URTs']
    LRTs = loadDict['LRTs']
    t_list = loadDict['t_list']
    a_inhale_list = loadDict['a_inhale_list']

    fig, axs = plt.subplots(2, 2, figsize=(24, 18))
    for ax, urt, lrt, a_conduct in zip(axs.flatten(), URTs, LRTs, a_inhale_list):
        ax.semilogy(t_list, urt, label='URT')
        ax.semilogy(t_list, lrt, label='LRT')
        ax.plot([-5, 35], [1e2, 1e2], ':')
        ax.set_xlim(-1, 31)
        ax.set_ylim(1, 1e14)
        ax.set_xlabel('time (day)')
        ax.set_ylabel('viral load')
        ax.set_title(f'conduct={a_conduct}')
        ax.legend(loc='best')
    plt.tight_layout()
    plt.savefig('tmp.png')
