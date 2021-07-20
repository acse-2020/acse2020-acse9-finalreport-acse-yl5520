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

    fig, axs = plt.subplots(2, 2, figsize=(16, 12))
    for ax, urt, lrt in zip(axs.flatten(), URTs, LRTs):
        ax.semilogy(t_list, urt, label='URT')
        ax.semilogy(t_list, lrt, label='LRT')
        ax.plot([-5, 35], [1e2, 1e2], ':')
        ax.set_xlim(-1, 31)
        ax.set_ylim(1, 1e14)
        ax.set_xlabel('time (day)')
        ax.set_ylabel('viral load')
        ax.legend(loc='best')
    axs.flatten()[0].set_title('patient')
    axs.flatten()[1].set_title('inhale=1.')
    axs.flatten()[2].set_title('inhale=.5')
    axs.flatten()[3].set_title('inhale=1.e-3')
    plt.tight_layout()
    plt.savefig('tmp.png')
