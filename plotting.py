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
    urt = loadDict['URTs'][0]
    lrt = loadDict['LRTs'][0]
    t_list = loadDict['t_list']

    fig, ax = plt.subplots(1, 1, figsize=(24, 18))
    ax.semilogy(t_list, urt, label='URT')
    ax.semilogy(t_list, lrt, label='LRT')
    ax.plot([-5, 35], [1e2, 1e2], ':')
    ax.set_xlim(-1, 31)
    ax.set_ylim(1, 1e14)
    ax.set_xlabel('time (day)')
    ax.set_ylabel('viral load')
    ax.legend(loc='best')
    plt.tight_layout()
    plt.savefig('tmp.png')
