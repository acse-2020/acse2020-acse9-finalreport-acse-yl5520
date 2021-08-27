#!/usr/bin/env python3
# encoding: utf-8
# versions: 3.8.2

# standard modules
import pickle

# third-party modules
import matplotlib.pyplot as plt
from tqdm import tqdm


data_file = 'data/unmask-unvaccine-cmp-uninfect-and-infect.pkl'
with open(data_file, 'rb') as fin:
    dic = pickle.load(fin)

ls = [
    'T1-URT', 'T2-URT', 'I-URT', 'V-URT',
    'T1-LRT', 'T2-LRT', 'I-LRT', 'V-LRT',
    'V-ENV'
]
ts = dic['t']
data = dic['data']
title = dic['title']
choice = [0, 1]

fig, axs = plt.subplots(2, 1, figsize=(8, 10))
for ax, i in tqdm(list(zip(axs.flatten(), choice))):
    for j in range(len(ls)):
        ax.plot(ts, data[:, i, j], label=f'{ls[j]}')
    ax.set_yscale('symlog')
    ax.set_ylim(bottom=1e-9)
    ax.set_xlabel('time (day)', size=16)
    ax.set_xticks([i for i in range(int(ts[-1]) + 1)])

plt.tight_layout()
plt.savefig('images/tmp')
