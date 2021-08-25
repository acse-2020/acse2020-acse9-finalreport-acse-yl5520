#!/usr/bin/env python3
# encoding: utf-8
# versions: 3.8.2

# standard modules
import pickle

# third-party modules
import matplotlib.pyplot as plt
from tqdm import tqdm


data_file = 'data/9hosts-10vent-weekly-interact-30day-nounmaskinfect.pkl'
with open(data_file, 'rb') as fin:
    dic = pickle.load(fin)

ls = [
    'T1-URT', 'T2-URT', 'I-URT', 'V-URT',
    'T1-LRT', 'T2-LRT', 'I-LRT', 'V-LRT',
    'V-ENV'
]
ts = dic['t']
data = dic['data']
nums = [1, 4]

titles = ['effective infection', 'uneffective infection']

fig, axs = plt.subplots(1, 2, figsize=(16, 5))
for ax, num, title in tqdm(zip(axs.flatten(), nums, titles)):
    for j in [0, 1, 2, 3, 8]:
        ax.plot(ts, data[:, num, j], label=f'{ls[j]}')
    ax.set_yscale('symlog')
    ax.set_ylim(bottom=1e-9)
    ax.set_xlabel('time (day)')
    ax.set_xticks([i for i in range(int(ts[-1]) + 1)])
    ax.set_title(title)
    ax.legend(loc='best')

plt.tight_layout()
plt.show()
