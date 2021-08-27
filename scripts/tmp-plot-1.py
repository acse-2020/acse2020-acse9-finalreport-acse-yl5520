#!/usr/bin/env python3
# encoding: utf-8
# versions: 3.8.2

# third-party modules
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


values = np.array([
    [ 14, 141, 150, 252],
    [ 77, 217, 230, 252],
    [ 21, 120, 137, 252],
    [ 77, 208, 224, 252],
    [ 91, 221, 240, 252],
    [164, 218, 230, 252],
    [174, 230, 242, 252],
])
label = ['Total win', 'URT win', 'LRT win', 'Lose']

fig, ax = plt.subplots(1, 1, figsize=(8, 6))
idx = np.arange(7)
for i in range(3, -1, -1):
    ax.bar(idx, values[:, i], .5, label=label[i])
ax.set_xlim(-1, 9)
ax.set_xticks(list(range(7)))
ax.set_xticklabels([
    f'type {i} in type {j}' for i, j in zip(
        [1, 1, 1, 1, 3, 5, 7],
        [0, 2, 4, 6, 6, 6, 6]
    )
], rotation=45, ha='right')
ax.legend(loc='best')
plt.tight_layout()
plt.savefig('images/bar-0.png')
