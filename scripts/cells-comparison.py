#!/usr/bin/env python3
# encoding: utf-8
# versions: 3.8.2

# standard modules
import math

# third-party modules
import matplotlib.pyplot as plt
import numpy as np


host_nums = [i for i in range(1, 101)]
square_cells = []
long_cells = []
for num in host_nums:
    # square case
    length = 1 + 2 * math.ceil(math.sqrt(num))
    square_cells.append(3 * length * length)

    length = 1 + 2 * num
    long_cells.append(3 * 3 * length)

fig, ax = plt.subplots(1, 1, figsize=(8, 6))
ax.plot(host_nums, square_cells, label='square total')
ax.plot(host_nums, np.array(square_cells) - np.array(host_nums), label='square empty')
ax.plot(host_nums, long_cells, label='long total')
ax.plot(host_nums, np.array(long_cells) - np.array(host_nums), label='long empty')
ax.set_xlabel('number of students')
ax.set_ylabel('number of empty cells')
ax.set_title('empty cells comparison')
ax.legend(loc='best')
plt.show()
