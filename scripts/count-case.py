#!/usr/bin/env python3
# encoding: utf-8
# versions: 3.8.2

# standard modules
import argparse
import pickle

# third-party modules
from tqdm import tqdm


parser = argparse.ArgumentParser()
parser.add_argument('pkl')
args = parser.parse_args()

with open(args.pkl, 'rb') as fin:
    dic = pickle.load(fin)

data = dic['data']
title = dic['title']
idx = [1, 5]

size = data.shape[1]
total_win = urt_win = lrt_win = 0

for i in range(size):
    if '-infect' in title[i]:
        pass

    state = [False, False]
    for j in range(2):
        state[j] = data[-1, i, idx[j]] != data[-100, i, idx[j]]

    if state[0] and state[1]:
        total_win += 1
    elif state[0]:
        urt_win += 1
    elif state[1]:
        lrt_win += 1

print(f'T: {total_win}, U: {urt_win}, L:{lrt_win}')
