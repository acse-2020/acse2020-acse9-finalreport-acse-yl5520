#!/usr/bin/env python3
# encoding: utf-8
# versions: 3.8.2

# standard modules
from random import gauss


mean = 1.04e-5
std = 1.35e-5
for _ in range(10):
    print(gauss(mean, std))
