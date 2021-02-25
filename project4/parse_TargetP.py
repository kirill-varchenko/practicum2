#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json
import sys

if len(sys.argv) < 2:
    sys.exit(0)

with open(sys.argv[1], 'r') as fi:
    data = json.load(fi)

for k in data['SEQUENCES']:
    print(k, data['SEQUENCES'][k]['Prediction'], sep='\t')