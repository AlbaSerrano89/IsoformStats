#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  5 16:14:47 2018

@author: aserrano
"""

import argparse
import Functions

parser = argparse.ArgumentParser()

parser.add_argument("statsdir", help = "A folder where all the tissue statistics files are located.")
parser.add_argument("minexp", help = "A threshold used to put the least expressed isoform as 'other'.")
parser.add_argument("minsamps", help = "A threshold to determine the minimum of samples to take a gene into account.")

args = parser.parse_args()

statsdir = args.statsdir
minexp = args.minexp
minsamps = args.minsamps

minexp = float(minexp)
minsamps = int(minsamps)

Functions.all_tissues_barplot(statsdir, minexp, minsamps)