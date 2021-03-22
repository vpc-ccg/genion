

import collections as co
import numpy as np
import sys
import argparse as ap



def main():
    parser = ap.ArgumentParser()

    HHHH = 0.0001
    ops = {
        "gte": lambda x,y: x >= y,
        "gt":  lambda x,y: x > y,
        "lte": lambda x,y: x <= y,
        "lt":  lambda x,y: x < y,
        "eq":  lambda x,y: abs(x-y) < HHHH
    }

    parser.add_argument("ops", choices=list(ops.keys()), help="Comparison operation")
    parser.add_argument("alig", help="Path to alig.paf, use minimap2 paf with cigar, only primary aligs")
    parser.add_argument("output", help="Path to output paf")
    parser.add_argument("-i","--identity-ratio", type=float, default=0.90)
    
    args = parser.parse_args()

    op = ops[args.ops]
    target = args.identity_ratio

    if target == 1:
        target -= HHHH

    with open(args.alig, 'r') as aand, open(args.output, 'w') as oand:
        for line in aand:
            fields = line.rstrip().split("\t")
            matches = float(fields[9])
            qlen    = float(fields[1])
            if op(matches/qlen, target):
                print(line.rstrip(), file=oand)


if __name__ == "__main__":
    exit(main())
