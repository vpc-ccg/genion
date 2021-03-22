#!/usr/bin/env python3
import argparse
import re
from itertools import groupby,chain
import sys
rev_comp = dict(
    A='T',
    C='G',
    G='C',
    T='A'
)
dna_id = dict(
    A='A',
    C='C',
    G='G',
    T='T'
)

def find_longest_poly(seq, match_score=1, mismatch_score=-2, char='A'):
    if len(seq)==0:
        return
    if seq[0]==char:
        scores=[match_score]
    else:
        scores=[0]
    for m in (match_score if c==char else mismatch_score for c in seq[1:]):
        scores.append(max(0,scores[-1]+m))
    for k,g in groupby(enumerate(scores),lambda x:x[1]>0):
        if not k:
            continue
        i,s = list(zip(*g))
        max_s,max_i=max(zip(s,i))
        l = max_i+1-i[0]
        yield i[0], l, seq[i[0]:i[0]+l].count(char)/l


def main(args = []):
    valid_ids = set()

    with open(args[1], 'r') as hand:
        for line in hand :
            rid, seq = line.rstrip().split("\t")
            s_polys = []

            for char in ['A','T']:
                for i,l,p in find_longest_poly(seq[0:], char=char):
                    if l < 20 or p < 0.85:
                        continue
                    s_polys.append((i,l,p,char))
            if len(s_polys) > 0:
                #sorted_polys = sorted(s_polys,key = lambda x:x[0])

                i,l,p,char = max(s_polys, key = lambda x: x[2])
                print(rid, char, i,l,p, sep="\t")
            else:
                print(rid, "NONE",sep="\t")


if __name__ == "__main__":
    exit(main(sys.argv))
