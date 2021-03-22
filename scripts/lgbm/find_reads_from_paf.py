

import collections as co
import numpy as np
import sys
import argparse as ap


def use_dict_filter(args, is_fastq):
    with open(args.input, 'r') as hand, open(args.alig, 'r') as aand, open(args.output, 'w') as oand:
        read_ids = set()
        for line in aand:
            first_space = line.find("\t")
            rid = line[:first_space]
            read_ids.add(rid)

        line = hand.readline()
        while line:
            fields = line.split()
            rid = fields[0][1:]
            if rid in read_ids:
                print(line.rstrip(), file=oand)
                for i in range(3):
                    line = hand.readline()
                    print(line.rstrip(), file=oand)
            else:
                for i in range(3):
                    line = hand.readline()

            line = hand.readline()

def same_order_filter(args):
    pass

def main():
    parser = ap.ArgumentParser()
    parser.add_argument("input",  help="Input.fastq")
    parser.add_argument("alig", help="Path to alig.paf, use minimap2 paf with cigar, only primary aligs")
    parser.add_argument("output", help="Output.fastq")
    parser.add_argument("-s", "--same_order", action="store_true")
    
    args = parser.parse_args()

    is_fastq = args.input.endswith(".fastq") or args.input.endswith(".fq")
    if not args.same_order:
        return use_dict_filter(args,is_fastq)
    else:
        return same_order_filter(args,is_fastq)

if __name__ == "__main__":
    exit(main())
