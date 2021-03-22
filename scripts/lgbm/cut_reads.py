

import os
import warnings
import pickle
import sys
import argparse
import numpy as np 


def get_next_pos(index_hand,rid):
    line = index_hand.readline()
    line = line.rstrip()
    try:
        sep = line.index("\t")
    except Exception:
        print("linebad",line,sep="\t",file=sys.stderr)
        exit(-1)
    if line[sep+1:]!= rid:
        print(line[sep+1:], rid, sep="\t")
        return -1
    return int(line[:sep])

def process_fq_2(seq_file, pred_file, index_file, outfile, min_polyA_probability = 0.5, min_cut=100, nice = 10, flank = 16):
    with open(seq_file, 'r') as hand, open(pred_file, 'r') as pand, open(index_file, 'r') as iand, open(outfile, 'w') as whand:
        cpos = 0
        line = hand.readline()
        while line:
            rid_fields = line.rstrip().split()
            rid = rid_fields[0]
            line = hand.readline()
            seq = line.rstrip()
            line = hand.readline()
            sep = line.rstrip()
            line = hand.readline()
            qual = line.rstrip()
            
            npos = get_next_pos(iand,rid)
            if npos == -1:
                line = hand.readline()
                continue
            rlen = len(seq)
            weights = np.ones((rlen,3))
            
         #   i = flank -1
            i = 0
            print(rlen,cpos,npos,weights.shape,sep='\n')
            while cpos < npos:
                wel = pand.readline()
                ww = [float(x) for x in wel.split("\t")]
                for j,w in enumerate(ww):
                    weights[i,j] = w
                i+=1
                cpos+=1
            #cpos+=1
            cuts = [ 0 ]
            for i in range(flank, weights.shape[0]):
                if weights[i, 2] > min_polyA_probability:
                    j = i
                    flag = False
                    while weights[j,2] > min_polyA_probability:
                        j-=1
                    cnice = nice
                    while cnice > 0:
                        while  weights[j,1] > min_polyA_probability:
                            j-=1
                            flag = True
                            cnice = 0
                        cnice-=1
                    if cuts[-1] !=j and j - cuts[-1] > min_cut  and flag:
                        cuts.append(j)
           
            cuts.append(rlen)


            for i,x in enumerate(zip(cuts[:-1],cuts[1:])):
                a, b = x
                if len(cuts) > 2:
                    print("{}-{}_{} {}".format(rid,a,b," ".join(rid_fields[1:])), file=whand)
                else:
                    print("{} {}".format(rid," ".join(rid_fields[1:])), file=whand)
                print(seq[a:b], file=whand)
                print(sep, file=whand)
                print(qual[a:b], file=whand)
            line = hand.readline()

def process_fq(seq_file, pred_file, index_file, outfile, min_polyA_probability = 0.5, min_cut=100, nice = 10, flank = 16):
    with open(seq_file, 'r') as hand, open(pred_file, 'r') as pand, open(index_file, 'r') as iand, open(outfile, 'w') as whand:
        cpos = 0
        line = hand.readline()
        while line:
            rid_fields = line.rstrip().split()
            rid = rid_fields[0]
            line = hand.readline()
            seq = line.rstrip()
            line = hand.readline()
            sep = line.rstrip()
            line = hand.readline()
            qual = line.rstrip()
            
            npos = get_next_pos(iand,rid)
            if npos == -1:
                line = hand.readline()
                continue
            rlen = len(seq)
            weights = np.ones((rlen,3))
            
         #   i = flank -1
            i = 0
            print(rlen,cpos,npos,weights.shape,sep='\n')
            while cpos < npos:
                wel = pand.readline()
                ww = [float(x) for x in wel.split("\t")]
                for j,w in enumerate(ww):
                    weights[i,j] = w
                i+=1
                cpos+=1
            #cpos+=1
            cuts = [ 0 ]
            for i in range(flank, weights.shape[0]):
                if weights[i, 2] > min_polyA_probability:
                    j = i
                    flag = False
                    while weights[j,2] > min_polyA_probability:
                        j-=1
                    cnice = nice
                    while cnice > 0:
                        while  weights[j,1] > min_polyA_probability:
                            j-=1
                            flag = True
                            cnice = 0
                        cnice-=1
                    if cuts[-1] !=j and j - cuts[-1] > min_cut  and flag:
                        cuts.append(j)
           
            cuts.append(rlen)


            for i,x in enumerate(zip(cuts[:-1],cuts[1:])):
                a, b = x
                if len(cuts) > 2:
                    print("{}-{}_{} {}".format(rid,a,b," ".join(rid_fields[1:])), file=whand)
                else:
                    print("{} {}".format(rid," ".join(rid_fields[1:])), file=whand)
                print(seq[a:b], file=whand)
                print(sep, file=whand)
                print(qual[a:b], file=whand)
            line = hand.readline()

def process_fa(seq_file, model):
    raise NotImplementedError

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Path to reads.fq")
    parser.add_argument("index", help="index of reads.fq")
    parser.add_argument("predict", help="Path to lgbm predictions")
    parser.add_argument("output", help="Path to output.fq")
    args = parser.parse_args()

    seqf = args.input
    predict = args.predict
    index = args.index
    if ".fq" in seqf or ".fastq" in seqf:
        process_fq(seqf, predict, index, args.output)
    elif ".fa" in seqf:
        process_fa(seqf, model)
    else:
        print("I don't know this file format: " + seqf , file=sys.stderr)
        return -1

    return 0


if __name__ == "__main__":
    exit(main())
