
import sys
import argparse as ap

#Find aggresively clipped reads, remap.


import gzip

def main():
    parser = ap.ArgumentParser()
    parser.add_argument("input", help="Path to reads.fq")
    parser.add_argument("paf", help="Path to reads.paf, sorted by name as it is reported by mapper")
    parser.add_argument("output", help="Path to output.fq")
    parser.add_argument("-q","--min-qual", type=int, default=20, help="Minimum avg quality of the clipped region to remap")

    parser.add_argument("-t","--tail-min-len", type=int, default=100, help="Minimum length of the clipped region tail to remap")
    parser.add_argument("-f","--front-min-len", type=int, default=200, help="Minimum length of the clipped region front to remap")
    args = parser.parse_args()

    reads2clip = {}

    with open(args.paf, 'r') as hand:
        line = hand.readline()
        prev = ""
        aligs = []
        rlen = -1
        while line:
            if "tp:A:P" not in line:
                line = hand.readline()
                continue
            fields = line.rstrip().split("\t")
            if fields[0] != prev:
                if prev != "":
                    aligs = sorted(aligs,key=int)
                    front_clip_len = aligs[0]
                    tail_clip_len  = rlen - aligs[-1]

                    if front_clip_len > args.front_min_len:
                        if prev not in reads2clip:
                            reads2clip[prev] = []
                        reads2clip[prev].append((0,front_clip_len,'f'))
                    if tail_clip_len > args.tail_min_len:
                        if prev not in reads2clip:
                            reads2clip[prev] = []
                        reads2clip[prev].append((aligs[-1],rlen,'t'))

                prev = fields[0]
                aligs = []
                rlen = int(fields[1])
               
            line = hand.readline()
            aligs.append(int(fields[2]))
            aligs.append(int(fields[3]))

    opener = gzip.open if ".gz" in args.input else open

    with opener(args.input,'rt') as hand, open(args.output, 'w') as oand:
        
        line = hand.readline()
        while line:
            for i,c in enumerate(line):
                if c.isspace():
                   
                    rid = line[1:i]
                    break
            to_clip = reads2clip.get(rid, [])

            seq = hand.readline().rstrip()
            plus = hand.readline().rstrip()
            qual = hand.readline().rstrip()
            for i,clip in enumerate(to_clip):
                if True: #avg_qual(qual[clip[0]:clip[1]]) > args.min_qual:

                    print("@{}/{} {} {}".format(rid,i,clip, len(seq)), file=oand)
                    print("{}{}{}".format('N'*clip[0],seq[clip[0]:clip[1]],'N'*(len(seq)-clip[1])), file=oand)
                    print(plus, file=oand)
                    print(qual, file=oand)
            line = hand.readline()

    return 0



if __name__ == "__main__":
    exit(main())
