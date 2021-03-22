


import sys
import gzip




valid_ids = set()

opener = gzip.open if ".gz" in sys.argv[2] else  open

with opener(sys.argv[2], 'rt' if ".gz" in sys.argv[2] else 'r') as hand:
    for line in hand:
        rid = line.rstrip().split("\t")[0]
        valid_ids.add(rid)


opener = gzip.open if ".gz" in sys.argv[1] else  open
with opener( sys.argv[1], 'rt' if ".gz" in sys.argv[1] else 'r') as hand:
    line = hand.readline()

    while line:
        rid = line.split()[0][1:]
        if rid not in valid_ids:
            line = hand.readline()
            line = hand.readline()
            line = hand.readline()
            line = hand.readline()
            continue
        line = hand.readline()
        print(rid,line.rstrip(),sep="\t")


        line = hand.readline()
        line = hand.readline()
        line = hand.readline()
