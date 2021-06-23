import sys

import gzip

def low_complexity_sequence(seq, start, end, max_percent = 0.65):
    counts = {}
    for i in range(start,end):
        if seq[i] not in counts:
            counts[seq[i]] = 0
        counts[seq[i]]+=1
    max_count = (end - start) * max_percent
    for s,c in counts.items():
        if c > max_count:
            return True
    return False

def main(argc, argv):
    fastq_file = argv[1]
    chain_file = argv[2]
    out_file = argv[3]
    rid_to_blocks = {}
    with open(chain_file, 'r') as hand:
        line = hand.readline()
        while line:
            fields = line.split("\t")
            rid = fields[0]
            block_count = int(fields[1])
            rid_to_blocks[rid] = []
            for i in range(block_count):
                line = hand.readline()
                rid_to_blocks[rid].append(line.rstrip())
                
            line = hand.readline()

    
    opener = gzip.open if ".gz" in fastq_file else open
    
    with opener(fastq_file, 'rt') as hand, open(out_file, 'w') as ohand:
        line = hand.readline()
        while line:
            fields = line.split()
            if fields[0][1:] not in rid_to_blocks:
                hand.readline()
                hand.readline()
                hand.readline()
                line = hand.readline()
            else:
                seq = hand.readline()
                blocks = rid_to_blocks[fields[0][1:]]
                blocks_to_print = []
                gene_set = set()
                for i,b in enumerate(blocks):
                    bf = b.split("\t")
                    start = int(bf[4])
                    end = int(bf[5])
                    if not low_complexity_sequence(seq,start,end):
                        blocks_to_print.append(i)
                        gene_set.add(bf[11])
                        #print(b.rstrip(), file=ohand)
                if len(blocks_to_print) > 0 and len(gene_set) > 1:
                    print(fields[0][1:],len(blocks_to_print),sep="\t",file=ohand)
                    for i in blocks_to_print:
                        print(blocks[i].rstrip(),file=ohand)
                hand.readline()
                hand.readline()
                line = hand.readline()

if __name__ == "__main__":
    exit(main(len(sys.argv),sys.argv))
