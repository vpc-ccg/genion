

import collections as co

import numpy as np
import sys
import argparse as ap
#resource.setrlimit(resource.RLIMIT_STACK, [0x10000000, resource.RLIM_INFINITY])
sys.setrecursionlimit(10000)
class kmer_machine:
    def __init__(self, K, file_path, merged=True, poly_t_range = (30,70), poly_t_flank = 4, train_mode=True, kmer_left_ratio=0.5):
        self.K = K
        self.fs = open(file_path, 'r')
        self.buffer_index = 0
        self.line_index = K
        self.merged = merged
        self.poly_t_range = poly_t_range
        self.poly_t_flank = poly_t_flank
        self.train_mode = train_mode
        self.buffer = []
        self.label_buffer = co.deque()
        self.label = np.zeros(1)
        self.line_size = 0
        self.new_read = False
        self.get_next_read()
        self.left_bases = int(K * kmer_left_ratio)
    def get_next_read(self):
        self.cr_rname = self.fs.readline() #rid
        self.line_buffer = self.fs.readline()
        if not self.line_buffer:
            self.fs.close()
            print("lol")
            return False

        self.line_buffer = self.line_buffer.rstrip()
        self.line_index = 0
        

        poly_t_begin = 0
        poly_t_flank = self.poly_t_flank
        poly_t_end   = poly_t_flank
        poly_T = "T" * poly_t_flank
        read = self.line_buffer
        while (poly_t_begin != -1 and
                poly_t_begin < self.poly_t_range[1] and 
                read.count('T',poly_t_begin, poly_t_begin + 3 * poly_t_flank) < 2 *  poly_t_flank): #TODO: Is this good, looking for 2/3 T's in poly-T
            poly_t_begin = read.find(poly_T, poly_t_begin + 1)
            poly_t_end   = poly_t_begin + self.poly_t_flank

        if self.train_mode and (poly_t_begin < self.poly_t_range[0] or poly_t_begin > self.poly_t_range[1]):
            #print(poly_t_begin)
            return self.get_next_read()
        self.line_size = len(read)
        self.label = np.zeros(self.line_size, dtype=int)
        poly_t_end = poly_t_begin + poly_t_flank
        while read.count('T',poly_t_end - poly_t_flank, poly_t_end + poly_t_flank) > poly_t_flank:
            poly_t_end +=1
        #poly_t_end +=window_flank
        self.label[poly_t_begin:poly_t_end] = 2
        self.label[:poly_t_begin] = 1
        
        if not self.merged:
            self.buffer = [str((ord(c) & 6) >> 1) for c in self.line_buffer[:self.K]]
            self.buffer_index = 0
            self.line_index = self.K
            self.label_buffer = co.deque(self.label[:self.K])
        if len(self.label_buffer) == 0:
            self.line_index = self.K
            self.buffer_index = 0
            self.buffer = [str((ord(c) & 6) >> 1) for c in self.line_buffer[:self.K]]
            self.label_buffer = co.deque(self.label[:self.K])
        self.fs.readline() #+
        self.fs.readline() #QUAL
        return True
    def get_kmer(self):

#   self.label[self.line_index-1]  #
        return (str(self.label_buffer[self.left_bases]), self.buffer[self.buffer_index:] + self.buffer[:self.buffer_index])

    def move(self):
        if self.buffer_index >= self.K:
            self.buffer_index=0
        if self.line_index >= self.line_size:
            #print(self.line_index, self.line_size, self.buffer_index, sep="\t")
            has_read = self.get_next_read()
            if not has_read:
                print("lol2")
                return False
            self.new_read = True
        else:
            self.new_read = False
        self.buffer[self.buffer_index] = str((ord(self.line_buffer[self.line_index]) & 6) >> 1)
        self.buffer_index+=1
        self.line_index+=1
        self.label_buffer.append(self.label[self.line_index-1])
        self.label_buffer.popleft()
        return True
        ##TODO if line is done, read until next seq line in fastq
        ## Read one character, put at buffer[index], increment index 

    

def main():
    parser = ap.ArgumentParser()
    parser.add_argument("input", help="Path to reads.fq")
    parser.add_argument("output", help="Path to output.fq")
    parser.add_argument("-v","--validate", action="store_true") ## TODO
    parser.add_argument("-r","--validate_ratio", type=float, default=0.1) ## TODO
    parser.add_argument("-t","--type", choices = ["train", "test"])
    parser.add_argument("-m","--merge", action='store_true', default=True)
    parser.add_argument("-f","--flank", type=int, default=16)
    parser.add_argument("-pi","--poly_a_beg", type=int, default=30)
    parser.add_argument("-pa","--poly_a_end", type=int, default=70)
    parser.add_argument("-pf","--poly_a_flank", type=int, default=4)
    parser.add_argument("-n","--first_n_base", type=int, default=-1)
    parser.add_argument("-kl","--kmer_left_ratio", type=float, default=0.5)
    args = parser.parse_args()
    if args.kmer_left_ratio == 1:
        args.kmer_left_ratio = 0.9999
    input_fastq = args.input
    output_tsv  = args.output
    for_train = args.type == "train"
    index_file = output_tsv + ".index"


    K_MER_FLANK = args.flank
    FIRST_N_BASE = args.first_n_base
    min_beg = args.poly_a_beg
    max_beg = args.poly_a_end

    poly_t_flank = args.poly_a_flank
    poly_T = 'T' * poly_t_flank
    
    index_ptr = 0
    
    km = kmer_machine(args.flank,input_fastq,
            merged=args.merge,
            poly_t_range=(min_beg,max_beg), 
            poly_t_flank=poly_t_flank,
            train_mode=for_train,
            kmer_left_ratio=args.kmer_left_ratio)

    current_base = 0
    with open(output_tsv, 'w') as dat, open(index_file, 'w') as hind:
        b = km.get_kmer()
        current_read = km.cr_rname

        while km.move():
            if km.new_read:
                print(current_base,current_read.split()[0],sep="\t",file=hind)
                current_read = km.cr_rname
            l,k = b
            print(l,"\t".join(k), sep="\t", file = dat)
            b = km.get_kmer()
            current_base+=1

        print(current_base,current_read.split()[0],sep="\t",file=hind)
#        
#        print(current_base,km.cr_rname.split()[0],sep="\t",file=hind)
        print(b)

    return 0
""" # Old Version #
    with open(input_fastq, 'r') as reads, open(output_tsv, 'w') as dat, open(index_file, 'w') as hind:

        line = reads.readline()
        while line:
            rid = line.rstrip()
            ridids = rid.split()
            read = reads.readline()

            label = np.zeros(len(read))
            poly_t_begin = read.find(poly_T)
            if for_train and (poly_t_begin < min_beg or poly_t_begin > max_beg):
                line = reads.readline()
                line = reads.readline()
                line = reads.readline()
                continue
            poly_t_end = poly_t_begin + poly_t_flank
            while read.count('T',poly_t_end - poly_t_flank, poly_t_end + poly_t_flank) > poly_t_flank:
                poly_t_end +=1
            #poly_t_end +=window_flank
            label[poly_t_begin:poly_t_end] = 2
            label[:poly_t_begin] = 1

            enc = [str((ord(c) & 6) >> 1) for c in read]

            for i, l in enumerate(label):
                if i < K_MER_FLANK:
                    continue
                print(int(l), end="\t", file=dat)
                print("\t".join(enc[i-K_MER_FLANK:i+K_MER_FLANK]),file=dat)
                index_ptr+=1
            print(index_ptr, ridids[0], sep="\t", file=hind)
            line = reads.readline()
            line = reads.readline()
            line = reads.readline()
    return 0
"""
if __name__ == "__main__":
    exit(main())
