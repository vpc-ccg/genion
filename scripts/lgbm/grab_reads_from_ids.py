import sys
import gzip

ids = {}
with open(sys.argv[1],"r") as id_file:
    for line in id_file:
        fields = line.rstrip().split("\t")
        ids[fields[0]] = 1


opener = gzip.open if ".gz" in sys.argv[2] else  open
#with opener( sys.argv[1], 'rt' if ".gz" in sys.argv[1] else 'r') as hand:
with opener(sys.argv[2],"rt" if ".gz" in sys.argv[2] else "r") as fq_file, open(sys.argv[3], "w") as out_file:
    line = fq_file.readline()
    while line:
        qname = line.rstrip().split()[0]
        if(qname[1:] in ids):
            print(line.rstrip(),file=out_file)
            for i in range(3):
                line = fq_file.readline()
                print(line.rstrip(),file=out_file)
        else:
            for i in range(3):
                line = fq_file.readline()
        
        line = fq_file.readline()
