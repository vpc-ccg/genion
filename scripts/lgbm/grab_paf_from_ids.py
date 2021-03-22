import sys


ids = {}
with open(sys.argv[1],"r") as id_file:
    for line in id_file:
        fields = line.rstrip().split("\t")
        ids[fields[0]] = 1


with open(sys.argv[2],"r") as fq_file, open(sys.argv[3], "w") as out_file:
    line = fq_file.readline()
    while line:
        qname = line.rstrip().split()[0]
        if(qname in ids):
            print(line.rstrip(),file=out_file)
        line = fq_file.readline()
