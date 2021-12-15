

import sys



################# UTIL FUNCTIONS ################

cfg_text = ""
def tee( st, file=sys.stderr,end="\n"):
    global cfg_text
    cfg_text = cfg_text + st + end
    print(st, file=file, end=end)

def get_tool_path(name):
    from shutil import which
    return which(name)


def assert_input_tool(tool):
    tp = get_tool_path(tool)
    assert tp is not None, tool + " not found in the PATH"
    tee("Using {} from {}".format(tool,tp), file=sys.stderr)
    return tp 

def assert_tool(tool):
    tp = get_tool_path(tool)
    assert tp is not None, tool + " not found in the PATH"
    tee("Using {} from {}".format(tool,tp), file=sys.stderr)

def cfg_default( key, val):
    if key not in config:
        config[key] = val
    tee( "{}: {}".format(key,config[key]), file=sys.stderr)


def cfg_mandatory( key):   
    assert key in config, "{} is a mandatory field in config".format(key)
    tee( "{}: {}".format(key,config[key]), file=sys.stderr)
    
def cfg_optional( key):
    val = "[unset]"
    if key in config:
        val = config[key]
    tee( "{}: {}".format(key,value), file=sys.stderr)

################# UTIL FUNCTIONS ################

no_gtf_bias = True

tee("### RUN CONFIGURATION ###", file =sys.stderr)
cfg_mandatory("path")

cfg_mandatory("reference-dna")
cfg_mandatory("reference-cdna")
cfg_mandatory("annotation-gtf")
cfg_mandatory("input")

cfg_default("wg-aligner", "minimap2")
cfg_default("rawdata-base","/raw-data")
cfg_default("analysis-base","/analysis")
cfg_default("results-base","/results")
assert_tool("minimap2")
assert_tool("paftools.js")
assert_tool("deSALT")

cfg_default("chimeric-correction", False)
cfg_default("minimap_threads", 63)


cfg_default("raw-data", config["path"] + config["rawdata-base"])
cfg_default("analysis", config["path"] + config["analysis-base"])
cfg_default("results", config["path"] + config["results-base"])

cfg_default("ext", "fastq.gz")
ext  = config["ext"]
tee("#########################", file =sys.stderr)


path_names = {
    "linkeddata"  : "linked-data",
    "mapping"     : "mapping",
    "fusion"      : "fusion",
    "annot"      : "annotate",
    "lgbm"        : "lgbm",
    "cutting"     : "cutting",
    "plots"     : "plots",
}


index = 0
for k,v in path_names.items():
    path_names[k] = config["analysis"] + "/{:0>3d}-{}".format(index,v)
    index+=1



onsuccess:
    cfg_log_path = config["analysis"] + "/config.log"
    with open(cfg_log_path, 'w') as hand:
        print(cfg_text, file=hand)
onerror:
    cfg_log_path = config["analysis"] + "/config.log"
    with open(cfg_log_path, 'w') as hand:
        print(cfg_text, file=hand)

def get_fq_name(wildcards):
    fq_name = config["input"][wildcards.sample]["fastq"][0]
    return config["raw-data"] +"/"+fq_name

def get_sample_name(wildcards):
    assert wildcards.id in config["input"], "Sample {} not in config!".format(wildcards.id)
    return config["input"][wildcards.id]["sample"]

if "chimeric-correction" in config and config["chimeric-correction"] == True:
    rule all_cut:
        input:
            expand( path_names["cutting"] + "/{sample}/candidate.cut.reads.fq" ,sample=config["input"].keys()),
            expand(config["results"] + "/{sample}.fusions.tsv",sample=config["input"].keys())
            #expand( path_names["annot"] + "/{sample}/annotation.tsv" ,sample=config["input"].keys())


rule all_results:
    input:
        expand(config["results"] + "/{sample}.fusions.tsv",sample=config["input"].keys())


rule move_to_results:
    input:
        fusions = path_names["fusion"] + "/{sample}.annot.tsv",
    output:
        fus = config["results"] + "/{sample}.fusions.tsv",
        rth = config["results"] + "/{sample}.readthrough.tsv", 
    run:
        with open(input[0], 'r') as ihand,\
             open(output.fus, 'w') as gfhand,\
             open(output.rth, 'w') as rthand:
            for hand in [gfhand, rthand]:
                print("ID","SYMBOL","FFIGF","FiN","COUNT","PASSTAG",sep="\t",file=hand)
            for line in ihand:
                if "PASS:GF" in line:
                    print(line.rstrip(), file=gfhand)
                elif "PASS:RT" in line:
                    print(line.rstrip(), file=rthand)


rule download_and_extract_dups:
    output:
        config["analysis"] + "/genomicSuperDups.txt"
    shell:
        "cd "+ config["analysis"] +  " && wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz && gzip -d genomicSuperDups.txt.gz"


rule find_fusion:
    input:
        fq=path_names["linkeddata"]+"/{sample}." + config["ext"],
        gpaf=path_names["mapping"]+"/{sample}.wg.final.paf",
        self_align=path_names["mapping"]+"/cdna.self.tsv",
        dups=config["analysis"]+"/genomicSuperDups.txt",
        binary = "genion" if "genion-binary" not in config else config["genion-binary"],
    output:
        fusions = path_names["fusion"] + "/{sample}.annot.tsv",
        normals = path_names["fusion"] + "/{sample}.annot.tsv.fail",
    params:
        gtf=config["annotation-gtf"],
    threads:
        1
    shell:
        "{input.binary} -i {input.fq} -d {input.dups} --gtf {params.gtf} -g {input.gpaf} -s {input.self_align} -t {threads} -o {output.fusions}"


def get_read_type(wildcards):
    assert wildcards.id in config["input"], "Sample {} not in config!".format(wildcards.id)
    return config["input"][wildcards.id]["type"]
    
if config["wg-aligner"] == "deSALT":
    if no_gtf_bias:

        if "desalt-index" in config:
            rule desalt_index_link:
                input:
                    config["desalt-index"]
                output:
                   folder=directory(config["analysis"] + "/desaltIndex/"),
                params:
                   folder=config["analysis"] + "/desaltIndex/",
                shell:
                    "ln -s {input} {output.folder}"
        else:
            rule desalt_index:
                input:
                   ancient(config["reference-dna"]),
                output:
                   seqfile=config["analysis"] + "/desaltIndex/ref.seq",
                   folder=directory(config["analysis"] + "/desaltIndex/"),
                params:
                   folder=config["analysis"] + "/desaltIndex/",
                shell:
                    "deSALT index {input} {params.folder}"
        rule maptowg_deSALT:
            input:
                index=config["analysis"]+"/desaltIndex/",
                fastq=path_names["linkeddata"] + "/{id}." + ext,
            output:
                sam=path_names["mapping"]+"/{id}.wg.sam"
            params:
                p="-N 10",
                typ=get_read_type,
            threads:
                48
            benchmark:
                path_names["mapping"]+"/{id}.benchmark.txt"
            shell:
                "deSALT aln -x {params.typ} -f {output.sam}.tmp -t {threads} {params.p}  {input.index} {input.fastq} -o {output.sam}"
        rule masked_maptowg_deSALT:
            input:
                index=config["analysis"]+"/desaltIndex/",
                fastq=path_names["mapping"] + "/{id}.masked." + ext,
            output:
                sam=path_names["mapping"]+"/{id}.wg.masked.sam"
            params:
                p="-N 10",
                typ=get_read_type,
            threads:
                48
            benchmark:
                path_names["mapping"]+"/{id}.masked.benchmark.txt"
            shell:
                "deSALT aln -x {params.typ} -f {output.sam}.tmp -t {threads} {params.p}  {input.index} {input.fastq} -o {output.sam}"
    rule samtopaf:
        input:
            "{sample}.sam"
        output:
            "{sample}.paf"
        shell:
            "paftools.js sam2paf -L {input} > {output}"
    rule mask_mapped:
        input:
            paf=path_names["mapping"]+"/{id}.wg.paf",
            fq =path_names["linkeddata"] + "/{id}." + ext,
        output:
            fq =path_names["mapping"]+"/{id}.masked." + ext,
        params:
            mintail=50,
            minfront=120,
        run:
            with open(input.paf, 'r') as hand:
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

                            if front_clip_len > params.minfront:
                                if prev not in reads2clip:
                                    reads2clip[prev] = []
                                reads2clip[prev].append((0,front_clip_len,'f'))
                            if tail_clip_len > params.mintail:
                                if prev not in reads2clip:
                                    reads2clip[prev] = []
                                reads2clip[prev].append((aligs[-1],rlen,'t'))

                        prev = fields[0]
                        aligs = []
                        rlen = int(fields[1])
                       
                    line = hand.readline()
                    aligs.append(int(fields[2]))
                    aligs.append(int(fields[3]))
                
            opener = gzip.open if ".gz" in input.fq else open

            with opener(input.fq,'rt') as hand, open(output.fq, 'w') as oand:
                
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
                        print("@{}/{} {} {}".format(rid,i,clip, len(seq)), file=oand)
                        print("{}{}{}".format('N'*clip[0],seq[clip[0]:clip[1]],'N'*(len(seq)-clip[1])), file=oand)
                        print(plus, file=oand)
                        print(qual, file=oand)
                    line = hand.readline()


    rule merge_pafs:
        input:
            p1=path_names["mapping"]+"/{id}.wg.paf",
            p2=path_names["mapping"]+"/{id}.wg.masked.paf",
        output:
            p1=path_names["mapping"]+"/{id}.wg.final.paf",
        shell:
            "cat {input.p1} <(sed 's/\/[0-9]//' {input.p2}) | sort -k1,1 -k2,2 -k3n,3 -k4n,4 > {output}"
else:
 rule maptowg:### USE INDEX
        input:
            path_names["linkeddata"]+"/{sample}." + ext,
        output:
            path_names["mapping"]+"/{sample}.wg.final.paf"
        params:
            ref=config["reference"]+"/0.dna.3.14.mmi",
            opt="-c -y -x splice -k 12 -w 3 --hard-mask-level -N 100 "
        threads:
            64
        shell:
            "minimap2 {params.opt} -t {threads}  {params.ref}  {input} -o {output}"



"""DEPRECATED
rule maptots:### USE INDEX
    input:
        path_names["linkeddata"]+"/{sample}." + ext,
    output:
        path_names["mapping"]+"/{sample}.ts.paf"
    params:
        ref=config["reference"]+"/reference.cdna.3.14.mmi",
        opt="-c -y -x splice -k 12 -w 3 --hard-mask-level -N 100 "
    threads:
        64
    shell:
        "minimap2 {params.opt} -t {threads}  {params.ref}  {input} -o {output} "
"""

rule self_align:
    input:    
        config["reference-cdna"]
    output:
        config["analysis"] + "/{sample}.self.paf"
    threads:
        64
    shell:
        "minimap2 {input} {input} -X -t {threads} -2 -c -o {output} "


rule self_align_compute: 
    input:
        "{ref}.self.paf"
    output:
        "{ref}.self.tsv"
    shell:
        "cat {input} | cut -f1,6 | sed 's/_/\t/g' | awk 'BEGIN{{OFS=\"\\t\";}}{{print substr($1,1,15),substr($2,1,15),substr($3,1,15),substr($4,1,15);}}' | awk '$1!=$3' | sort | uniq > {output}"

rule link_fq:
    input:
        get_fq_name
    output:
         path_names["linkeddata"]+"/{sample}." + ext,
    shell:
        "ln -s {input} {output}"
