# Genion
## Installation


## Installation from Source
```bash
git clone https://github.com/vpc-ccg/genion
cd genion
make
```

## Prerequisites

Can be installed using conda/mamba
```bash
conda create --file genion.env --name genion-env
conda activate genion-env
```

|Dependencies | Version |
|-------- |-----|
|Python   | 3.x |
|c++ | std++17 |
|[deSALT](https://github.com/ydLiu-HIT/deSALT) | >= 1.5.5 |
|[minimap2](https://github.com/lh3/minimap2/tree/master/misc) | >= 2.17 |
|[paftools](https://github.com/lh3/minimap2/tree/master/misc) |  |
|[snakemake](https://snakemake.readthedocs.io/en/stable/) | >= 5.3.0 |


lightgbm is required for optional chimeric read correction step.

## Running Genion

```bash
    ./fusion run
        --reference /path/to/ref/fasta
        --gtf       /path/to/annot/gtf
        --gpaf      /path/to/genomic/mapping/paf 
        -s          /path/to/gene/homology/tsv
        -o          /path/to/output/tsv
        -d          /path/to/genomicSuperDups.txt
        -i          /path/to/input/fastq
```

## Required References 

GTF annotation and Whole genome reference sequence can be downloaded from https://uswest.ensembl.org/info/data/ftp/index.html

genomicSuperDups.txt can bve downloaded from  ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz (Should be extracted using gzip)

Homology tsv file can be produced using ENSEMBL cDNA reference and following:

```bash
        minimap2 [cdna.fa] [cdna.fa] -X -t [threads] -2 -c -o [cdna.selfalign.paf]
        cat [cdna.selfalign.paf] | cut -f1,6 | sed 's/_/\t/g' | awk 'BEGIN{OFS=\"\\t\";}{print substr($1,1,15),substr($2,1,15),substr($3,1,15),substr($4,1,15);}' | awk '$1!=$3' | sort | uniq > [cdna.selfalign.tsv]
```

# Genion Snakemake
We provide a snakemake file to help running genion. This snakemake also assists two step deSALT mapping.
## Project Configuration
In order to run Genion, you need to create a project configuration file namely ``config.yaml``. 
This configuration consists of a number mandatory settings and some optional advance settings. 
Below is the list of the all the settings that you can set in your project.

|config-paramater-name | Type | Description|
|------------------------------|-----------|--------------------------------------------------------------------------------------------------------------------------------------|
| path                          | Mandatory | Full path to project directory.  |
| reference-dna                 | Mandatory | Full path to the dna reference                |
| reference-cdna                | Mandatory | Full path to the cdna reference                |
| annotation-gtf                 | Mandatory | Full path to the gtf annotation                |
| rawdata-base                  | Mandatory | Location of the input fastq files relative to ``path``.                                                         |
| or rawdata                       | Mandatory | Full path to the location of the input fastq files                                                        |
| input                        | Mandatory | A list of input files per sample. See the following example     |
| ext                          | Optional  | (Important, if this is wrong snakemake won't run correctly) extension of the fastq files used in input (``fastq``,``fastq.gz``,``fq``,``fq.gz``) default:``fastq`` |
| genion-binary                | Optional  | Path to genion binary, should be set if genion is not in $PATH ||
| desalt-index                 | Optional  | If not provided, reference will be indexed on the run ||
| analysis-base                | Optional  | Location of intermediate files relative to ``path``. default: ``{path}/analysis``|
| or analysis                  | Optional  | Full path to the location of intermediate files. default: ``{path}/analysis``|
| results-base                 | Optional  | Location of final results relative to the ``path``. default: ``{path}/results``  |
| or results                   | Optional  | Full path to the location of final results. default: ``{path}/results``  |
| wg-aligner                   | Optional  | Mapper to use (``deSALT``, ``minimap2``) default: ``deSALT``                                                                      |

### Input formatting in the config file
Each input requires a fastq file and type. Type is used to configure parameters by the mapper.
Following are the available types of input:

|type   |   Technology                                  |
|-------|-----------------------------------------------|
|ccs    |   PacBio SMRT CCS reads: error rate 1%        |
|clr    |   PacBio SMRT CLR reads: error rate 15%       |
|ont1d  |   Oxford Nanopore 1D reads: error rate > 20%  |
|ont2d  |   Oxford Nanopore 2D reads: error rate > 12%  |

The following a an example of ``config-yaml`` for Nanopore and Pacbio runs for a sample
```yaml
path:
    /path/to/project/directory

annotation-gtf/:
    /path/to/annotation/gtf
reference-dna:
    /path/to/reference/dna
reference-cdna:
    /path/to/reference/cdna
desalt-index:
    /path/to/index/dir/
ext:
    fastq.gz
wg-aligner:
    deSALT
input:
    "A_clr":
        type:
            clr
        fastq:
            - A_clr.fastq.gz
    "A_ont":
        type:
            ont1d
        fastq:
            - A_ont.fastq.gz
```

# Simulated Dataset

Simulated sequences are uploaded in 4 gzipped parts.
simulation.reference.tar.gz contains reference of the simulated sequences in the simulation. Details are written in the README document.

## Download
Simulated gene fusion dataset can be downloaded from:
https://figshare.com/articles/dataset/Simulated_RNA_Long_Reads_with_ONT_error_profile_and_gene_fusions_/14265554


## Setup
```bash
zcat simulation.fastq.part00.gz  simulation.fastq.part01.gz  simulation.fastq.part02.gz  simulation.fastq.part03.gz > simulation.fastq
tar xzvf simulation.reference.tar.gz 
```







