# Genion
An accurate tool to detect gene fusion from long transcriptomics reads. 

Genion comes with a stand-alone binary and helper Snakemake to assist mapping and preparing reference files.

## Installation
You can install genion through conda, docker or from source.
### Installation with [bioconda](http://bioconda.github.io/recipes/genion/README.html)
```bash
conda install -c bioconda genion
```
Check out [mamba](https://github.com/mamba-org/mamba). A faster conda reimplementation.

### Installation with Docker
```bash
git clone https://github.com/vpc-ccg/genion
cd genion/docker
docker build . -t genion:latest
```
After genion is built with Docker, you can run it with the following command 
```bash
docker run --user=$UID -v /path/to/inputs:/input -v /path/to/outputdir:/output genion [args]
```
### Installation from Source
|Dependencies | Version |
|-------- |-----|
|c++ | gcc >= 9 or clang >= 8|
|zlib| >= 1.2.11 |

```bash
git clone https://github.com/vpc-ccg/genion
cd genion
make
```

## Running Genion


### Input 
Genion requires following input files to run:
* [Mapping file of transcriptomics long reads (paf)](https://github.com/lh3/miniasm/blob/master/PAF.md): Genion does not do mapping. It accepts mappings in paf format. You can use any splice-aware long read to whole genome mapper (and convert sam to paf using paftools if mapper doesn't output paf).
* Long reads file(fast{a,q}): These are used for filter low complexity sequence filtering
* [Gene annotation file (GTF)](https://m.ensembl.org/info/website/upload/gff.html)
* Sequence similarity file: This is used to filter candidates from genes with similar sequences. This file is produced by all to all mapping cDNA reference file with itself. It can be created using genion snakemake or command line given in the Required References section. This file is a tab separated 2 column file containing transcript pairs.
* [Duplication annotation](http://genome.ucsc.edu/cgi-bin/hgTables?hgta_doSchemaDb=hg18&hgta_doSchemaTable=genomicSuperDups): Genomic segmental duplication annotation. This used to filter out candidates that come from copies of the same segmental duplication. For hg38, it can be downloaded from ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz.

### Output
* [output]: Contains called gene fusions and readthrough. It is a tab separated sheet with the following columns.
    *  `gene1.id::gene2.id gene1.name::gene2.name  ffigf-score   FiN-score supporting-reads    pass-fail-code`
       *  ffigf-score of fusion A::B : Number of supporting A::B fusion reads divided by number of fusion reads mapping to gene A or gene B but not both.
       *  FiN score: Number of supporting A::B fusion reads divided by sum of number of normal A reads and B normal reads (normal being reads not supporting any gene fusions)
       *  pass-fail-code: PASS:GF for gene fusions, PASS:RT for readthroughs, FAIL::reason::... for filtered candidates.
* [output].fail: Contains called filtered fusion candidates in the same column format.

```bash
./genion run
    -i          /path/to/input/fastq           
    --gtf       /path/to/annotation/gtf        
    --gpaf      /path/to/genomic/mapping/paf    
    -s          /path/to/gene/homology/tsv      
    -d          /path/to/genomicSuperDups.txt   
    -o          /path/to/output/tsv            
```

## Files required by Genion 

GTF annotation, cDNA reference sequence and Whole genome reference sequence can be downloaded from https://ensembl.org/info/data/ftp/index.html

`genomicSuperDups.txt` can be downloaded from ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz (Should be extracted using gzip)

Homology tsv file can be produced using ENSEMBL cDNA reference and following:

```bash
minimap2 [cdna.fa] [cdna.fa] -X -t [threads] -2 -c -o [cdna.selfalign.paf]
cat [cdna.selfalign.paf] | cut -f1,6 | sed 's/_/\t/g' | awk 'BEGIN{OFS=\"\\t\";}{print substr($1,1,15),substr($2,1,15),substr($3,1,15),substr($4,1,15);}' | awk '$1!=$3' | sort | uniq > [cdna.selfalign.tsv]
```

# Genion Snakemake
We provide a snakemake file to help running genion.
* Maps Long reads
* Downloads the duplication annotation
* Prepares the Sequence similarity file
* Runs Genion

## Genion Snakemake dependencies
|Dependencies | Version |
|-------- |-----|
|c++ | gcc >= 9 or clang >= 8|
|zlib| >= 1.2.11 |
|Python   | >= 3.7 |
|[snakemake](https://snakemake.readthedocs.io/en/stable/) | >= 5.3.0 |
|[deSALT](https://github.com/ydLiu-HIT/deSALT) | >= 1.5.5 |
|[minimap2](https://github.com/lh3/minimap2/tree/master/misc) | >= 2.17 |
|[paftools](https://github.com/lh3/minimap2/tree/master/misc) |  |

Snakemake dependencies can be installed using conda/mamba
```bash
conda create --file genion.env --name genion-env
conda activate genion-env
```

## Running Genion Snakemake
After preparing a config file following the Project Configuration section, you can run snakemake with the following command.

```bash
snakemake -j [number-of-threads] --config-file [path-to-config-file]
```

## Snakemake Project Configuration
In order to run Genion, you need to create a project configuration file namely ``config.yaml``. 
This configuration consists of a number mandatory settings and some optional advance settings. 
Below is the list of the all the settings that you can set in your project.

|config-paramater-name | Type | Description|
|------------------------------|-----------|--------------------------------------------------------------------------------------------------------------------------------------|
| path                          | Mandatory | Full path to project directory.  |
| reference-dna                 | Mandatory | Full path to the DNA reference                |
| reference-cdna                | Mandatory | Full path to the cDNA reference                |
| annotation-gtf                 | Mandatory | Full path to the GTF annotation                |
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

## Snakemake Input/Output file structure

```
[path]/
├── rawdata                       
│   ├── A_clr.fastq.gz
|   └── A_ont.fastq.gz
├── analysis  (intermediate files)
│   ├── A_clr   
|   └── A_ont
└── results                       
    ├── A_clr.fusions.tsv
    ├── A_clr.readthrough.tsv
    ├── A_ont.fusions.tsv
    └── A_ont.readthrough.tsv
```

For the input/output file structure description, snakemake configuration comes with two options each for rawdata, analysis and results.
You can use `-base` suffix (like `rawdata-base`). This way snakemake will know that given path is relative to the project path.
Or you can directly use `rawdata` to enter absolute path. This may be helpful if input files are not in the project directory. 

# Simulated Dataset

## Download
Simulated gene fusion dataset can be downloaded from:
https://figshare.com/articles/dataset/Small_gene_fusion_simulated_long_read_dataset/17253821

## Contents

```
example.fastq                  #Simulated reads fastq
example.paf                    #Mapping in paf format
genomicSuperDups.txt           #Genomic segmental duplication annotation
Homo_sapiens.GRCh38.97.gtf     #Gene annotation
cdna.self.tsv                  #Homology information
```

## Setup
```bash
tar xzvf small_example.tar.gz
```

## Example genion run

```bash
cd small_example
genion -i example.fastq -d genomicSuperDups.txt --gtf Homo_sapiens.GRCh38.97.gtf -g example.paf -s cdna.self.tsv -t 1 -o output.tsv 
```






