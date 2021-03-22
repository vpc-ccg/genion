# Genion
## Installation


## Installation from Source
```bash
git clone git@bitbucket.org:CharlesTaylor/genion.git
cd genion
make
```

## Prerequisites


Dependencies | Version 
-------- |-----|
Python   | 3.x |
[deSALT](https://github.com/ydLiu-HIT/deSALT) | >= 1.5.5 |
[minimap2](https://github.com/lh3/minimap2/tree/master/misc) | >= 2.17 |
[paftools](https://github.com/lh3/minimap2/tree/master/misc) |  |
[snakemake](https://snakemake.readthedocs.io/en/stable/) | >= 5.3.0 |
[lightgbm*](https://lightgbm.readthedocs.io/en/latest/) | >= 2.3.2 |

lightgbm is required for optional chimeric read correction step.

## Running Genion

Genion is runned using run.sh

```bash
run.sh --configfile /path/to/config.yaml -j [num-threads]
```
./run.sh is a snakemake container script. Any snakemake command can be used with it.

## Reference building
Reference for Genion can be build using: 
```bash
./ref-build.sh --configfile=config.yaml [ref-name]/1.gtf
```
and passing following in a config file
```yaml
dna_ref:
    Homo_sapiens.GRCh38.dna.primary_assembly.fa
cdna_ref:
    Homo_sapiens.GRCh38.cdna.all.fa
gtf:
    Homo_sapiens.GRCh38.97.gtf
```
All of these files can be downloaded from https://uswest.ensembl.org/info/data/ftp/index.html.

## Project Configuration
In order to run Genion, you need to create a project configuration file namely ``config.yaml``. 
This configuration consists of a number mandatory settings and some optional advance settings. 
Below is the list of the all the settings that you can set in your project.

|config-paramater-name | Type | Description|
|------------------------------|-----------|--------------------------------------------------------------------------------------------------------------------------------------|
| path                         | Mandatory | Full path to project directory.  |
| rawdata-base                 | Mandatory | Location of the input fastq files relative to ``path``.                                                         |
| reference                    | Mandatory | Full path to the reference build by ./ref-build.sh.                                                               |
| input                        | Mandatory | A list of input files per sample.       |
| duplications                 | Mandatory | Path to genomicSuperDups.txt file. Can be downloaded from ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz (Should be extracted using gzip) |
| analysis-base                | Optional  | Location of intermediate files relative to ``path``. default: ``{path}/analysis``|
| results-base                 | Optional  | Location of final results relative to the ``path``. default: ``{path}/results``  |
| wg-aligner                   | Optional  | Mapper to use (``deSALT``, ``minimap2``) default: ``deSALT``                                                                      |
| ext                          | Optinal   | extension of the fastq files used in input (``fastq``,``fastq.gz``,``fq``,``fq.gz``) default:``fastq`` |

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
reference:
    /path/to/reference/
duplications:
    /path/to/genomicSuperDups.txt
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
