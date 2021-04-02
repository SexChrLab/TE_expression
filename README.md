# TE_expression


## Background
This README will serve as documentation of the steps in the pipeline to quantify transposable element (TE) expression in different tissues. After I process the data, differential TE expression analysis will be performed between males and females for 1st trimester placenta, full-term placenta, healthy liver tissue, and tumor tissue from liver cancer patients. I will first quantify TE expression at the family level. Eventually I will quantify TE expression at the locus level.


## Create conda environment and install necessary software
```
# Create a conda environment called TE_expression
conda create --name TE_expression

# To activate this environment, use
#
#     $ conda activate TE_expression
#
# To deactivate an active environment, use
#
#     $ conda deactivate

# Add software to environment
# First, activate the environment
conda activate TE_expression

conda install -c bioconda fastqc
conda install -c bioconda multiqc
conda install -c agbiome bbtools
conda install -c bioconda star
conda install -c bioconda snakemake

conda install -c bioconda tetranscripts
```


## Creating the pipeline - Full-term placentas
Since the data need to be aligned with parameters optimized for multi mapped reads, I will start from FASTQ files. Quality of FASTQ reads will be vizualized using fastqc and multiqc, trimming will be performed using bbduk, and quality of trimmed reads will be vizualized using fastqc and multiqc. Next, I will align the trimmed FASTQ files using STAR. TE quantification will be performed using TEtranscripts. With the TE count data, differential expression analysis will be performed using the limma/voom workflow.

**NOTE**: Scripts for the lab generated full-term placenta data can be found `analysis/placenta/scripts`.

### Step 1. Create snakefile and associated files
```
cd /scratch/amtarave/TE_expression/analysis_test/placenta/scripts
touch placenta_rna_TE_processing.snakefile
touch placenta_rna_TE_processing.sh
touch placenta_rna_TE_processing.config.json

```
### Step 2: Generate STAR indexes
This should only need to be done one time. However, `--sjdbOverhang ` option may need to be changed depending on the read length of the sequence data. For example, the placenta data was sequences in 100 bp reads. This may differ for the other data sets.


Example command:
```
STAR --runThreadN 8 \
--runMode genomeGenerate \
--genomeDir STAR_index \ # add full path to STAR_index_XX or STAR_index_XY (/scratch/amtarave/TE_expression/STAR_index_GRCh38_p12_genome_XXonly)
--genomeFastaFiles /data/CEM/shared/public_data/references/GENCODE/GRCh38.p12.genome.fa \ # or sex chromosome complement reference
#--sjdbGTFfile data/CEM/shared/public_data/references/GENCODE/gencode.v29.annotation.gtf \ # this is optional, so for first pass I am not going to use this. I will just use this in the quantification step
--sjdbOverhang 99 # to determine after fastqc (read length - 1). Read length is 100 so set to 99
```

```
# Make star index in separate script since this only needs to happen once
touch make_STAR_index.sh

cd /scratch/amtarave/TE_expression
mkdir STAR_index_GRCh38_p12_genome_XXonly
mkdir STAR_index_GRCh38_p12_genome_XY

```

### Step 3: Run script on Agave
```
sbatch placenta_rna_TE_processing.sh
Submitted batch job 8679518
```


### Notes
#### Placenta samples that failed
For this pipeline, I did not process the following samples. These are samples that were noted to have failed in a previous analysis in sex differences in gene expression in the placenta.

```
OBG0174-1-020
OBG0175-1-020
OBG0015-AMP-RNA1
OBG0015-AMP-RNA2
OBG0065-BMP-RNA1
OBG0065-BMP-RNA2
OBG0188-BFP-RNA1
OBG0188-BFP-RNA2
OBG0014_WMP-RNA1
OBG0014_WMP-RNA2
OBG0026-WFP-RNA1
OBG0026-WFP-RNA2
YPOPS0007M-AMP-RNA1
YPOPS0007M-AMP-RNA2
OBG00021_BMP # could not find in config
OBG0019-WMP-RNA1
OBG0019-WMP-RNA2

```


## Creating the pipeline - TCGA liver cancer tumors
TBD

For gtex and tcga need to start with bams. See: https://gatk.broadinstitute.org/hc/en-us/articles/360036485372-SamToFastq-Picard-. So go from bams to fastqs to bams https://www.biostars.org/p/326714/

The pipeline using this data can be found `analysis/tcag_lihc_tumor/scripts`.

```
# Make directory for scripts and results
cd /scratch/amtarave/TE_expression/analysis
mkdir tcag_lihc_tumor
```

*TODO: For TCGA (and GTEx), I need to check if the samples are unstranded, forward or reverse. Can figure this out using Salmon. See: https://salmon.readthedocs.io/en/latest/salmon.html#what-s-this-libtype*

Path to RNA-seq bam files:
`/data/CEM/shared/controlled_access/TCGA/LIHC/RNA/RNAseq/bam`

These bams will need to be converted back to FASTQs and then re-aligned with STAR. All other steps should be the same.


## Creating the pipeline - GTEx normal liver
TBD


## Creating the pipeline - 1st term placentas
TBD


## Resources
### Tools used:

Software | Website
--- | ---
FastQC | https://www.bioinformatics.babraham.ac.uk/projects/download.html
MultiQC | https://multiqc.info/
bbduk | https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/
STAR | https://github.com/alexdobin/STAR
TEtranscripts | https://github.com/mhammell-laboratory/TEtranscripts
limma | http://www.bioconductor.org/packages/release/bioc/html/limma.html
edgeR | https://bioconductor.org/packages/release/bioc/html/edgeR.html
