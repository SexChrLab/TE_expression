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
The pipeline using this data can be found `analysis/tcag_lihc_tumor/scripts`.

```
# Make directory for scripts and results
cd /scratch/amtarave/TE_expression/analysis
mkdir tcag_lihc_tumor
```

Path to RNA-seq bam files: `/data/CEM/shared/controlled_access/TCGA/LIHC/RNA/RNAseq/bam/`

Info on samples to process (primary tumor RNA seq bams): `/scratch/amtarave/introgression_pilot/elife/TCGA_exome_processing/bam.primary.tumor.prefix.txt`.


*TODO: For TCGA (and GTEx), I need to check if the samples are unstranded, forward or reverse. Can figure this out using Salmon. See: https://salmon.readthedocs.io/en/latest/salmon.html#what-s-this-libtype*

https://littlebitofdata.com/en/2017/08/strandness_in_rnaseq/

Install salmon to environment:
```
conda install -c bioconda salmon

# Test command on one bam file
salmon quant -t transcripts.fa -l A -a aln.bam -o salmon_quant

salmon quant -t /data/CEM/shared/public_data/references/GENCODE/gencode.v29.annotation.gtf -l A -a /data/CEM/shared/controlled_access/TCGA/LIHC/RNA/RNAseq/bam/2V-A95S-01A-XY_non_maskRNAs_sorted.bam -o salmon_quant_test

salmon quant -t /data/CEM/shared/public_data/references/GENCODE/gencode.v29.annotation.gtf -l A -a /data/CEM/shared/controlled_access/TCGA/LIHC/RNA/RNAseq/bam/2V-A95S-01A-XY_non_maskRNAs_sorted.bam -o salmon_quant_test

# Step 1: Build index for transcriptome
salmon index -t athal.fa.gz -i athal_index

# Step 2: Run salmon quant to get library type
salmon quant -i <transcriptome index> --libType A \
             -o <out dir and prefix> -1 test_1.fq.gz -2 test_2.fq.gz
             -p <assigned threads>

cd /scratch/amtarave/TE_expression/analysis/tcag_lihc_tumor
mkdir library_type
cd library_type/
cp /data/CEM/shared/public_data/references/GENCODE/gencode.v29.annotation.gtf .
# I think I need the fa not gtf
cd /data/CEM/shared/public_data/references/GENCODE
cp gencode.v29.transcripts.fa /scratch/amtarave/TE_expression/analysis/tcag_lihc_tumor/library_type/

#salmon index -t gencode.v29.annotation.gtf -i gencode.v29.annotation_index
salmon index -t gencode.v29.transcripts.fa -i gencode.v29.annotation_fa_index

#salmon quant -t /scratch/amtarave/TE_expression/analysis/tcag_lihc_tumor/library_type/gencode.v29.annotation_index -l A -a /data/CEM/shared/controlled_access/TCGA/LIHC/RNA/RNAseq/bam/2V-A95S-01A-XY_non_maskRNAs_sorted.bam -o 2V-A95S-01A-XY_salmon_quant_test

salmon quant -i /scratch/amtarave/TE_expression/analysis/tcag_lihc_tumor/library_type/gencode.v29.annotation_index --libType A -o 2V-A95S-01A-XY_salmon_quant_test -1 /scratch/amtarave/TE_expression/analysis/tcag_lihc_tumor/fastqs_converted/2V-A95S-01A-XY_1.fq.gz -2 /scratch/amtarave/TE_expression/analysis/tcag_lihc_tumor/fastqs_converted/2V-A95S-01A-XY_2.fq.gz

# From Salmon Log
[2021-04-29 09:20:12.295] [jointLog] [info] Automatically detected most likely library type as IU
# so it is unstranded
```

Notes:
For gtex and tcga need to start with bams. See: https://gatk.broadinstitute.org/hc/en-us/articles/360036485372-SamToFastq-Picard-. So go from bams to fastqs to bams https://www.biostars.org/p/326714/

bam to fastq: https://github.com/broadinstitute/picard/issues/715
add `VALIDATION_STRINGENCY=SILENT` to command

```
# Example commands
java -jar picard.jar SamToFastq
     I=input.bam
     FASTQ=output.fastq
java -jar picard.jar SamToFastq I=<file_alnMAP.bam> FASTQ=<filemap_1.fq> SECOND_END_FASTQ=<filemap_2.fq>

#or

samtools fastq -1 <filemap_1.fq.gz> -2 <filemap_1.fq.gz> <file_alnMAP.bam>
picard SamToFastq \
  -Xmx2g \ #THIS
  I=input.bam \
  F=out_unmapped_R1.fastq \
  F2=out_unmapped_R2.fastq \
  VALIDATION_STRINGENCY=SILENT


```

Reads were sequenced to 49 bps, so need to make new STAR indexes for this data
```
# Make star index in separate script since this only needs to happen once
cd /scratch/amtarave/TE_expression/analysis/tcag_lihc_tumor/scripts
touch make_STAR_index_lihc.sh

cd /scratch/amtarave/TE_expression
mkdir STAR_index_GRCh38_p12_genome_XXonly_lihc
mkdir STAR_index_GRCh38_p12_genome_XY_lihc

```

Submit pipeline
```
sbatch lihc_rna_TE_processing.sh
Submitted batch job 9383516

# memory failure for some star runs. Edit sh to allocate more mem per cpu
sbatch lihc_rna_TE_processing.sh
Submitted batch job 9384529

```

## Creating the pipeline - GTEx normal liver
The pipeline using this data can be found `analysis/gtex_liver/scripts`.

Converted bams and multiqc reports can be found: `/scratch/amtarave/gtex/`.

Reads seem to be sequenced to 75 bps, so will have to make another STAR index.

```
# Make star index in separate script since this only needs to happen once
cd /scratch/amtarave/TE_expression/analysis/gtex_liver/scripts
touch make_STAR_index_gtex_liver.sh

cd /scratch/amtarave/TE_expression
mkdir STAR_index_GRCh38_p12_genome_XXonly_gtex_liver
mkdir STAR_index_GRCh38_p12_genome_XY_gtex_liver

```

Sample information. Downloaded sample information from: https://www.gtexportal.org/home/datasets. Go to GTEx Analysis V8 (dbGaP Accession phs000424.v8.p2), then under "Annotations" there are 4 files:
- GTEx_Analysis_v8_Annotations_SampleAttributesDD.xlsx
- GTEx_Analysis_v8_Annotations_SubjectPhenotypesDD.xlsx
- GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt
- GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt

`GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt` has sex information. 1 = Male, 2 = Female

```
mkdir sample_info
cd sample_info
wget https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDD.xlsx
wget https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDD.xlsx
wget https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt
wget https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt

head GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt
SUBJID  SEX     AGE     DTHHRDY
GTEX-1117F      2       60-69   4
GTEX-111CU      1       50-59   0
GTEX-111FC      1       60-69   1
GTEX-111VG      1       60-69   3
GTEX-111YS      1       60-69   0
GTEX-1122O      2       60-69   0
GTEX-1128S      2       60-69   2
GTEX-113IC      1       60-69   
GTEX-113JC      2       50-59   2

# note: SUBJID is different from the ID for the RNA-seq data
# Make list of liver samples and get SUBJID from that. Take liver samples from config
vi liver_samples.txt
head liver_samples.txt
GTEX-11DXY-0526-SM-5EGGQ
GTEX-11DXZ-0126-SM-5EGGY
GTEX-11EQ9-0526-SM-5A5JZ
GTEX-11GSP-0626-SM-5986T
GTEX-11NUK-1226-SM-5P9GM
GTEX-11NV4-1326-SM-5HL6V
GTEX-11OF3-0726-SM-5BC4Z
GTEX-11TT1-1726-SM-5EQLJ
GTEX-11TUW-1726-SM-5BC5C
GTEX-11WQC-0726-SM-5EQMR

# to get SUBJID just cut the first two entries sep by '-'
cut -f1,2 -d'-' liver_samples.txt > liver_SUBJID.txt

grep -f liver_SUBJID.txt GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt | grep "        1       " | cut -f1 > liver_SUBJID_males.txt
grep -f liver_SUBJID.txt GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt | grep "        2       " | cut -f1 > liver_SUBJID_females.txt

grep -f liver_SUBJID_males.txt liver_samples.txt > liver_samples_males.txt
grep -f liver_SUBJID_females.txt liver_samples.txt > liver_samples_females.txt

```


Get library type for this data. Test on one sample.

```
cd /scratch/amtarave/TE_expression/analysis/gtex_liver
mkdir library_type
cd library_type/
# use index created when I ran this for tcga

salmon quant -i /scratch/amtarave/TE_expression/analysis/tcag_lihc_tumor/library_type/gencode.v29.annotation_fa_index --libType A -o GTEX-11DXY-0526-SM-5EGGQ_salmon_quant_test -1 /scratch/amtarave/gtex/trimmed_fastqs/GTEX-11DXY-0526-SM-5EGGQ_1_trimmed.fq.gz -2 /scratch/amtarave/gtex/trimmed_fastqs/GTEX-11DXY-0526-SM-5EGGQ_2_trimmed.fq.gz

# From output
[2021-06-09 10:50:57.253] [jointLog] [info] Automatically detected most likely library type as IU
# IU is unstranded
```

Run pipeline.

```
sbatch gtex_rna_TE_processing.sh
Submitted batch job 9682191
```

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
