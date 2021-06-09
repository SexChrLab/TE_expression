import os

# Environment: TE_expression

configfile: "gtex_rna_TE_processing.config.json"


rule all:
    input:
        expand(os.path.join(config["scratchdir"], "stats/{sample_name}_XX.picard_mkdup_metrics.txt"), sample_name=config["females"]),
        expand(os.path.join(config["scratchdir"], "stats/{sample_name}_XY.picard_mkdup_metrics.txt"), sample_name=config["males"]),

        expand(os.path.join(config["scratchdir"], "bams/{sample_name}_XX_sorted_mkdups_rdgrp.bam"), sample_name=config["females"]),
        expand(os.path.join(config["scratchdir"], "bams/{sample_name}_XX_sorted_mkdups_rdgrp.bam.bai"), sample_name=config["females"]),
        expand(os.path.join(config["scratchdir"], "bams/{sample_name}_XY_sorted_mkdups_rdgrp.bam"), sample_name=config["males"]),
        expand(os.path.join(config["scratchdir"], "bams/{sample_name}_XY_sorted_mkdups_rdgrp.bam.bai"), sample_name=config["males"]),

        expand(os.path.join(config["scratchdir"], "results_TEcounts_multi/{sample_name}_XX_sorted_mkdups_rdgrp_TEcount_multi.cntTable"), sample_name=config["females"]),
        expand(os.path.join(config["scratchdir"], "results_TEcounts_multi/{sample_name}_XY_sorted_mkdups_rdgrp_TEcount_multi.cntTable"), sample_name=config["males"])

#expand(os.path.join(config["scratchdir"], "bams/{sample_name}_XX_sorted_Aligned.sortedByCoord.out.bam"), sample_name=config["females"]),
#expand(os.path.join(config["scratchdir"], "bams/{sample_name}_XY_sorted_Aligned.sortedByCoord.out.bam"), sample_name=config["males"]),

#expand(os.path.join(config["scratchdir"], "bams/{sample_name}_XX_sorted_mkdups.bam"), sample_name=config["females"]),
#expand(os.path.join(config["scratchdir"], "stats/{sample_name}_XX.picard_mkdup_metrics.txt"), sample_name=config["females"]),
#expand(os.path.join(config["scratchdir"], "bams/{sample_name}_XY_sorted_mkdups.bam"), sample_name=config["males"]),
#expand(os.path.join(config["scratchdir"], "stats/{sample_name}_XY.picard_mkdup_metrics.txt"), sample_name=config["males"]),

#expand(os.path.join(config["scratchdir"], "bams/{sample_name}_XX_sorted_mkdups_rdgrp.bam"), sample_name=config["females"]),
#expand(os.path.join(config["scratchdir"], "bams/{sample_name}_XX_sorted_mkdups_rdgrp.bam.bai"), sample_name=config["females"]),
#expand(os.path.join(config["scratchdir"], "bams/{sample_name}_XY_sorted_mkdups_rdgrp.bam"), sample_name=config["males"]),
#expand(os.path.join(config["scratchdir"], "bams/{sample_name}_XY_sorted_mkdups_rdgrp.bam.bai"), sample_name=config["males"]),

#expand(os.path.join(config["scratchdir"], "results_TEcounts_multi/{sample_name}_XX_sorted_mkdups_rdgrp_TEcount_multi.cntTable"), sample_name=config["females"]),
#expand(os.path.join(config["scratchdir"], "results_TEcounts_multi/{sample_name}_XY_sorted_mkdups_rdgrp_TEcount_multi.cntTable"), sample_name=config["males"])

#-------------------------------------------------------------------------------
# Main step 1 #
# STAR: Align data
#-------------------------------------------------------------------------------
# Indexing only occurs once so I will do this outside snakefile.
# See: /scratch/amtarave/TE_expression/analysis/gtex_liver/scripts/make_STAR_index_gtex_liver.sh
# Unzip fastq files. I was getting an error with STAR before reading the fastq
# files and i think its because they shouldnt be zipped

# add --readFilesCommand zcat  to read in zipped fastqs so you dont have to unzip
# --readStrand Reverse flag is no longer avaliable in STAR
rule aling_STAR_females:
    input:
        f1 = os.path.join(config["fq_path"], "{sample_name}_1_trimmed.fq.gz"),
        f2 = os.path.join(config["fq_path"], "{sample_name}_2_trimmed.fq.gz")
    output:
        temp(os.path.join(config["scratchdir"], "bams/{sample_name}_XX_sorted_Aligned.sortedByCoord.out.bam"))
    params:
        STAR_index_XX = config["STAR_index_XX"],
        STAR_index_XY = config["STAR_index_XY"],
        gtf = config["gene_gtf"],
        stem = os.path.join(config["scratchdir"], "bams/{sample_name}_XX_sorted_")
    shell:
        """
        STAR --runMode alignReads --genomeDir {params.STAR_index_XX} --sjdbGTFfile {params.gtf} --winAnchorMultimapNmax 100 --outFilterMultimapNmax 100 --readFilesIn {input.f1} {input.f2} --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix {params.stem} --runThreadN 8
        """

rule aling_STAR_males:
    input:
        f1 = os.path.join(config["fq_path"], "{sample_name}_1_trimmed.fq.gz"),
        f2 = os.path.join(config["fq_path"], "{sample_name}_2_trimmed.fq.gz")
    output:
        temp(os.path.join(config["scratchdir"], "bams/{sample_name}_XY_sorted_Aligned.sortedByCoord.out.bam"))
    params:
        STAR_index_XX = config["STAR_index_XX"],
        STAR_index_XY = config["STAR_index_XY"],
        gtf = config["gene_gtf"],
        stem = os.path.join(config["scratchdir"], "bams/{sample_name}_XY_sorted_")
    shell:
        """
        STAR --runMode alignReads --genomeDir {params.STAR_index_XY} --sjdbGTFfile {params.gtf} --winAnchorMultimapNmax 100 --outFilterMultimapNmax 100 --readFilesIn {input.f1} {input.f2} --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix {params.stem} --runThreadN 8
        """


#-------------------------------------------------------------------------------
# Main step 2 #
# Mark duplicates
#-------------------------------------------------------------------------------
rule MarkDups_females:
    input:
        os.path.join(config["scratchdir"], "bams/{sample_name}_XX_sorted_Aligned.sortedByCoord.out.bam")
    output:
        bam = temp(os.path.join(config["scratchdir"], "bams/{sample_name}_XX_sorted_mkdups.bam")),
        metrics = os.path.join(config["scratchdir"], "stats/{sample_name}_XX.picard_mkdup_metrics.txt")
    shell:
        """
        picard -Xmx14g MarkDuplicates I={input} O={output.bam} M={output.metrics} VALIDATION_STRINGENCY=LENIENT
        """

rule MarkDups_males:
    input:
        os.path.join(config["scratchdir"], "bams/{sample_name}_XY_sorted_Aligned.sortedByCoord.out.bam")
    output:
        bam = temp(os.path.join(config["scratchdir"], "bams/{sample_name}_XY_sorted_mkdups.bam")),
        metrics = os.path.join(config["scratchdir"], "stats/{sample_name}_XY.picard_mkdup_metrics.txt")
    shell:
        """
        picard -Xmx14g MarkDuplicates I={input} O={output.bam} M={output.metrics} VALIDATION_STRINGENCY=LENIENT
        """


#-------------------------------------------------------------------------------
# Main step 3 #
# Add read groups to bams
#-------------------------------------------------------------------------------
# I am getting a java error because java version is too old for picard
# run module load java/latest before submitting job
rule add_read_grps_females:
    input:
        os.path.join(config["scratchdir"], "bams/{sample_name}_XX_sorted_mkdups.bam")
    output:
        os.path.join(config["scratchdir"], "bams/{sample_name}_XX_sorted_mkdups_rdgrp.bam")
    params:
        id = lambda wildcards: config[wildcards.sample_name]["ID"],
        sm = lambda wildcards: config[wildcards.sample_name]["SM"],
        lb = lambda wildcards: config[wildcards.sample_name]["LB"],
        pu = lambda wildcards: config[wildcards.sample_name]["PU"],
        pl = lambda wildcards: config[wildcards.sample_name]["PL"]
    shell:
        """
        picard -Xmx14g AddOrReplaceReadGroups I={input} O={output} RGID={params.id} RGPU={params.pu} RGSM={params.sm} RGPL={params.pl} RGLB={params.lb} VALIDATION_STRINGENCY=LENIENT
        """

rule index_bams_females:
    input:
        os.path.join(config["scratchdir"], "bams/{sample_name}_XX_sorted_mkdups_rdgrp.bam")
    output:
        os.path.join(config["scratchdir"], "bams/{sample_name}_XX_sorted_mkdups_rdgrp.bam.bai")
    shell:
        "samtools index {input}"

rule add_read_grps_males:
    input:
        os.path.join(config["scratchdir"], "bams/{sample_name}_XY_sorted_mkdups.bam")
    output:
        os.path.join(config["scratchdir"], "bams/{sample_name}_XY_sorted_mkdups_rdgrp.bam")
    params:
        id = lambda wildcards: config[wildcards.sample_name]["ID"],
        sm = lambda wildcards: config[wildcards.sample_name]["SM"],
        lb = lambda wildcards: config[wildcards.sample_name]["LB"],
        pu = lambda wildcards: config[wildcards.sample_name]["PU"],
        pl = lambda wildcards: config[wildcards.sample_name]["PL"]
    shell:
        """
        picard -Xmx14g AddOrReplaceReadGroups I={input} O={output} RGID={params.id} RGPU={params.pu} RGSM={params.sm} RGPL={params.pl} RGLB={params.lb} VALIDATION_STRINGENCY=LENIENT
        """

rule index_bams_males:
    input:
        os.path.join(config["scratchdir"], "bams/{sample_name}_XY_sorted_mkdups_rdgrp.bam")
    output:
        os.path.join(config["scratchdir"], "bams/{sample_name}_XY_sorted_mkdups_rdgrp.bam.bai")
    shell:
        "samtools index {input}"


#-------------------------------------------------------------------------------
# Main step 4 #
# Generate TE count table in uniq mode (this will cound all reads - both unique
# and multi mapped reads)
#-------------------------------------------------------------------------------
# TEtranscripts does the DE analysis for you. See command below. I will run
# TEcount which justs outputs the count file
rule TE_counts_multi_females:
    input:
        os.path.join(config["scratchdir"], "bams/{sample_name}_XX_sorted_mkdups_rdgrp.bam")
    params:
        tegtf = config["TE_gtf"],
        genegtf = config["gene_gtf"],
        outdir = os.path.join(config["scratchdir"], "results_TEcounts_multi/{sample_name}_XX_sorted_mkdups_rdgrp_TEcount_multi")
    output:
        os.path.join(config["scratchdir"], "results_TEcounts_multi/{sample_name}_XX_sorted_mkdups_rdgrp_TEcount_multi.cntTable")
    shell:
        """
        TEcount --sortByPos --format BAM --GTF {params.genegtf} --TE {params.tegtf} --mode multi -b {input} --stranded no --project {params.outdir}
        """

rule TE_counts_multi_males:
    input:
        os.path.join(config["scratchdir"], "bams/{sample_name}_XY_sorted_mkdups_rdgrp.bam")
    params:
        tegtf = config["TE_gtf"],
        genegtf = config["gene_gtf"],
        outdir = os.path.join(config["scratchdir"], "results_TEcounts_multi/{sample_name}_XY_sorted_mkdups_rdgrp_TEcount_multi")
    output:
        os.path.join(config["scratchdir"], "results_TEcounts_multi/{sample_name}_XY_sorted_mkdups_rdgrp_TEcount_multi.cntTable")
    shell:
        """
        TEcount --sortByPos --format BAM --GTF {params.genegtf} --TE {params.tegtf} --mode multi -b {input} --stranded no --project {params.outdir}
        """
