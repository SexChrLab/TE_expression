import os

# Environment: TE_expression

configfile: "lihc_rna_TE_processing.config.json"


rule all:
    input:
        expand(os.path.join(config["scratchdir"], "fastqs_converted/{sample_name}_1.fq.gz"), sample_name=config["sample_name"]),
        expand(os.path.join(config["scratchdir"], "fastqs_converted/{sample_name}_2.fq.gz"), sample_name=config["sample_name"]),
        expand(os.path.join(config["scratchdir"], "fastqc/{sample_name}_1_fastqc.html"), sample_name=config["sample_name"]),
        expand(os.path.join(config["scratchdir"], "fastqc/{sample_name}_2_fastqc.html"), sample_name=config["sample_name"]),
        os.path.join(config["scratchdir"], "multiqc/multiqc_report.html"),

        expand(os.path.join(config["scratchdir"], "bams/{sample_name}_sorted_Aligned.sortedByCoord.out.bam"), sample_name=config["sample_name"]),

        expand(os.path.join(config["scratchdir"], "bams/{sample_name}_sorted_mkdups.bam"), sample_name=config["sample_name"]),
        expand(os.path.join(config["scratchdir"], "stats/{sample_name}.picard_mkdup_metrics.txt"), sample_name=config["sample_name"]),

        expand(os.path.join(config["scratchdir"], "bams/{sample_name}_sorted_mkdups_rdgrp.bam"), sample_name=config["sample_name"]),
        expand(os.path.join(config["scratchdir"], "bams/{sample_name}_sorted_mkdups_rdgrp.bam.bai"), sample_name=config["sample_name"]),

        expand(os.path.join(config["scratchdir"], "results_TEcounts_multi/{sample_name}_sorted_mkdups_rdgrp_TEcount_multi.cntTable"), sample_name=config["sample_name"])


#-------------------------------------------------------------------------------
# Main step 1 #
# SamToFastq: Convert bam files to fastq
#-------------------------------------------------------------------------------
rule bamToFq:
    input:
        os.path.join(config["bam_path"], "{sample_name}_non_maskRNAs_sorted.bam")
    output:
        fq1 = os.path.join(config["scratchdir"], "fastqs_converted/{sample_name}_1.fq.gz"),
        fq2 = os.path.join(config["scratchdir"], "fastqs_converted/{sample_name}_2.fq.gz")
    shell:
        """
        picard -Xmx14g SamToFastq I={input} FASTQ={output.fq1} SECOND_END_FASTQ={output.fq2} VALIDATION_STRINGENCY=SILENT
        """

#-------------------------------------------------------------------------------
# Main step 2 #
# Initial QC vizualization
#-------------------------------------------------------------------------------
rule fastqc_analysis:
    input:
        f1 = os.path.join(config["scratchdir"], "fastqs_converted/{sample_name}_1.fq.gz"),
        f2 = os.path.join(config["scratchdir"], "fastqs_converted/{sample_name}_2.fq.gz")
    output:
        f1 = os.path.join(config["scratchdir"], "fastqc/{sample_name}_1_fastqc.html"),
        f2 = os.path.join(config["scratchdir"], "fastqc/{sample_name}_2_fastqc.html")
    params:
        fastqc = os.path.join(config["scratchdir"], "fastqc/")
    shell:
        """
        PERL5LIB="";
        fastqc -o {params.fastqc} {input.f1};
        fastqc -o {params.fastqc} {input.f2}
        """

rule multiqc_analysis:
    input:
        f1 = expand(os.path.join(config["scratchdir"], "fastqc/{sample_name}_1_fastqc.html"),
        sample_name=config["sample_name"]),
        f2 = expand(os.path.join(config["scratchdir"], "fastqc/{sample_name}_2_fastqc.html"),
        sample_name=config["sample_name"])
    output:
        os.path.join(config["scratchdir"], "multiqc/multiqc_report.html")
    params:
        fastqc = os.path.join(config["scratchdir"], "fastqc/"),
        multiqc = os.path.join(config["scratchdir"], "multiqc/")
    shell:
        "export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 &&"
        "multiqc --interactive -o {params.multiqc} {params.fastqc}"

# All looks good so move on to alignment. One sample may need to be removed though
# because of low sequence quality (ES-A2HS-01A-XY_2)

#-------------------------------------------------------------------------------
# Main step 3 #
# STAR: Align data
#-------------------------------------------------------------------------------
# Indexing only occurs once so I will do this outside snakefile.
# See: /scratch/amtarave/TE_expression/analysis_test/tcag_lihc_tumor/scripts/make_STAR_index_lihc.sh
# Unzip fastq files. I was getting an error with STAR before reading the fastq
# files and i think its because they shouldnt be zipped

# add --readFilesCommand zcat  to read in zipped fastqs so you dont have to unzip
# --readStrand Reverse flag is no longer avaliable in STAR
rule aling_STAR:
    input:
        f1 = os.path.join(config["scratchdir"], "fastqs_converted/{sample_name}_1.fq.gz"),
        f2 = os.path.join(config["scratchdir"], "fastqs_converted/{sample_name}_2.fq.gz")
    output:
        os.path.join(config["scratchdir"], "bams/{sample_name}_sorted_Aligned.sortedByCoord.out.bam")
    params:
        STAR_index_XX = config["STAR_index_XX"],
        STAR_index_XY = config["STAR_index_XY"],
        gtf = config["gene_gtf"],
        stem = os.path.join(config["scratchdir"], "bams/{sample_name}_sorted_")
    run:
        if "XY" in wildcards.sample_name:
            shell("STAR --runMode alignReads --genomeDir {params.STAR_index_XY} --sjdbGTFfile {params.gtf} --winAnchorMultimapNmax 100 --outFilterMultimapNmax 100 --readFilesIn {input.f1} {input.f2} --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix {params.stem} --runThreadN 8")
        if "XX" in wildcards.sample_name:
            shell("STAR --runMode alignReads --genomeDir {params.STAR_index_XX} --sjdbGTFfile {params.gtf} --winAnchorMultimapNmax 100 --outFilterMultimapNmax 100 --readFilesIn {input.f1} {input.f2} --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix {params.stem} --runThreadN 8")


#-------------------------------------------------------------------------------
# Main step 4 #
# Mark duplicates
#-------------------------------------------------------------------------------
rule MarkDups:
    input:
        os.path.join(config["scratchdir"], "bams/{sample_name}_sorted_Aligned.sortedByCoord.out.bam")
    output:
        bam = os.path.join(config["scratchdir"], "bams/{sample_name}_sorted_mkdups.bam"),
        metrics = os.path.join(config["scratchdir"], "stats/{sample_name}.picard_mkdup_metrics.txt")
    shell:
        """
        picard -Xmx14g MarkDuplicates I={input} O={output.bam} M={output.metrics} VALIDATION_STRINGENCY=LENIENT
        """


#-------------------------------------------------------------------------------
# Main step 5 #
# Add read groups to bams
#-------------------------------------------------------------------------------
# I am getting a java error because java version is too old for picard
# run module load java/latest before submitting job
rule add_read_grps:
    input:
        os.path.join(config["scratchdir"], "bams/{sample_name}_sorted_mkdups.bam")
    output:
        os.path.join(config["scratchdir"], "bams/{sample_name}_sorted_mkdups_rdgrp.bam")
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

rule index_bams:
    input:
        os.path.join(config["scratchdir"], "bams/{sample_name}_sorted_mkdups_rdgrp.bam")
    output:
        os.path.join(config["scratchdir"], "bams/{sample_name}_sorted_mkdups_rdgrp.bam.bai")
    shell:
        "samtools index {input}"

#-------------------------------------------------------------------------------
# Main step 6 #
# Generate TE count table in uniq mode (this will cound all reads - both unique
# and multi mapped reads)
#-------------------------------------------------------------------------------
# TEtranscripts does the DE analysis for you. See command below. I will run
# TEcount which justs outputs the count file
# TO DO: Determine library type using salmon
rule TE_counts_multi:
    input:
        os.path.join(config["scratchdir"], "bams/{sample_name}_sorted_mkdups_rdgrp.bam")
    params:
        tegtf = config["TE_gtf"],
        genegtf = config["gene_gtf"],
        outdir = os.path.join(config["scratchdir"], "results_TEcounts_multi/{sample_name}_sorted_mkdups_rdgrp_TEcount_multi")
    output:
        os.path.join(config["scratchdir"], "results_TEcounts_multi/{sample_name}_sorted_mkdups_rdgrp_TEcount_multi.cntTable")
    shell:
        """
        TEcount --sortByPos --format BAM --GTF {params.genegtf} --TE {params.tegtf} --mode multi -b {input} --stranded no --project {params.outdir}
        """
