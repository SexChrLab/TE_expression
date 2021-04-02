import os

# Environment: TE_expression

configfile: "placenta_rna_TE_processing.config.json"

fastq_prefixes = [config[x]["fq1_sy"][:-12] for x in config["sample_name"]] + [config[x]["fq2_sy"][:-12] for x in config["sample_name"]]

rule all:
    input:
        expand(os.path.join(config["scratchdir"], "fastq_files/{sample_name}_R1.fastq.gz"), sample_name=config["sample_name"]),
        expand(os.path.join(config["scratchdir"], "fastq_files/{sample_name}_R2.fastq.gz"), sample_name=config["sample_name"]),
#        expand(os.path.join(config["scratchdir"], "fastqc/{sample_name}_R1_fastqc.html"), sample_name=config["sample_name"]),
#        expand(os.path.join(config["scratchdir"], "fastqc/{sample_name}_R2_fastqc.html"), sample_name=config["sample_name"]),
#        os.path.join(config["scratchdir"], "multiqc/multiqc_report.html"),
#        #expand(os.path.join(config["scratchdir"], "trimmed_fastqs/{sample_name}_R1_trimmed.fastq.gz"), sample_name=config["sample_name"]),#
#        #expand(os.path.join(config["scratchdir"], "trimmed_fastqs/{sample_name}_R2_trimmed.fastq.gz"), sample_name=config["sample_name"]),
#        #expand(os.path.join(config["scratchdir"], "fastqc_trimmed/{sample_name}_R1_trimmed_fastqc.html"), sample_name=config["sample_name"]),
#        #expand(os.path.join(config["scratchdir"], "fastqc_trimmed/{sample_name}_R2_trimmed_fastqc.html"), sample_name=config["sample_name"]),
#        #os.path.join(config["scratchdir"], "multiqc_trimmed/multiqc_report.html"),
#        #expand(os.path.join(config["scratchdir"], "trimmed_fastqs/{sample_name}_R1_trimmed.fastq"), sample_name=config["sample_name"]),
#        #expand(os.path.join(config["scratchdir"], "trimmed_fastqs/{sample_name}_R2_trimmed.fastq"), sample_name=config["sample_name"]),
        expand(os.path.join(config["scratchdir"], "bams/{sample_name}_XX_sorted_Aligned.sortedByCoord.out.bam"), sample_name=config["females"]),
        expand(os.path.join(config["scratchdir"], "bams/{sample_name}_XY_sorted_Aligned.sortedByCoord.out.bam"), sample_name=config["males"]),
        expand(os.path.join(config["scratchdir"], "bams/{sample_name}_XX_sorted_mkdups.bam"), sample_name=config["females"]),
        expand(os.path.join(config["scratchdir"], "stats/{sample_name}.XX.picard_mkdup_metrics.txt"), sample_name=config["females"]),
        expand(os.path.join(config["scratchdir"], "bams/{sample_name}_XY_sorted_mkdups.bam"), sample_name=config["males"]),
        expand(os.path.join(config["scratchdir"], "stats/{sample_name}.XY.picard_mkdup_metrics.txt"), sample_name=config["males"]),
        expand(os.path.join(config["scratchdir"], "bams/{sample_name}_XX_sorted_mkdups_rdgrp.bam"), sample_name=config["females"]),
        expand(os.path.join(config["scratchdir"], "bams/{sample_name}_XY_sorted_mkdups_rdgrp.bam"), sample_name=config["males"]),
        expand(os.path.join(config["scratchdir"], "bams/{sample_name}_XX_sorted_mkdups_rdgrp.bam.bai"), sample_name=config["females"]),
        expand(os.path.join(config["scratchdir"], "bams/{sample_name}_XY_sorted_mkdups_rdgrp.bam.bai"), sample_name=config["males"]),
        expand(os.path.join(config["scratchdir"], "results_TEcounts_multi/{sample_name}_XX_sorted_mkdups_rdgrp_TEcount_multi.cntTable"), sample_name=config["females"]),
#        #expand(os.path.join(config["scratchdir"], "results_TEcounts_uniq/{sample_name}_XX_sorted_mkdups_rdgrp_TEcount_uniq.cntTable"), sample_name=config["females"]),
        expand(os.path.join(config["scratchdir"], "results_TEcounts_multi/{sample_name}_XY_sorted_mkdups_rdgrp_TEcount_multi.cntTable"), sample_name=config["males"]),
#        #expand(os.path.join(config["scratchdir"], "results_TEcounts_uniq/{sample_name}_XY_sorted_mkdups_rdgrp_TEcount_uniq.cntTable"), sample_name=config["males"])




#-------------------------------------------------------------------------------
# Main step 1 #
# Prep files: Make make symbolic link for fastqs
#-------------------------------------------------------------------------------
rule make_symbolic_link_for_fastqs:
    input:
        f1 = lambda wildcards: config[wildcards.sample_name]["fq_path"] + config[wildcards.sample_name]["fq1"],
        f2 = lambda wildcards: config[wildcards.sample_name]["fq_path"] + config[wildcards.sample_name]["fq2"]
    output:
        f1 = os.path.join(config["scratchdir"], "fastq_files/{sample_name}_R1.fastq.gz"),
        f2 = os.path.join(config["scratchdir"], "fastq_files/{sample_name}_R2.fastq.gz")
    shell:
        """
        ln -s {input.f1} {output.f1} && touch -h {output.f1};
        ln -s {input.f2} {output.f2} && touch -h {output.f2}
        """

#-------------------------------------------------------------------------------
# Main step 2 #
# Initial QC vizualization
#-------------------------------------------------------------------------------
rule fastqc_analysis:
    input:
        f1 = os.path.join(config["scratchdir"], "fastq_files/{sample_name}_R1.fastq.gz"),
        f2 = os.path.join(config["scratchdir"], "fastq_files/{sample_name}_R2.fastq.gz")
    output:
        f1 = os.path.join(config["scratchdir"], "fastqc/{sample_name}_R1_fastqc.html"),
        f2 = os.path.join(config["scratchdir"], "fastqc/{sample_name}_R2_fastqc.html")
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
        f1 = expand(os.path.join(config["scratchdir"], "fastqc/{sample_name}_R1_fastqc.html"),
        sample_name=config["sample_name"]),
        f2 = expand(os.path.join(config["scratchdir"], "fastqc/{sample_name}_R2_fastqc.html"),
        sample_name=config["sample_name"])
    output:
        os.path.join(config["scratchdir"], "multiqc/multiqc_report.html")
    params:
        fastqc = os.path.join(config["scratchdir"], "fastqc/"),
        multiqc = os.path.join(config["scratchdir"], "multiqc/")
    shell:
        "export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 &&"
        "multiqc --interactive -o {params.multiqc} {params.fastqc}"


#-------------------------------------------------------------------------------
# Main step 3 #
# Trim FASTQ files
#-------------------------------------------------------------------------------
rule trim_bbduk:
    input:
        f1 = os.path.join(config["scratchdir"], "fastq_files/{sample_name}_R1.fastq.gz"),
        f2 = os.path.join(config["scratchdir"], "fastq_files/{sample_name}_R2.fastq.gz")
    output:
        f1 = os.path.join(config["scratchdir"], "trimmed_fastqs/{sample_name}_R1_trimmed.fastq.gz"),
        f2 = os.path.join(config["scratchdir"], "trimmed_fastqs/{sample_name}_R2_trimmed.fastq.gz")
    params:
        adapter = config["adapters"]
    shell:
        "bbduk.sh -Xmx3g in1={input.f1} in2={input.f2} "
        "out1={output.f1} out2={output.f2} "
        "ref={params.adapter} "
        "qtrim=rl trimq=30 minlen=75 maq=20"


#-------------------------------------------------------------------------------
# Main step 4 #
# Post-QC vizualization
#-------------------------------------------------------------------------------
rule fastqc_analysis_trimmed:
    input:
        f1 = os.path.join(config["scratchdir"], "trimmed_fastqs/{sample_name}_R1_trimmed.fastq.gz"),
        f2 = os.path.join(config["scratchdir"], "trimmed_fastqs/{sample_name}_R2_trimmed.fastq.gz")
    output:
        f1 = os.path.join(config["scratchdir"], "fastqc_trimmed/{sample_name}_R1_trimmed_fastqc.html"),
        f2 = os.path.join(config["scratchdir"], "fastqc_trimmed/{sample_name}_R2_trimmed_fastqc.html")
    params:
        fastqc = os.path.join(config["scratchdir"], "fastqc_trimmed/")
    shell:
        """
        PERL5LIB="";
        fastqc -o {params.fastqc} {input.f1};
        fastqc -o {params.fastqc} {input.f2}
        """


rule multiqc_analysis_trimmed:
    input:
        f1 = expand(os.path.join(config["scratchdir"], "fastqc_trimmed/{sample_name}_R1_trimmed_fastqc.html"),
        sample_name=config["sample_name"]),
        f2 = expand(os.path.join(config["scratchdir"], "fastqc_trimmed/{sample_name}_R2_trimmed_fastqc.html"),
        sample_name=config["sample_name"])
    output:
        os.path.join(config["scratchdir"], "multiqc_trimmed/multiqc_report.html")
    params:
        fastqc = os.path.join(config["scratchdir"], "fastqc_trimmed/"),
        multiqc = os.path.join(config["scratchdir"], "multiqc_trimmed/")
    shell:
        "export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 &&"
        "multiqc --interactive -o {params.multiqc} {params.fastqc}"


#-------------------------------------------------------------------------------
# Main step 5 #
# STAR: Align data
#-------------------------------------------------------------------------------
# Indexing only occurs once so I will do this outside snakefile.
# See: /scratch/amtarave/TE_expression/analysis_test/placenta/scripts/make_STAR_index.sh
# Unzip fastq files. I was getting an error with STAR before reading the fastq
# files and i think its because they shouldnt be zipped
'''
rule unzip_fastq:
    input:
        f1 = os.path.join(config["scratchdir"], "trimmed_fastqs/{sample_name}_R1_trimmed.fastq.gz"),
        f2 = os.path.join(config["scratchdir"], "trimmed_fastqs/{sample_name}_R2_trimmed.fastq.gz")
    output:
        f1 = os.path.join(config["scratchdir"], "trimmed_fastqs/{sample_name}_R1_trimmed.fastq"),
        f2 = os.path.join(config["scratchdir"], "trimmed_fastqs/{sample_name}_R2_trimmed.fastq")
    shell:
        """
        gunzip {input.f1};
        gunzip {input.f2}
        """
'''
# add --readFilesCommand zcat  to read in zipped fastqs so you dont have to unzip
rule aling_STAR_femlaes:
    input:
        f1 = os.path.join(config["scratchdir"], "trimmed_fastqs/{sample_name}_R1_trimmed.fastq.gz"),
        f2 = os.path.join(config["scratchdir"], "trimmed_fastqs/{sample_name}_R2_trimmed.fastq.gz")
    output:
        os.path.join(config["scratchdir"], "bams/{sample_name}_XX_sorted_Aligned.sortedByCoord.out.bam")
    params:
        STAR_index = config["STAR_index_XX"],
        gtf = config["gene_gtf"],
        stem = os.path.join(config["scratchdir"], "bams/{sample_name}_XX_sorted_")
    shell:
        "STAR --runMode alignReads --genomeDir {params.STAR_index} "
        "--sjdbGTFfile {params.gtf}  --winAnchorMultimapNmax 100 "
        "--outFilterMultimapNmax 100 "
        "--readFilesIn {input.f1} {input.f2} --readFilesCommand zcat "
        "--outSAMtype BAM SortedByCoordinate "
        "--outFileNamePrefix {params.stem} --runThreadN 8"


# --readStrand Reverse flag is no longer avaliable in STAR
rule aling_STAR_males:
    input:
        f1 = os.path.join(config["scratchdir"], "trimmed_fastqs/{sample_name}_R1_trimmed.fastq.gz"),
        f2 = os.path.join(config["scratchdir"], "trimmed_fastqs/{sample_name}_R2_trimmed.fastq.gz")
    output:
        os.path.join(config["scratchdir"], "bams/{sample_name}_XY_sorted_Aligned.sortedByCoord.out.bam")
    params:
        STAR_index = config["STAR_index_XY"],
        gtf = config["gene_gtf"],
        stem = os.path.join(config["scratchdir"], "bams/{sample_name}_XY_sorted_")
    shell:
        "STAR --runMode alignReads --genomeDir {params.STAR_index} "
        "--sjdbGTFfile {params.gtf}  --winAnchorMultimapNmax 100 "
        "--outFilterMultimapNmax 100 "
        "--readFilesIn {input.f1} {input.f2} --readFilesCommand zcat "
        "--outSAMtype BAM SortedByCoordinate "
        "--outFileNamePrefix {params.stem} --runThreadN 8"


#-------------------------------------------------------------------------------
# Main step 6 #
# Mark duplicates
#-------------------------------------------------------------------------------
rule MarkDups_females:
    input:
        os.path.join(config["scratchdir"], "bams/{sample_name}_XX_sorted_Aligned.sortedByCoord.out.bam")
    output:
        bam = os.path.join(config["scratchdir"], "bams/{sample_name}_XX_sorted_mkdups.bam"),
        metrics = os.path.join(config["scratchdir"], "stats/{sample_name}.XX.picard_mkdup_metrics.txt")
    shell:
        """
        picard -Xmx14g MarkDuplicates I={input} O={output.bam} M={output.metrics} VALIDATION_STRINGENCY=LENIENT
        """

rule MarkDups_males:
    input:
        os.path.join(config["scratchdir"], "bams/{sample_name}_XY_sorted_Aligned.sortedByCoord.out.bam")
    output:
        bam = os.path.join(config["scratchdir"], "bams/{sample_name}_XY_sorted_mkdups.bam"),
        metrics = os.path.join(config["scratchdir"], "stats/{sample_name}.XY.picard_mkdup_metrics.txt")
    shell:
        """
        picard -Xmx14g MarkDuplicates I={input} O={output.bam} M={output.metrics} VALIDATION_STRINGENCY=LENIENT
        """


#-------------------------------------------------------------------------------
# Main step 7 #
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

rule index_bam_females:
    input:
        os.path.join(config["scratchdir"], "bams/{sample_name}_XX_sorted_mkdups_rdgrp.bam")
    output:
        os.path.join(config["scratchdir"], "bams/{sample_name}_XX_sorted_mkdups_rdgrp.bam.bai")
    shell:
        "samtools index {input}"

rule index_bam_males:
    input:
        os.path.join(config["scratchdir"], "bams/{sample_name}_XY_sorted_mkdups_rdgrp.bam")
    output:
        os.path.join(config["scratchdir"], "bams/{sample_name}_XY_sorted_mkdups_rdgrp.bam.bai")
    shell:
        "samtools index {input}"


#-------------------------------------------------------------------------------
# Main step 7 #
# Generate TE count table in uniq mode (this will cound all reads - both unique
# and multi mapped reads)
#-------------------------------------------------------------------------------
# TEtranscripts does the DE analysis for you. See command below. I will run
# TEcount which I believe just outputs the count file
# TEtranscripts --format BAM --mode multi -t RNAseq1.bam RNAseq2.bam -c CtlRNAseq1.bam CtlRNAseq.bam --project sample_nosort_test
# TEcount --sortByPos --format BAM --GTF /data/CEM/shared/public_data/references/GENCODE/gencode.v29.annotation.gtf \
# --TE /scratch/amtarave/TE_expression/te_gtf/GRCh38_GENCODE_rmsk_TE.gtf --mode multi -b bams/OBG0115-2-011_XX_sorted_rdgrp.bam \
# --stranded reverse --project results_TEcounts_multi/OBG0115-2-011_XX
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
        TEcount --sortByPos --format BAM --GTF {params.genegtf} --TE {params.tegtf} --mode multi -b {input} --stranded reverse --project {params.outdir}
        """
'''
rule TE_counts_uniq_females:
    input:
        os.path.join(config["scratchdir"], "bams/{sample_name}_XX_sorted_mkdups_rdgrp.bam")
    params:
        tegtf = config["TE_gtf"],
        genegtf = config["gene_gtf"],
        outdir = os.path.join(config["scratchdir"], "results_TEcounts_uniq/{sample_name}_XX_sorted_mkdups_rdgrp_TEcount_uniq")
    output:
        os.path.join(config["scratchdir"], "results_TEcounts_uniq/{sample_name}_XX_sorted_mkdups_rdgrp_TEcount_uniq.cntTable")
    shell:
        """
        TEcount --sortByPos --format BAM --GTF {params.genegtf} --TE {params.tegtf} --mode uniq -b {input} --stranded reverse --project {params.outdir}
        """
'''


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
        TEcount --sortByPos --format BAM --GTF {params.genegtf} --TE {params.tegtf} --mode multi -b {input} --stranded reverse --project {params.outdir}
        """

'''
rule TE_counts_uniq_males:
    input:
        os.path.join(config["scratchdir"], "bams/{sample_name}_XY_sorted_mkdups_rdgrp.bam")
    params:
        tegtf = config["TE_gtf"],
        genegtf = config["gene_gtf"],
        outdir = os.path.join(config["scratchdir"], "results_TEcounts_uniq/{sample_name}_XY_sorted_mkdups_rdgrp_TEcount_uniq")
    output:
        os.path.join(config["scratchdir"], "results_TEcounts_uniq/{sample_name}_XY_sorted_mkdups_rdgrp_TEcount_uniq.cntTable")
    shell:
        """
        TEcount --sortByPos --format BAM --GTF {params.genegtf} --TE {params.tegtf} --mode uniq -b {input} --stranded reverse --project {params.outdir}
        """
'''
