# -*- coding: utf-8 -*-

rule star_index_gdc38_gencode38:
    input:
        fasta = 'databases/sequences/gdc38.fna',
        gtf = 'databases/annotations/gencode.v38.annotation.gtf'
    output:
        directory("databases/indexes/STAR_gdc38_gencode38"),
        "databases/indexes/STAR_gdc38_gencode38/SAindex"
    params:
        sjdbOverhang = 74
    conda: "../envs/star.yaml"
    threads: snakemake.utils.available_cpu_count()
    shell:
        '''
STAR\
 --runThreadN {threads}\
 --runMode genomeGenerate\
 --genomeDir {output[0]}\
 --outFileNamePrefix {output[0]}\
 --genomeFastaFiles {input.fasta}\
 --sjdbGTFfile {input.gtf}\
 --sjdbOverhang {params.sjdbOverhang}
        '''


rule star_index_ncbi38_gencode38:
    input:
        fasta = 'databases/sequences/ncbi38.fna',
        gtf = 'databases/annotations/gencode.v38.annotation.gtf'
    output:
        directory("databases/indexes/STAR_ncbi38_gencode38"),
        "databases/indexes/STAR_ncbi38_gencode38/SAindex"
    params:
        sjdbOverhang = 74
    conda: "../envs/star.yaml"
    threads: snakemake.utils.available_cpu_count()
    shell:
        '''
STAR\
 --runThreadN {threads}\
 --runMode genomeGenerate\
 --genomeDir {output[0]}\
 --outFileNamePrefix {output[0]}\
 --genomeFastaFiles {input.fasta}\
 --sjdbGTFfile {input.gtf}\
 --sjdbOverhang {params.sjdbOverhang}
        '''
