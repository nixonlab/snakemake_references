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

"""
rule star_index_generic:
    input:
        genome_gz = 'databases/sequences/{genome_name}/genome.fna.gz',
        annot_gz = 'databases/annotations/{annot_name}/transcripts.gtf.gz'
    output:
        'databases/indexes/STAR_{genome_name}_{annot_name}/SA',
        'databases/indexes/STAR_{genome_name}_{annot_name}/SAindex',
        'databases/indexes/STAR_{genome_name}_{annot_name}/Genome'
    log:
        'databases/indexes/STAR_{genome_name}_{annot_name}/Log.out'
    params:
        sjdbOverhang=100,
        genomeDir = lambda wc, out: os.path.dirname(out[0])
    conda: "../envs/star.yaml"
    threads: snakemake.utils.available_cpu_count()
    shell:
        '''
mkdir -p {params.genomeDir}
STAR\
 --runThreadN {threads}\
 --runMode genomeGenerate\
 --genomeDir {params.genomeDir}\
 --genomeFastaFiles {input.genome_gz}\
 --sjdbGTFfile {input.annot_gz}\
 --sjdbOverhang {params.sjdbOverhang}
        '''
"""


rule star_index_local:
    input:
        genome_gz = 'databases/sequences/{genome_name}/genome.fna.gz',
        annot_gz = 'databases/annotations/{annot_name}/transcripts.gtf.gz'
    output:
        directory("databases/indexes/STAR_{genome_name}_{annot_name}"),
        'databases/indexes/STAR_{genome_name}_{annot_name}/SA',
        'databases/indexes/STAR_{genome_name}_{annot_name}/SAindex',
        'databases/indexes/STAR_{genome_name}_{annot_name}/Genome'
    log:
        'databases/indexes/STAR_{genome_name}_{annot_name}/Log.out'
    params:
        sjdbOverhang=100,
        genomeDir = lambda wc, out: os.path.dirname(out[1])
    conda: "../envs/star.yaml"
    threads: snakemake.utils.available_cpu_count()
    shell:
        '''
tdir=$(mktemp -d {config[local_tmp]}/{rule}.XXXXXX)
STAR\
 --runThreadN {threads}\
 --runMode genomeGenerate\
 --genomeDir $tdir\
 --genomeFastaFiles {input.genome_gz}\
 --sjdbGTFfile {input.annot_gz}\
 --sjdbOverhang {params.sjdbOverhang}
 
 cp $tdir {params.genomeDir}
 chmod 0444 {params.genomeDir}/*
 chmod 0555 {params.genomeDir}
        '''