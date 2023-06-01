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

star_index_files = [
 'chrLength.txt',
 'chrNameLength.txt',
 'chrName.txt',
 'chrStart.txt',
 'exonGeTrInfo.tab',
 'exonInfo.tab',
 'geneInfo.tab',
 'Genome',
 'genomeParameters.txt',
 'SA',
 'SAindex',
 'sjdbInfo.txt',
 'sjdbList.fromGTF.out.tab',
 'sjdbList.out.tab',
 'transcriptInfo.tab',
]


rule star_index_local:
    input:
        genome_gz = 'databases/sequences/{genome_name}/genome.fna.gz',
        annot_gz = 'databases/annotations/{annot_name}/transcripts.gtf.gz'
    output:
        protected(directory("databases/indexes/STAR_{genome_name}_{annot_name}")),
        expand('databases/indexes/STAR_{{genome_name}}_{{annot_name}}/{basename}',
            basename = star_index_files
        )
    params:
        sjdbOverhang=100,
        genomeDir = lambda wildcards, output: os.path.dirname(output[1])
    conda: "../envs/star.yaml"
    threads: snakemake.utils.available_cpu_count()
    shell:
        '''
tdir=$(mktemp -d {config[local_tmp]}/{rule}.XXXXXX)
echo "temporary directory: $tdir"
unpigz < {input.genome_gz} > $tdir/genome.fna &
unpigz < {input.annot_gz} > $tdir/transcripts.gtf &
wait

mkdir -p $tdir/starindex
STAR\
 --runThreadN {threads}\
 --runMode genomeGenerate\
 --genomeDir $tdir/starindex\
 --genomeFastaFiles $tdir/genome.fna\
 --sjdbGTFfile $tdir/transcripts.gtf\
 --sjdbOverhang {params.sjdbOverhang}

mkdir -p {output[0]}
cp $tdir/starindex/* {output[0]}
        '''
