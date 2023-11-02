# -*- coding: utf-8 -*-

rule extract_cellranger_ref:
    input:
        'databases/remotefiles/refdata-{refver}.tar.gz'
    output:
        'databases/indexes/cellranger/refdata-{refver}/reference.json',
        'databases/indexes/cellranger/refdata-{refver}/fasta/genome.fa',
        'databases/indexes/cellranger/refdata-{refver}/fasta/genome.fa.fai',
        'databases/indexes/cellranger/refdata-{refver}/genes/genes.gtf',
        'databases/indexes/cellranger/refdata-{refver}/star/SAindex'
    wildcard_constraints:
        refver = '(gex-GRCh38-2020-A|gex-mm10-2020-A|gex-GRCh38-and-mm10-2020-A|cellranger-vdj-GRCh38-alts-ensembl-7.1.0)'
    shell:
        '''
mkdir -p databases/indexes/cellranger
tar -xzvf {input[0]} -C databases/indexes/cellranger
        '''


rule star_index_GRCh38_2020A_F4HSA_gencode32:
    input:
        genome_fna = 'databases/indexes/cellranger/refdata-gex-GRCh38-2020-A/fasta/genome.fa',
        hiv_fna = 'databases/localfiles/HIV1.H4-FSA.fna',
        genome_gtf = 'databases/indexes/cellranger/refdata-gex-GRCh38-2020-A/genes/genes.gtf'
    output:
        protected(directory("databases/indexes/STAR_GRCh38_2020A_F4HSA_gencode32")),
        expand('databases/indexes/STAR_GRCh38_2020A_F4HSA_gencode32/{basename}',
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
cat {input.genome_fna} {input.hiv_fna} > $tdir/genome.fna
cat {input.genome_gtf} > $tdir/transcripts.gtf

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
