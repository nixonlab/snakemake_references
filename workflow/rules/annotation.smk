include: "rules/gencode.smk"
include: "rules/repeat_annotation.smk"

rule extract_gtf_gz:
    input:
        'databases/remotefiles/{f}.gtf.gz'
    output:
        'databases/annotations/{f}.gtf'    
    shell:
        '''
gunzip -c {input} > {output}
        '''


rule tabix_gtf:
    input:
        'databases/annotations/{f}.gtf'
    output:
        'databases/annotations/{f}.gtf.gz',
        'databases/annotations/{f}.gtf.gz.tbi'
    conda:
        '../envs/utils.yaml'
    shell:
        '''
(grep "^#" {input}; grep -v "^#" {input} | sort -t"`printf '\t'`" -k1,1 -k4,4n) | bgzip > {output[0]}
tabix -p gff {output[0]}       
        '''


# rule bgzip_fasta:
#     input:
#         'databases/sequences/{prefix}.{ext, fna|fa}'
#     output:
#         'databases/sequences/{prefix}.{ext}.gz',
#         'databases/sequences/{prefix}.{ext}.gz.fai'
#     conda:
#         '../envs/utils.yaml'
#     shell:
#         '''bgzip < {input} > {output[0]} && samtools faidx {output[0]}'''


# rule extract_gdc38:
#     input:
#         'databases/remotefiles/GRCh38.d1.vd1.fa.tar.gz'
#     output:
#         'databases/sequences/GRCh38.d1.vd1.fa',
#         'databases/sequences/GRCh38.d1.vd1.fa.fai',
#         'databases/sequences/GRCh38.d1.vd1.dict'
#     conda:
#         '../envs/utils.yaml'
#     shell:
#         '''
# mkdir -p $(dirname {output[0]})
# tar -Oxzf {input} > {output[0]}
# samtools faidx {output[0]}
# picard CreateSequenceDictionary R={output[0]} O={output[2]}
#         '''


# rule extract_ncbi38:
#     input:
#         'databases/remotefiles/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz'
#     output:
#         'databases/sequences/ncbi38.fna',
#         'databases/sequences/ncbi38.fna.fai',
#         'databases/sequences/ncbi38.dict'
#     conda:
#         '../envs/utils.yaml'
#     shell:
#         '''
# mkdir -p $(dirname {output[0]})
# gunzip -c {input} > {output[0]}
# samtools faidx {output[0]}
# picard CreateSequenceDictionary R={output[0]} O={output[2]}
#         '''


