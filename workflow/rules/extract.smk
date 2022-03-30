rule extract_gtf_gz:
    input:
        'databases/remotefiles/{f}.gtf.gz'
    output:
        'databases/annotations/{f}.gtf'    
    shell:
        '''
echo {resources.tmpdir}
gunzip -c {input} > {output}
        '''


rule extract_gdc38:
    input:
        'databases/remotefiles/GRCh38.d1.vd1.fa.tar.gz'
    output:
        'databases/sequences/gdc38.fna',
        'databases/sequences/gdc38.fna.fai',
        'databases/sequences/gdc38.dict'
    conda:
        '../envs/utils.yaml'
    shell:
        '''
mkdir -p $(dirname {output[0]})
tar -Oxzf {input} > {output[0]}
samtools faidx {output[0]}
picard CreateSequenceDictionary R={output[0]} O={output[2]}
        '''


rule extract_ncbi38:
    input:
        'databases/remotefiles/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz'
    output:
        'databases/sequences/ncbi38.fna',
        'databases/sequences/ncbi38.fna.fai',
        'databases/sequences/ncbi38.dict'
    conda:
        '../envs/utils.yaml'
    shell:
        '''
mkdir -p $(dirname {output[0]})
gunzip -c {input} > {output[0]}
samtools faidx {output[0]}
picard CreateSequenceDictionary R={output[0]} O={output[2]}
        '''
