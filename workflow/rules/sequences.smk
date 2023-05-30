# -*- coding: utf-8 -*-

from os.path import basename

"""
Add GCA_000001405.15_GRCh38 analysis sets to config['remotefiles']

See https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/README_analysis_sets.txt
"""
_baseftp = '''ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids'''
for l in open('resources/checksums/GCA_000001405.15_GRCh38.md5checksums.txt', 'r'):
    _md5,_fn = l.split()
    _fn = basename(_fn)
    config['remotefiles'][_fn] = {
        'url': f'{_baseftp}/{_fn}',
        'md5': _md5,
    }

"""
Add T2T analysis set to config['remotefiles']
"""
_basehtml = '''https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set'''
for l in open('resources/checksums/T2T.CHM13.assemblies.analysis_set.md5checksums.txt', 'r'):
    _md5,_fn = l.split()
    _fn = basename(_fn)
    config['remotefiles'][_fn] = {
        'url': f'{_basehtml}/{_fn}',
        'md5': _md5,
    }

# for k in config['remotefiles']:
#    print(f'{k}\n{config["remotefiles"][k]["url"]}')


rule setup_ncbi_analysis_set:
    input:
        'databases/remotefiles/GCA_000001405.15_GRCh38_{aset}_analysis_set.fna.gz'
    output:
        'databases/sequences/GCA_000001405.15_GRCh38_{aset}_analysis_set/genome.fna.gz',
        'databases/sequences/GCA_000001405.15_GRCh38_{aset}_analysis_set/genome.fna.gz.gzi',
        'databases/sequences/GCA_000001405.15_GRCh38_{aset}_analysis_set/genome.fna.gz.fai',
        'databases/sequences/GCA_000001405.15_GRCh38_{aset}_analysis_set/genome.dict'
    conda:
        '../envs/utils.yaml'
    shell:
        '''
mkdir -p $(dirname {output[0]})        
gunzip -c {input} | bgzip > {output[0]}
samtools faidx {output[0]}
picard CreateSequenceDictionary -R {output[0]}
        '''


rule setup_gdc_analysis_set:
    input:
        'databases/remotefiles/GRCh38.d1.vd1.fa.tar.gz'
    output:
        'databases/sequences/GRCh38.d1.vd1/genome.fna.gz',
        'databases/sequences/GRCh38.d1.vd1/genome.fna.gz.gzi',
        'databases/sequences/GRCh38.d1.vd1/genome.fna.gz.fai',
        'databases/sequences/GRCh38.d1.vd1/genome.dict'
    conda:
        '../envs/utils.yaml'
    shell:
        '''
mkdir -p $(dirname {output[0]})
tar -Oxzf {input} | bgzip > {output[0]}
samtools faidx {output[0]}
picard CreateSequenceDictionary -R {output[0]}
        '''

rule setup_T2T_analysis_set:
    input:
        'databases/remotefiles/chm13v2.0{aset}.fa.gz'
    output:
        'databases/sequences/chm13v2.0{aset}/genome.fna.gz',
        'databases/sequences/chm13v2.0{aset}/genome.fna.gz.gzi',
        'databases/sequences/chm13v2.0{aset}/genome.fna.gz.fai',
        'databases/sequences/chm13v2.0{aset}/genome.dict'
    wildcard_constraints:
        aset = '(_maskedY|_maskedY_rCRS|_noY)*'
    conda:
        '../envs/utils.yaml'
    shell:
        '''
mkdir -p $(dirname {output[0]})
gunzip -c {input} | bgzip > {output[0]}
samtools faidx {output[0]}
picard CreateSequenceDictionary -R {output[0]}
        '''
