#! /usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd

gdc_viruses = pd.read_csv('resources/GRCh83.d1.vd1_virus_decoy.txt', sep='\t', header=0, index_col="Abbreviation", dtype=pd.StringDtype())
gdc_viruses['acc'] = np.where(pd.notnull(gdc_viruses["RefSeq"]), gdc_viruses["RefSeq"], gdc_viruses["GenBank"])
gdc_viruses['acc'] = gdc_viruses['acc'].str.split(pat='.', expand=True)[0]

# Get a subset
# gdc_viruses = pd.concat([gdc_viruses.iloc[:8], gdc_viruses.iloc[186:190]])


rule get_genbank_gtf:
    output:
        'databases/remotefiles/viruses/{accession}_genomic.gtf.gz'
    params:
        newref = lambda wc: gdc_viruses.index[gdc_viruses['acc'] == wc.accession].tolist()[0]        
    conda: '../envs/entrez.yaml'
    shell:
        '''
# Determine if accession is in "Assembly" database
res1=$(esearch -db assembly -query {wildcards.accession})
asmcount=$(echo $res1 | xtract -pattern ENTREZ_DIRECT -element Count)

# Download from Genomes FTP
if [[ "$asmcount" -gt "0" ]]; then
    baseurl=$(echo $res1 | esummary | xtract -pattern DocumentSummary -element FtpPath_GenBank)
    [[ -z "$baseurl" ]] && baseurl=$(echo $res1 | esummary | xtract -pattern DocumentSummary -element FtpPath_RefSeq)
    url=$baseurl/$(basename $baseurl)_genomic.gtf.gz
    curl -o {output[0]} $url || echo "curl fail $?"
fi

if [[ ! -e {output[0]} ]]; then
    rawgff=$(dirname {output[0]})/{wildcards.accession}_genomic.gff3
    url="https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id={wildcards.accession}"
    curl -o $rawgff $url
    
    echo "#gtf-version 2.2" | gzip > {output[0]}
    head $rawgff | grep "^#" | sed "s/^##/#!/" | gzip >> {output[0]}
    gffread --keep-comments --force-exons -F -E -T $rawgff | gzip >> {output[0]}
fi
        '''

rule get_genbank_gff3:
    output:
        'databases/remotefiles/viruses/{accession}.gff3'
    shell:
        '''
curl -o {output} "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id={wildcards.accession}"
        '''

rule gff3_to_gtf:
    input:
        'databases/remotefiles/viruses/{accession}.gff3'    
    output:
        'databases/annotations/viruses/{accession}.gtf'
    params:
        newref = lambda wc: gdc_viruses.index[gdc_viruses['acc'] == wc.accession].tolist()[0]
    conda: '../envs/gffread.yaml'
    shell:
        '''
gffread --force-exons -FET {input} | workflow/scripts/rename_reference_gtf.py {params.newref} > {output}
        '''



#         '''
# res1=$(esearch -db assembly -query {wildcards.accession})
# if echo $res1 | grep '<Count>0</Count>'; then
#     echo "no match for {wildcards.accession} {params.newref}"
#     wget -O databases/remotefiles/viruses/{wildcards.accession}_genomic.gff3 "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id={wildcards.accession}"
#     gffread --force-exons -F -E -T databases/remotefiles/viruses/{wildcards.accession}_genomic.gff3 | gzip > {output[0]}
# else
#     baseurl=$(echo $res1 | esummary | xtract -pattern DocumentSummary -element FtpPath_GenBank)
#     [[ -z "$baseurl" ]] && baseurl=$(echo $res1 | esummary | xtract -pattern DocumentSummary -element FtpPath_RefSeq)
#     url=$baseurl/$(basename $baseurl)_genomic.gtf.gz
#     wget -O {output[0]} $url
# fi
#         '''


# efetch -db nuccore -id NC_006273 -format gb


rule rename_gdc_virus:
    input:
        'databases/remotefiles/viruses/{accession}_genomic.gtf.gz'
    output:
        'databases/annotations/viruses/{accession}_genomic.gtf'
    params:
        newref = lambda wc: gdc_viruses.index[gdc_viruses['acc'] == wc.accession].tolist()[0]
    shell:
        'gunzip -c {input} | workflow/scripts/rename_reference_gtf.py {params.newref} > {output}'


rule gdc_viruses_annotation:
    output:
        "databases/annotations/gdc38_viruses.gtf"
    input:
        expand('databases/annotations/viruses/{acc}_genomic.gtf', acc=gdc_viruses['acc'])
    shell:
        'cat {input} > {output}'


rule get_gdc_virus_fasta:
    output:
        "databases/sequences/viruses/{newref}.fna"
    input:
        "databases/sequences/gdc38.fna"
    conda: '../envs/utils.yaml'        
    shell:
        'samtools faidx {input} {wildcards.newref} > {output}'


rule gdc_viruses_fasta:
    output:
        "databases/sequences/gdc38_viruses.fna"
    input:
        expand("databases/sequences/viruses/{newref}.fna", newref=gdc_viruses.index)
    shell:
        'cat {input} > {output}'


# rule gdc_viruses_annotation_from_assembly:
#     output:
#         "databases/annotations/virusesASM.annotation.gtf.gz"
#     input:
#         expand('databases/remotefiles/viruses/{acc}_genomic.gtf.gz', acc=gdc_viruses['acc'])
#     shell:
#         'cat {input} > {output}'

# one reference per GFF
# for f in databases/remotefiles/viruses/*.gff; do
#     grep -v '^#' $f | sed -r '/^\s*$/d' | cut -f1 | sort | uniq | wc -l
# done




# rule get_genbank_fasta:
#     output:
#         'databases/remotefiles/viruses/{accession}.fna'
#     shell:
#         '''
# wget -O {output} "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=fasta&id={wildcards.accession}"
#         '''

# rule gff_to_gtf:
#     input:
#         'databases/remotefiles/viruses/{accession}.gff'
#     output:
#         'databases/annotations/viruses/{accession}.gtf'
#     params:
#         newref = lambda wc: gdc_viruses.index[gdc_viruses['acc'] == wc.accession].tolist()[0]
#     conda: '../envs/gffread.yaml'
#     shell:
#         '''
# echo -e "$(grep -v '^#' {input} | head -n1 | cut -f1)\t{params.newref}" > {input}.namemap.txt
# gffread -m {input}.namemap.txt --force-exons -F -E -T {input} > {output}
# rm {input}.namemap.txt
#         '''



