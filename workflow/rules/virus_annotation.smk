#! /usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd

gdc_viruses = pd.read_csv(
    'resources/GRCh83.d1.vd1_virus_decoy.txt',
    sep='\t',
    header=0,
    index_col="Abbreviation",
    dtype=pd.StringDtype()
)
gdc_viruses['acc'] = np.where(pd.notnull(gdc_viruses["RefSeq"]), gdc_viruses["RefSeq"], gdc_viruses["GenBank"])
gdc_viruses['acc'] = gdc_viruses['acc'].str.split(pat='.', expand=True)[0]

rule build_vd1_table:
    output:
        'resources/GRCh83.d1.vd1_virus_decoy.txt'
    params:
        table_url = 'https://gdc.cancer.gov/files/public/file/GRCh83.d1.vd1_virus_decoy.txt'
    shell:
        '''
curl -s -L {params.table_url} |\
  sed 's/$/\n/' | sed 's/\r/\n/g' | sed -r '/^\s*$/d' | \
  sed \
    -e 's/CMV ( HHV-5)/CMV/' \
    -e 's/EBV (HHV-4)/chrEBV/' \
    -e 's/"HHV-8, KSHV"/KSHV/' \
  > GRCh83.d1.vd1_virus_decoy.txt        
        '''


# rule virus_genome_from_GRCh38_d1_vd1:
#     """ Extract genome sequence (FASTA) for one virus from GRCh38_d1_vd1"""
#     output:
#         "databases/sequences/viruses/{newref}.fna.gz",
#         "databases/sequences/viruses/{newref}.fna.gz.gzi",
#         "databases/sequences/viruses/{newref}.fna.gz.fai",
#         "databases/sequences/viruses/{newref}.dict"
#     input:
#         "databases/sequences/GRCh38.d1.vd1.fa.gz",
#         "databases/sequences/GRCh38.d1.vd1.fa.gz.fai"
#     conda: '../envs/utils.yaml'
#     shell:
#         '''
# mkdir -p $(dirname {output[0]})
# samtools faidx {input[0]} {wildcards.newref} | bgzip > {output[0]}
# samtools faidx {output[0]}
# picard CreateSequenceDictionary -R {output[0]}
#         '''

rule setup_vd1_fasta:
    """ Extract genome sequences (FASTA) for all viruses from GRCh38.d1.vd1"""
    output:
        "databases/sequences/vd1/genome.fna.gz",
        "databases/sequences/vd1/genome.fna.gz.gzi",
        "databases/sequences/vd1/genome.fna.gz.fai",
        "databases/sequences/vd1/genome.dict"
    input:
        "databases/sequences/GRCh38.d1.vd1/genome.fna.gz",
        "databases/sequences/GRCh38.d1.vd1/genome.fna.gz.fai",
        'resources/GRCh83.d1.vd1_virus_decoy.txt'
    params:
        virus_reg = lambda wc: ' '.join(gdc_viruses.index)
    conda: '../envs/utils.yaml'
    shell:
        '''
mkdir -p $(dirname {output[0]})
samtools faidx {input[0]} {params.virus_reg} | bgzip > {output[0]}
samtools faidx {output[0]}
picard CreateSequenceDictionary -R {output[0]}        
        '''

rule remote_gff3_genbank:
    """ Download GFF3 annotation for accession """
    output:
        'databases/remotefiles/genbank/{accession}.gff3'
    shell:
        '''
curl --create-dirs -o {output} "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id={wildcards.accession}"
        '''

rule gdc_virus_gtf:
    """ Create GTF for a virus in GDC virus decoy """
    output:
        'databases/annotations/gdc.vd1/{newref}.gtf'
    input:
        lambda wc: f'databases/remotefiles/viruses/{gdc_viruses.loc[wc.newref]["acc"]}.gff3'
    conda: '../envs/gffread.yaml'
    shell:
        '''
gffread --force-exons -FET {input} | workflow/scripts/rename_reference.py {wildcards.newref} > {output}
        '''

rule gdc_viruses_gtf:
    """ GTF for all GDC viruses """
    output:
        "databases/annotations/gdc.vd1.gtf"
    input:
        expand('databases/annotations/gdc.vd1/{newref}.gtf', newref=gdc_viruses.index)
    shell:
        '''
cat {input} > {output}        
        '''


rule virus_assembly_gtf:
    """ 
        Use remote_gff3_genbank instead of this because not all viruses are 
        in the assembly database and I am unsure of the GTF. 
    """
    output:
        'databases/remotefiles/viruses/{accession}_genomic.gtf.gz'
    conda: '../envs/entrez.yaml'
    shell:
        '''
# Determine if accession is in "Assembly" database
searchXML=$(esearch -db assembly -query {wildcards.accession})
asmcount=$(echo $searchXML | xtract -pattern ENTREZ_DIRECT -element Count)

# Download from "Assembly" FTP
if [[ "$asmcount" -gt "0" ]]; then
    echo "Found $asmcount result(s) for {wildcards.accession} in assembly database"
    summaryXML=$(echo $searchXML | esummary)
    gbkURL=$(echo $summaryXML | xtract -pattern DocumentSummary -element FtpPath_GenBank)
    refURL=$(echo $summaryXML | xtract -pattern DocumentSummary -element FtpPath_RefSeq)
    if [[ ! -z "$refURL" ]]; then
        url=$refURL/$(basename $refURL)_genomic.gtf.gz
        echo "Downloading from $url"
        curl -o {output[0]} $url || echo "curl fail $?" 
    elif [[ ! -z "$gbkURL" ]]; then
        url=$gbkURL/$(basename $gbkURL)_genomic.gtf.gz
        echo "Downloading from $url"
        curl -o {output[0]} $url || echo "curl fail $?"
    fi    
fi
        '''

# if [[ ! -e {output[0]} ]]; then
#     rawgff=$(dirname {output[0]})/{wildcards.accession}_genomic.gff3
#     url="https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id={wildcards.accession}"
#     curl -o $rawgff $url
#
#     echo "#gtf-version 2.2" | gzip > {output[0]}
#     head $rawgff | grep "^#" | sed "s/^##/#!/" | gzip >> {output[0]}
#     gffread --keep-comments --force-exons -F -E -T $rawgff | gzip >> {output[0]}
# fi


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
# efetch -db nuccore -id NC_006273 -format gff3




