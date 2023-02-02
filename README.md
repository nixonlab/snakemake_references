# Snakemake Reference Database


[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)


The goal of this workflow is to have reproducible reference data for use in
analysis pipelines. The generated reference data - such as indexes and data
summaries - will depend on versioned, publicly available data sources
and commonly used open-access software tools. This workflow can be rapidly 
deployed to different computing platforms.

In the future it might be good to host pre-built versions of the reference
data on a cloud platform.

## `remotefiles`

Contains original files downloaded from remote locations. The 
[`download_remote`](workflow/rules/remote.smk)
rule downloads targets from a URL and checks the file's MD5SUM against 
a value stored in `config['remotefiles']`. At runtime, targets are added with
the destination file name as the key and the value being a dictionary with the
URL and MD5SUM for the file.
For example:

```python
config['remotefiles'] = {
    'ACC12345.genome_sequence.fasta.gz': {
        'url': 'ftp://ftp.public_resource.org/ACC12345/genome_sequence.fasta.gz',
        'md5': 'ea9832679ae2a42b6d645de352913307',
    },
    'ACC987654.genes.gff.gz': {
        'url': 'ftp://ftp.other_resource.gov/ACC987654/genes.gff.gz',
        'md5': 'a743a98262a42245645de529c13ba6a4',
    }
}
```

Remote files can be added to [config.yaml](config/config.yaml) or can be
loaded from an MD5 checksum file. Typically, full URL paths will need
to be constructed from the base URL - see example in
[`download_remote`](workflow/rules/remote.smk).
MD5 checksum files can be added to this repository and stored in 
`resources/checksums`; documentation for downloading checksum files should be
added to [`resources/README.md`](resources/README.md).

## `sequences`

Contains reference sequence data.

Genomic reference sequences are nucleic acid (DNA) sequences in FASTA format.
Files have the extension `*.fasta.gz`, `*.fa.gz`, or `*.fna.gz`  and are 
compressed using [`bgzip`](http://www.htslib.org/doc/bgzip.html).
It is recommended that you index the compressed file using 
[`samtools faidx`](http://www.htslib.org/doc/faidx.html)
and create a sequence dictionary using
[`picard CreateSequenceDictionary`](https://gatk.broadinstitute.org/hc/en-us/articles/9570414410267-CreateSequenceDictionary-Picard-).

## `annotations`

Annotation data for sequences.

## `indexes`

Indexes for aligning sequencing reads to genomes.

Aligners include STAR, bowtie2, hisat2


### Install Snakemake

```bash
mamba create -c bioconda -c conda-forge --name snakemake snakemake snakedeploy
```

