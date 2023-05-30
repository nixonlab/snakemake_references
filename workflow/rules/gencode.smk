wildcard_constraints:
    gencode_release="M?\d+"

def gencode_url(wildcards):
    base_url = 'http://ftp.ebi.ac.uk/pub/databases/gencode'
    if re.match('^M\d+', wildcards.gencode_release):
        return f'{base_url}/Gencode_mouse/release_{wildcards.gencode_release}'
    elif re.match('\d+', wildcards.gencode_release):
        return f'{base_url}/Gencode_human/release_{wildcards.gencode_release}'

"""
Add GENCODE file checksums to config['remotefiles']
"""
# from glob import glob
from snakemake.io import Wildcards
# gencode_checksum_files = glob('resources/checksums/gencode.*.MD5SUMS.txt')
gencode_checksum_files = [f'resources/checksums/gencode.v{v}.MD5SUMS.txt' for v in config['gencode_versions']]
for f in gencode_checksum_files:
    rel = re.search('gencode.v(M?\d+)', f).group(1)
    _baseurl = gencode_url(Wildcards(fromdict={'gencode_release': rel}))
    for l in open(f, 'r'):
        _md5, _fn = l.split()
        config['remotefiles'][_fn] = {
            'url': f'{_baseurl}/{_fn}',
            'md5': _md5,
        }

rule gencode_checksum_file:
    output:
        'resources/checksums/gencode.v{gencode_release}.MD5SUMS.txt'
    params:
        urlbase = gencode_url
    shell:
        '''
curl --create-dirs -o {output[0]} -L {params.urlbase}/MD5SUMS
        '''


rule gencode_annotation_gtf:
    output:
        refgtf = 'databases/annotations/gencode.v{gencode_release}/gencode.v{gencode_release}.REF.annotation.gtf.gz',
        reftbi = 'databases/annotations/gencode.v{gencode_release}/gencode.v{gencode_release}.REF.annotation.gtf.gz.tbi',
        allgtf = 'databases/annotations/gencode.v{gencode_release}/gencode.v{gencode_release}.ALL.annotation.gtf.gz',
        alltbi = 'databases/annotations/gencode.v{gencode_release}/gencode.v{gencode_release}.ALL.annotation.gtf.gz.tbi',
        chroms = 'databases/annotations/gencode.v{gencode_release}/chroms.txt'

    input:
        refgtf = 'databases/remotefiles/gencode.v{gencode_release}.annotation.gtf.gz',
        allgtf = 'databases/remotefiles/gencode.v{gencode_release}.chr_patch_hapl_scaff.annotation.gtf.gz',
    conda: '../envs/utils.yaml'
    shell:
        '''
gunzip -c {input.allgtf} | grep -v '^#' | cut -f1 | uniq > {output.chroms}

(
    grep -Z "^#" {input.refgtf}
    echo "## PG: bedtools sort -g {output.chroms} -i -"
    echo "## PG: bgzip"
    grep -Zv "^#" {input.refgtf} | bedtools sort -g {output.chroms} -i -
) | bgzip > {output.refgtf}
tabix -p gff {output.refgtf}

(
    grep -Z "^#" {input.allgtf}
    echo "## PG: bedtools sort -g {output.chroms} -i -"
    echo "## PG: bgzip"    
    grep -Zv "^#" {input.allgtf} | bedtools sort -g {output.chroms} -i -
) | bgzip > {output.allgtf}
tabix -p gff {output.allgtf}
        '''

rule gencode_annotation_metadata:
    input:
        allgtf = rules.gencode_annotation_gtf.output.allgtf
    output:
        expand('databases/annotations/gencode.v{{gencode_release}}/metadata.{table}.txt.gz',
            table = ['gene_features', 'tx_features', 'exon_features',
                     'gid_gname', 'gid_gtype', 'gid_hgnc', 'gid_tid',
                     'tid_tname', 'tid_ttype', 'tid_hgnc', 'tid_gid',
                    ]
        )
    params:
        out_prefix = lambda wc: f'databases/annotations/gencode.v{wc.gencode_release}/'
    shell:
        '''
workflow/scripts/gtf_metadata.py --intfeats exon {input.allgtf} {params.out_prefix}
        '''

rule gencode_annotation_rds:
    input:
        rules.gencode_annotation_metadata.output
    output:
        expand('databases/annotations/gencode.v{{gencode_release}}/metadata.{table}.rds',
            table=['gene_features', 'tx_features', 'exon_features',
                   'gid_gname', 'gid_gtype', 'gid_hgnc', 'gid_tid',
                   'tid_tname', 'tid_ttype', 'tid_hgnc', 'tid_gid',
                   ]
        )

    conda: '../envs/rbase.yaml'
    script:
        '../scripts/metadata_rds.R'

rule gencode_annotation_complete:
    input:
        rules.gencode_annotation_gtf.output.allgtf,
        rules.gencode_annotation_gtf.output.alltbi,
        rules.gencode_annotation_rds.output
    output:
        'databases/annotations/gencode.v{gencode_release}/transcripts.gtf.gz',
        'databases/annotations/gencode.v{gencode_release}/transcripts.gtf.gz.tbi'
    shell:
        '''
chmod 0444 $(dirname {output[0]})/*
chmod 0555 $(dirname {output[0]})/*
ln -s ./$(basename {input[0]}) {output[0]}
ln -s ./$(basename {input[1]}) {output[1]}
        '''
