# -*- coding: utf-8 -*-

rule repeat_annotation_retro38:
    input:
        'databases/remotefiles/retro.hg38.v1.gtf',
        'databases/remotefiles/retro.hg38.v1.tsv.gz'
    output:
        'databases/annotations/retro.hg38.v1.gtf',
        'databases/annotations/retro.hg38.v1.tsv.gz'
    shell:
        '''
cp {input[0]} {output[0]}
cp {input[1]} {output[1]}
        '''
