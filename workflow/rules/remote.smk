config['remotefiles']['gencode.v40.annotation.gtf.gz'] = {
    'url': 'http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/gencode.v40.annotation.gtf.gz',
    'md5': '14a867b82917c8c3006838c3a5053a3e'
}

rule download_remote:
    """ Downloads a remote file and checks the md5sum.
        Filenames, URLs and md5 checksums are configured in the `remotefiles` element of
        the configfile
    """
    output:
        'databases/remotefiles/{f}'
    params:
        url = lambda wildcards: config['remotefiles'][wildcards.f]['url'],
        md5 = lambda wildcards: config['remotefiles'][wildcards.f]['md5']
    shell:
        '''
curl -L {params.url} > {output[0]}
echo {params.md5}  {output[0]} | md5sum -c -
        '''
