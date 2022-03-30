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
