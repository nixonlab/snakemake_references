import json

localrules: build_remotefile_db
rule build_remotefile_db:
    output:
        'resources/remotefile_db.json'
    script:
        '../scripts/build_remotefile_db.py'


def getparams_download_remote_db(wildcards):
    with open(rules.build_remotefile_db.output[0], 'r') as dbfh:
        db = json.load(dbfh)
    return dict(db[wildcards.f])

rule download_remotefile:
    """ Downloads a remote file and checks the md5sum.
        Filenames, URLs and md5 checksums are stored in 
        resources/remotefiles_db.txt
    """
    output:
        'databases/remotefiles/{f}'
    input:
        rules.build_remotefile_db.output
    params:
        getparams_download_remote_db
    shell:
        '''
curl -L {params[0][url]} > {output[0]}
echo {params[0][md5]}  {output[0]} | md5sum -c -
        '''

