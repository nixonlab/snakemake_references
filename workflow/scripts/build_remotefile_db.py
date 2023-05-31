#! /usr/bin/env python
import re
import requests
import json

def gencode_url(gencode_release):
    base_url = 'http://ftp.ebi.ac.uk/pub/databases/gencode'
    if re.match('^M\d+', gencode_release):
        return f'{base_url}/Gencode_mouse/release_{gencode_release}'
    elif re.match('\d+', gencode_release):
        return f'{base_url}/Gencode_human/release_{gencode_release}'

def build_remote_db(
        outfn: str,
        from_config: dict[str, dict],
        gencode_releases: list[str],
        remote_checksums: dict[str, dict]
):

    db = {}

    """ Add items in config """
    db = {**db, **from_config}

    """ Add GENCODE versions """
    for gr in gencode_releases:
        print(f'Adding remote information for GENCODE v{gr}')
        _baseurl = gencode_url(gr)
        r = requests.get(f'{_baseurl}/MD5SUMS')
        _iter = (row.split() for row in r.text.split('\n') if row)
        for _md5, _fn in _iter:
            db[_fn] = {
                'url': f'{_baseurl}/{_fn}',
                'md5': _md5,
            }

    """ Add remote checksums """
    for n, d in remote_checksums.items():
        print(f'Adding remote information for {n}')
        _baseurl = d['baseurl']
        r = requests.get(d['checksum_url'])
        _iter = re.findall(r"(?P<_md5>[a-fA-F\d]{32})\s+(?P<_fn>\S+)", r.text)
        for _md5, _fn in _iter:
            _fn0 = re.sub(r"^\./", "", _fn) # removed leading dotslash
            _fnL = re.sub(r"/", ".", _fn0)  # replace slash with dot
            assert _fnL not in db, f'ERROR: filename conflict "{_fnL}"'
            db[_fnL] = {
                'url': f'{_baseurl}/{_fn0}',
                'md5': _md5,
            }

    # for k,v in db.items():
    #     print(f'{k}:\n    url: "{v["url"]}"\n    md5: "{v["md5"]}"')

    with open(outfn, 'w') as outh:
        json.dump(db, outh, indent=4)
    return



build_remote_db(
    snakemake.output[0],
    snakemake.config['remotefiles'],
    snakemake.config['gencode_versions'],
    snakemake.config['remote_checksums']
)
