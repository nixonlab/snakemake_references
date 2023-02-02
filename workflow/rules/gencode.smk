wildcard_constraints:
    gencode_release="M?\d+"

def gencode_url(wildcards):
    base_url = 'http://ftp.ebi.ac.uk/pub/databases/gencode'
    if re.match('^M\d+', wildcards.gencode_release):
        return f'{base_url}/Gencode_mouse/release_{wildcards.gencode_release}'
    elif re.match('\d+', wildcards.gencode_release):
        return f'{base_url}/Gencode_human/release_{wildcards.gencode_release}'


rule gencode_file_list:
    output:
        'databases/annotations/gencode.v{gencode_release}/MD5SUMS'
    params:
        urlbase = gencode_url
    shell:
        '''
curl --create-dirs -o {output[0]} -L {params.urlbase}/MD5SUMS
        '''


# def gencode_md5(wildcards, input):
#     with open(input[0]) as fh:
#         lkup = {v:k for k,v in (l.strip().split()[:2] for l in fh)}
#     if wildcards.f in lkup:
#         return lkup[wildcards.f]


rule gencode_remotefile:
    """ Downloads a remote file and checks the md5sum.
        Filenames and md5 checksums are in the output file created by 
        `gencode_file_list`
    """
    input:
        rules.gencode_file_list.output
    output:
        'databases/remotefiles/gencode.v{gencode_release}/{f}'
    params:
        urlbase = gencode_url
    shell:
        '''
curl --create-dirs -o {output[0]} -L {params.urlbase}/{wildcards.f}
md5=$(grep {wildcards.f} {input[0]} | cut -f1 -d' ')        
echo "$md5  {output[0]}" | md5sum -c -
        '''


#
#
#
# '''
# # curl -L {params.urlbase}/{wildcards.f} > {output[0]}
# echo {params.md5}  {output[0]} | md5sum -c -
# '''