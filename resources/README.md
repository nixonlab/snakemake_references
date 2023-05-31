# Resources

## TODO: Update this, we are not doing this anymore.

## Remote file checksum downloads


We use published MD5 checksums to validate the download of remote files. 
We can parse the checksum files to populate the `config['remotefiles']`
dictionary with URL and checksum information.

### GCA_000001405.15_GRCh38

Obtain the file list and md5sums for NCBI reference genome:

```bash
curl -o checksums/GCA_000001405.15_GRCh38.md5checksums.txt\
   ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/md5checksums.txt
```

### GCA_009914755.4_T2T-CHM13v2.0

```bash
curl -o checksums/GCA_009914755.4_T2T-CHM13v2.0.md5checksums.txt\
    ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/914/755/GCA_009914755.4_T2T-CHM13v2.0/md5checksums.txt
```

### T2T.CHM13.assemblies.analysis_set

The checksums for this genome are at the bottom of their README file. First,
download the README:

```bash
curl -o checksums/T2T.CHM13.assemblies.analysis_set.README.txt\
    https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/README.txt
```

Extract checksum lines from file (after line containing "md5sum")

```bash
l=$(($(grep -nm1 '^md5sum$' checksums/T2T.CHM13.assemblies.analysis_set.README.txt | cut -d':' -f1) + 1))
tail -n+$l checksums/T2T.CHM13.assemblies.analysis_set.README.txt\
    > checksums/T2T.CHM13.assemblies.analysis_set.md5checksums.txt  
```

### GENCODE

The following function returns the full URL path for a GENCODE annotation 
release. Releases without a prefix are human, and mouse annotations are 
prefixed with "M". 

```bash
gencode_url () {
    local ver=$1
    baseurl='ftp://ftp.ebi.ac.uk/pub/databases/gencode'
    if [[ $ver = M* ]]; then
        baseurl="$baseurl/Gencode_mouse/release_$ver"
    else
        baseurl="$baseurl/Gencode_human/release_$ver"
    fi
    echo $baseurl
}
```

The checksum file is called "MD5SUMS". Download checksums for a few of 
"important" GENCODE releases

+ GENCODE v22 Mar2015 *original GDC*
+ GENCODE v36 Oct2020 *updated GDC*
+ GENCODE v38 May2021 *Stellarscope development*
+ GENCODE v41 Jul2022 *Current with UCSC*
+ GENCODE v43 Feb2023 *Current as of Feb 8 2023*

```bash
for GENCODE_RELEASE in 22 36 38 41 43; do
    curl -o checksums/gencode.v${GENCODE_RELEASE}.MD5SUMS.txt\
        $(gencode_url ${GENCODE_RELEASE})/MD5SUMS
done
```

Mouse releases:

+ GENCODE vM25 Apr2020 *Last release on GRCm38*
+ GENCODE vM32 Feb2023 *Current as of Feb 8 2023 (on GRCm39)*


```bash
for GENCODE_RELEASE in M25 M32; do
    curl -o checksums/gencode.v${GENCODE_RELEASE}.MD5SUMS.txt\
        $(gencode_url ${GENCODE_RELEASE})/MD5SUMS
done
```


## `GRCh83.d1.vd1_virus_decoy.txt`

This text file contains the virus names, abbreviations, and accession numbers that are
included in the GDC reference genome. 

We download the file from the [GDC references](https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files)
site and perform the following modifications:

  - Add trailing newline `sed 's/$/\n/'`
  - Replace windows carriage returns with newline `sed 's/\r/\n/g'`
  - Remove blank lines `sed -r '/^\s*$/d'`

Also three records have reference name abbreviations that do not match
the reference file, so we make these corrections. The replacements are:

  - `CMV ( HHV-5)` &rarr; `CMV`
  - `EBV (HHV-4)` &rarr; `chrEBV`
  - `"HHV-8, KSHV"` &rarr; `KSHV`

```bash
wget -O - https://gdc.cancer.gov/files/public/file/GRCh83.d1.vd1_virus_decoy.txt |\
  sed 's/$/\n/' | sed 's/\r/\n/g' | sed -r '/^\s*$/d' | \
  sed 's/CMV ( HHV-5)/CMV/' | sed 's/EBV (HHV-4)/chrEBV/' | sed 's/"HHV-8, KSHV"/KSHV/' \
  > GRCh83.d1.vd1_virus_decoy.txt
```
