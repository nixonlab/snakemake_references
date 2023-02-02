# Resources

## Remote file checksum downloads

## GCA_000001405.15_GRCh38

Obtain the file list and md5sums for NCBI reference genome:

```bash
curl -o checksums/GCA_000001405.15_GRCh38.md5checksums.txt\
   ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/md5checksums.txt
```


```bash
curl -o checksums/GCA_009914755.4_T2T-CHM13v2.0.md5checksums.txt\
    ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/914/755/GCA_009914755.4_T2T-CHM13v2.0/md5checksums.txt
```


```bash
curl -o checksums/T2T.CHM13.assemblies.analysis_set.README.txt  https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/README.txt
l=$(($(grep -nm1 '^md5sum$' checksums/T2T.CHM13.assemblies.analysis_set.README.txt | cut -d':' -f1) + 1))
tail -n+$l checksums/T2T.CHM13.assemblies.analysis_set.README.txt > checksums/T2T.CHM13.assemblies.analysis_set.md5checksums.txt  
```


### `GRCh83.d1.vd1_virus_decoy.txt`

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
