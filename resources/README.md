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

