# Resources


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
