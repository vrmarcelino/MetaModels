# MetaModels
Scripts to analyse metagenome-wide metabolic models


### Data mining:
Train dataset: 1143 samples. Including 608 healthy and 535 diseased samples (running)
Test dataset: 554 samples. Including 280 healthy and 274 diseased samples (used for prelim results)

### Initial workflow:
Assembly (Megahit), Binning (VAMB), QC (CheckM), identification (GTDBtk), select one MAG per species, calculate bin abundances across samples.


### Generate GEMs (CarveMe)
For gap filling, we used the M8 media, which was based on Tramontano et al 2018 and is provided with the [metaGEM pipeline](https://github.com/franciscozorrilla/metaGEM).

`carve {input.bin} --gapfill M8 --mediadb {input.media} -u {params.domain} -v -o {output}`


### MICOM:


