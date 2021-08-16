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

After creating a manifest (one per sample), containing MAGs abundance and the path to GEMs, we used snakemake to run MICOM on multiple samples
The [MICOM_Snakefile.py](https://github.com/vrmarcelino/MetaModels/tree/main/Snakemake/5_MICOM) has two rules: "build_community" and "tradeoff", explained in more detail below:

#### 1. Build community models

Rule "build_community" calls the [MICOM_build_comm_models.py](https://github.com/vrmarcelino/MetaModels/blob/main/MICOM_build_comm_models.py) script. This script will, for each sample:

1: build community:
`com = Community(mag_tb)`

2. Add media (western diet), to the community file, after diluting some nutrients to account for absorption in the intestine:

```
# check if the names match & add to the community object:
ex_ids = [r.id for r in com.exchanges]

# replace _m with _e_m (see what needs to change with the --fbc2 models)
medium = medium.replace(regex=r'_m$', value='_e_m')

medium['reaction'].isin(ex_ids).sum() # must be a large number of reactions (>100)
medium.index = medium.index.str.replace(r'_m$','_e_m', regex=True)
med = medium[medium.index.isin(ex_ids)] # exclude medium items not used by the microbiome

diet = med.flux * med.dilution # dilute nutrients absorbed in the small intestine

# This only affects the com.cooperative_tradeoff (the grow workflow is not affected)
com.medium = diet
```

Then saves this community to file in pickle format.


#### 2. Build community models



