# MetaModels
Scripts to analyse metagenome-wide metabolic models


### Data mining:
Train dataset: 1143 samples. Including 608 healthy and 535 diseased samples (running).

Test dataset: 554 samples. Including 280 healthy and 274 diseased samples (used for preliminary results).

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

1.1. Build community:
`com = Community(mag_tb)`

1.2. Add media (western diet), to the community file, after diluting some nutrients to account for absorption in the intestine:

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

1.3. saves this community to file, in pickle format.


#### 2. Cooperative tradeoffs - calculate metabolic exchanges

The script [MICOM_coop_tradeoff.py](https://github.com/vrmarcelino/MetaModels/blob/main/MICOM_coop_tradeoff.py) optimizes the cooperative tradeoff, first using the western media for upper boundaries, then using the minimal media to get the metabolic exchanges. This script was largely based on the code developed for the MICOM paper:

Each sample will run thought the function:

```
def media_and_gcs(sam):
    com = load_pickle(pickles_path +"/"+ sam)

    # Get growth rates
    try:
        sol = com.cooperative_tradeoff(fraction=trade_off)
        rates = sol.members["growth_rate"].copy()
        rates["community"] = sol.growth_rate
        rates.name = sam
    except Exception:
        logger.warning("Could not solve cooperative tradeoff for %s." % sam)
        return None

    # Get the minimal medium
    med = minimal_medium(com, 0.95 * sol.growth_rate, exports=True)
    med.name = sam

    # Apply medium and reoptimize
    com.medium = med[med > 0]
    sol = com.cooperative_tradeoff(fraction=0.5, fluxes=True, pfba=False) # uses the 'classic' FBA instead of the parsimonious FBA
    fluxes = sol.fluxes
    fluxes["sample"] = sam
    return {"medium": med, "gcs": rates, "fluxes": fluxes}
```

This will generate growth rates, media (file minimal_imports_xx) and fluxes (minimal_fluxes_all_).

From the flux output, we extracted the exchange_reactions with the filtering:

```
ex_flux = fluxes.filter(regex='^EX_') # get only exchanges (starts with 'EX_')
ex_flux = ex_flux.filter(regex='e$') # remove media (media ends with 'e_m', so I want the ones that end with 'e' only)
ex_flux = ex_flux.fillna(0) #fill NANs with zeros

ex_flux = ex_flux.loc[:, (ex_flux != 0).any(axis=0)] #remove columns with all zeros
ex_flux['sample'] = fluxes['sample']
ex_flux = ex_flux.drop(index='medium') # remove medium
```
And save it as 'minimal_fluxes_exchange_xxx' file (one per sample)


### Post-processing:

1. Merge exchange tables for all samples into a single csv file with [MetModels_merge_exchange_tables.py](https://github.com/vrmarcelino/MetaModels/blob/main/MetModels_merge_exchange_tables.py)

The output of this file can be used to calculate interaction networks (where nodes are species and/or metabolites, and edges represent the flux of the metabolic exchanges).

Note - spp abundance not taken into consideration so far (except when building the community).


2. Summarize net production and consumption of specific metabolites per sample with [R script summarize_exchanges.R](https://github.com/vrmarcelino/MetaModels/blob/main/summarize_exchanges.R). This script will, for each sample, sum the total production of each metabolite, and subtract the consumption, in order to calculate the net exchange of metabolites per sample. The output is a csv table with one sample per row and one metabolite per column (negative indicating net consumption and positive indicating net production).

Note - spp abundance not taken into consideration here either, at least not yet.


### Questions / points for discussion:
1. Media (in carveme, build_community and tradeoffs)
2. Exchange reactions:

	Drop medium or meaningful info?

	com.cooperative_tradeoff -> seems less subjective to errors ("solver encountered an error infeasible") than grow workflow?
    
3. Adjust for spp. abundances

