#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Following Christian Diener's neat example of how to find the H2S genes
Thanks Christian!!

Created on 19/12/22

"""


from cobra.io import read_sbml_model

model = read_sbml_model("../1_GEMs/S975C400.batch2.xml") # P. dorei
#model = read_sbml_model("../1_GEMs/S367C239.batch2.xml") # Megamonas funiformis
#model = read_sbml_model("../1_GEMs/S702C844.batch2.xml") # Alistipes onderdonkii


model

# Okay as the first step let's find all possible producers.
# There are several ways to do it but we will use the reaction query functionality.
# #So we need to write a small function what will return True if a reaction can produce H2S and False otherwise.

h2s = model.metabolites.h2s_c

is_h2s_producer = lambda rxn: (h2s in rxn.reactants and rxn.lower_bound < 0) or (h2s in rxn.products and rxn.upper_bound > 0)

h2s_producers = model.reactions.query(is_h2s_producer)
h2s_producers

#Okay, we could also get the genes for those:
{r.name: r.genes for r in h2s_producers}


#And finally let's get the fluxes and quickly get the H2S contributions.

from micom.qiime_formats import load_qiime_medium
from cobra.flux_analysis import pfba

medium = load_qiime_medium("0_diet/not_used/western_diet_gut_carveme.qza")
medium.reaction = medium.reaction.replace("_m$", "_e", regex=True)
medium.index = medium.reaction
rids = [r.id for r in model.reactions]
model.medium = medium.loc[medium.reaction.isin(rids), "flux"]
sol = pfba(model)
sol

print(model.metabolites.h2s_c.summary(solution=sol))


# Okay seems to come from a single reaction. We can also create a dataframe with the results though it is a bit trivial here.
import pandas as pd

results = pd.DataFrame()
results["id"] = [r.id for r in h2s_producers]
results["reaction"] = [r.name for r in h2s_producers]
results["genes"] = [",".join(str(g) for g in r.genes) for r in h2s_producers]
results["flux"] = [r.metabolites[h2s] * sol.fluxes[r.id] for r in h2s_producers]

pd.set_option('display.max_columns', 5)

results

#save it:
results.to_csv("h2s_producer_genes_MAGxx.csv")

