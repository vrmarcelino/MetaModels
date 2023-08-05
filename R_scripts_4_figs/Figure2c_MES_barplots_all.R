#### Check if the metabolites most affected by redundancy loss are different across diseases
# only test metabolite if it is present in > 50 samples AND in at least 15 diseased individuals:
# then produce barplots


rm(list=ls(all=TRUE))

library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyverse)
library(tidytext)
library(MetBrewer)#colors

####################################################
### Function that select the top 10 metabolites with biggest difference between health and disease
####################################################
# also removes oxigen and water
# using top 10 as default
# keeps only metabolites where there was a decrease in redundancy

select_top_MES <- function(met_table, top_X = 5) {
#function to select the top #5 metabolites 
  #(more complicated than it needs to be because I was plotting these metabolites only, 
  #but it makes sense to plot all 27 top metabolites (5 from each disease))
  
  ### Select the metabolites with significant differences only:
  sig_differences <- met_table[met_table$sig_after_bonferroni == "*",]
  
  ### Keep only the metabolites for which there was a decrease in redundancy (i.e. positive MIS)
  # most redundunancy that is negative is very small, difficult to interpret what they mean.
  sig_differences_pos <- sig_differences[sig_differences$diff_MIS > 0,]
  
  ### Simplify graph:
  # remove water and oxigen:
  names2remove <- c("o2_e", "h2o_e")
  wanted_metabolites_pre <- sig_differences_pos[!sig_differences_pos$X %in% names2remove,]
  
  # keep only top X (if there are more than that)
  if (dim(wanted_metabolites_pre)[1] > top_X) {
    print ("...more than X significant results found, keeping only the top X.")
    wanted_metabolites <- wanted_metabolites_pre[1:top_X,]} else {
      wanted_metabolites <- wanted_metabolites_pre}
  return (wanted_metabolites)
}



select_wanted_met <- function(met_table, met_list) {
#function to filter tables so they have all metabolites selected as top 5 in all diseases
    
  ### Select the metabolites with significant differences only:
  sig_differences <- met_table[met_table$sig_after_bonferroni == "*",]
  
  ### Keep only the metabolites for which there was a decrease in redundancy (i.e. positive MIS)
  # most redundunancy that is negative is very small, difficult to interpret what they mean.
  sig_differences_pos <- sig_differences[sig_differences$diff_MIS > 0,]
  
  ### Simplify graph:
  # remove water and oxigen:
  names2remove <- c("o2_e", "h2o_e")
  wanted_metabolites_pre <- sig_differences_pos[!sig_differences_pos$X %in% names2remove,]
  
  # keep only the metabolites in the selected metabolite list
  wanted_metabolites <-  wanted_metabolites_pre[wanted_metabolites_pre$X %in% met_list,]
  return (wanted_metabolites)
}





####################################################
### Input and parse tables
####################################################

# these tables were produced with the code 0.1_diff_metabolite_exc_healthy_AND_non_healthy.R
ank_raw <- read.csv("0.1_MES_tables_per_disease/0.1_mean_importance_diff_ankylosing_spondylitis.csv")
ank_raw[,"disease"] <- "ankylosing spondylitis"
ank <- select_top_MES(ank_raw)

arth_raw <- read.csv("0.1_MES_tables_per_disease/0.1_mean_importance_diff_atherosclerosis.csv")
arth_raw[,"disease"] <- "atherosclerosis"
arth <- select_top_MES(arth_raw)

beh_raw <- read.csv("0.1_MES_tables_per_disease/0.1_mean_importance_diff_behcets_disease.csv")
beh_raw[,"disease"] <- "Behcets"
beh <- select_top_MES(beh_raw)

crc_raw <- read.csv("0.1_MES_tables_per_disease/0.1_mean_importance_diff_colorectal_cancer.csv")
crc_raw[,"disease"] <- "colorectal cancer"
crc <- select_top_MES(crc_raw)

ibd_raw <- read.csv("0.1_MES_tables_per_disease/0.1_mean_importance_diff_IBD.csv")
ibd_raw[,"disease"] <- "IBD"
ibd <- select_top_MES(ibd_raw)

cir_raw <- read.csv("0.1_MES_tables_per_disease/0.1_mean_importance_diff_liver_cirrhosis.csv")
cir_raw[,"disease"] <- "liver cirrhosis"
cir <- select_top_MES(cir_raw)

cfs_raw <- read.csv("0.1_MES_tables_per_disease/0.1_mean_importance_diff_ME_CFS.csv")
cfs_raw[,"disease"] <- "ME/CFS"
cfs <- select_top_MES(cfs_raw)

nafld_raw <- read.csv("0.1_MES_tables_per_disease/0.1_mean_importance_diff_NAFLD.csv")
nafld_raw[,"disease"] <- "NAFLD"
nafld <- select_top_MES(nafld_raw)

ra_raw <- read.csv("0.1_MES_tables_per_disease/0.1_mean_importance_diff_rheumatoid_arthritis.csv")
ra_raw[,"disease"] <- "rheumatoid arthritis"
ra <- select_top_MES(ra_raw)

t2d_raw <- read.csv("0.1_MES_tables_per_disease/0.1_mean_importance_diff_type2_diabetes.csv")
t2d_raw[,"disease"] <- "type2 diabetes"
t2d <- select_top_MES(t2d_raw)



# merge tables to get a list of the top 5 metabolites from all diseases
all_MES <- rbind(ank,arth,beh,crc,ibd,cir,cfs,nafld,ra,t2d)

# These are the metabolites that will be plotted:
top_5_metabolites_from_all_dis <- unique(all_MES$X)




##### Now get a table with all selected metabolites that made the wanted list 
#(i.e. the top 5 from each disease)

ank <- select_wanted_met(ank_raw, top_5_metabolites_from_all_dis)
arth <- select_wanted_met(arth_raw, top_5_metabolites_from_all_dis)
beh <- select_wanted_met(beh_raw, top_5_metabolites_from_all_dis)
crc <- select_wanted_met(crc_raw, top_5_metabolites_from_all_dis)
ibd <- select_wanted_met(ibd_raw, top_5_metabolites_from_all_dis)
cir <- select_wanted_met(cir_raw, top_5_metabolites_from_all_dis)
cfs <- select_wanted_met(cfs_raw, top_5_metabolites_from_all_dis)
nafld <- select_wanted_met(nafld_raw, top_5_metabolites_from_all_dis)
ra <- select_wanted_met(ra_raw, top_5_metabolites_from_all_dis)
t2d <- select_wanted_met(t2d_raw, top_5_metabolites_from_all_dis)


# merge tables
all_MES <- rbind(ank,arth,beh,crc,ibd,cir,cfs,nafld,ra,t2d)


# add metabolite class
metab_names_fp <- "0_input/bigg_models_w_classes_curated_simplified_names.tsv"
metab_names <- read.csv(metab_names_fp, sep = "\t")
rownames(metab_names) <- metab_names[,1]

# add metabolite names:
names(all_MES)[names(all_MES) == 'X'] <- 'metabolite'
metab_names$metabolite <- rownames(metab_names)

top_metabs_named <- merge(all_MES, metab_names[,c("metabolite","bigg_id","SubClass")], by="metabolite")



top_metabs_named$disease <- factor(top_metabs_named$disease, levels = c("IBD","liver cirrhosis","ankylosing spondylitis",
                           "NAFLD","Behcets","ME/CFS","type2 diabetes","atherosclerosis",
                           "colorectal cancer","rheumatoid arthritis"))


top_metabs_named$SubClass <- factor(top_metabs_named$SubClass, levels = c(
  "Alcohols and polyols", "Amines", "Amino acids, peptides, and analogues",
  "Beta hydroxy acids and derivatives", "Carbohydrates and carbohydrate conjugates",
  "Carbonyl compounds", "Fatty acids and conjugates", "Non-metal phosphates",
  "Other non-metal sulfides", "Purine 2'-deoxyribonucleosides", "Pyrimidines and pyrimidine derivatives", 
  "Quaternary ammonium salts","Tricarboxylic acids and derivatives", "no sub class"))



# order metab names so the classes are together:
top_metabs_named <- arrange(top_metabs_named,SubClass, name_simple)
top_metabs_named$name_simple <- factor(top_metabs_named$name_simple, levels = c(
  "Ethanol","Putrescine","L-Valine","Reduced glutathione","L-Malate",
  "6-Phospho-D-gluconate","D-Galactose","D-Glucose","D-Ribose","Glycerol",
  "N Ribosylnicotinamide","N-Acetyl-D-glucosamine 1-phosphate","Acetaldehyde",
  "Pentanoate","Phosphate","Hydrogen sulfide","Deoxyinosine","Uracil","Choline",
  "Citrate"," 3 Hydroxypentanoic acid"," 3 hydroxyhexanoic acid","Ammonium","Fe2+",
  "Pyoverdine","Thiamin","Xanthosine"))


####################################################
### Plot
####################################################


## get colors:
n_colors <- unique(top_metabs_named$SubClass)
n_colors
wanted_colors <- met.brewer("Juarez", n=14)

pdf(file = "MES_diff_by_disease.pdf", width = 14, height = 8)

ggplot(top_metabs_named, aes(x=diff_MIS, y=name_simple)) + geom_col(aes(fill=SubClass)) + 
  facet_wrap(~disease,nrow=1,scales="free_x") + theme_bw() +
  scale_y_discrete(limits = rev) + scale_fill_manual(values=wanted_colors) +
  theme(axis.text=element_text(size=6),panel.grid.minor = element_blank())

dev.off()







