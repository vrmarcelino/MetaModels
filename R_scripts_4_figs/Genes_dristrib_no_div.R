#### Test if certain genes are differentially present in healthy and CD microbiomes
#  remove genes with less than 10 occurrences:
# This script uses spp presence/abdsence only, probably informative analysis as we expect H2S production to be flexible at the community scale, 
# with less abundant taxa potentially able to produce a fair amount of it (S.G.)


#updated: 21 Jan 2023

rm(list=ls(all=TRUE))

#Import required library
library(ape) # for pcoa function
library(ggplot2)
library(funrar)
library(gridExtra)

####################################################
### Input and parse tables
genes_raw <- read.csv("0_input/all_genes.csv")
rownames(genes_raw) <- genes_raw[,1]
genes_raw$MAG <- NULL

# This table was generated with 0_make_table.R script, includes the total + net production & consumption of H2S
input_parsed_table_fp <- "0_input/metad_and_metab_h2s.csv"
metad_and_metab_study <- read.csv(input_parsed_table_fp)
rownames(metad_and_metab_study) <- metad_and_metab_study[,1]

# species composion (to map MAGs to samples)
spp_composition_fp <- "0_input/He_etal_CD_abundances.csv"
spp_conposition <- read.csv(spp_composition_fp)
rownames(spp_conposition) <- spp_conposition[,1]
spp_conposition$Taxonomy <- NULL
spp_conposition$X <- NULL

### test: make it presence/absence:
spp_conposition <- as.data.frame(spp_conposition>0)*1L

####################################################
### make a matrix of Sample x Gene abundance:

### simpler way (Thanks Jamie!):
t_gene_table <- t(genes_raw)
gene_by_sample <- t_gene_table %*% as.matrix(spp_conposition)


####################################
########### Test and Plot ###########
#####################################


### Merge genes and metadata
genes_t <- t(gene_by_sample)
genes_met <-merge(metad_and_metab_study, genes_t, by="row.names")
rownames(genes_met) <- genes_met[,1]

## Add spp diversity:
spp_conposition["spp_diversity",] <- colSums(spp_conposition != 0) # calculate the number of spp in each sample
spp_div <- t(spp_conposition["spp_diversity",])

genes_met_spp_div <-merge(genes_met, spp_div, by="row.names")
rownames(genes_met_spp_div) <- genes_met_spp_div[,1]

# drop the unwanted columns:
gene_df_0 <- genes_met_spp_div[-c(1:9)]
gene_df_1 <- gene_df_0[-c(2:6)]

# use only samples containing the min number of spp in Healthy cohort (99 species)
### test removing this step!!!!
gene_df2 <- gene_df_1[gene_df_1$spp_diversity > 99,]

# remove genes with less than 10 occurrences:
gene_df2["gene_occurrence",] <- colSums(gene_df2 != 0)
gene_df3 <- gene_df2[,gene_df2["gene_occurrence",] > 9]

# add diagnosis back:
gene_df3$Diagnosis <- gene_df2$Diagnosis


#remove occurrence count:
gene_df <- gene_df3[!(row.names(gene_df2) %in% c("gene_occurrence")),]

# save it:
#write.csv(gene_df, "4.0_gene_by_sample_filt.csv")


## We would probably find a significant difference for most genes due to spp diversity differences,
# so It makes sense to use a lm to account for the effect of spp diversity here:

genes <- names(gene_df)[1:40]

store_Ps <- data.frame(matrix(nrow = length(genes), ncol = 2))
rownames(store_Ps) <- genes
colnames(store_Ps) <- c("p_value", "estimate")

# calculate the bonferroni-corrected p-value:
bon_p <- 0.05 / length(genes)
bon_p 


#log transform
# add 0.1 to avoid inf values:
gene_df[1:(ncol(gene_df)-1)] <- gene_df[1:(ncol(gene_df)-1)] + 0.1
# then log transform:
gene_df[1:(ncol(gene_df)-1)] <- log(gene_df[1:(ncol(gene_df)-1)])



#loop through genes, get the p value of the correlation:
for (i in 1:length(genes)){
  target_gene <- genes[i]
  model1 <- lm(paste0(target_gene," ~ Diagnosis"), data=gene_df)
  summary(model1)
  targeted_p <- summary(model1)$coefficients[2, 4]
  targeted_estimate <- summary(model1)$coefficients[2, 1]
  
  store_Ps[target_gene,"p_value"] <- targeted_p
  store_Ps[target_gene,"estimate"] <- targeted_estimate
}

# add significance if kruskal-wallis p is < than bin_p:
store_Ps$sig_after_bonferroni <- "ns"
store_Ps$sig_after_bonferroni <- ifelse(as.numeric(store_Ps$p_value) < bon_p, "*", store_Ps$sig_after_bonferroni)

# add directio: if estimate is neg, then direction is "more abundant in healthy":
store_Ps$direction <- "more abundant in CD"
store_Ps$direction <- ifelse(as.numeric(store_Ps$estimate) < 0, "more abundant in healthy", store_Ps$direction)

print(store_Ps)


## save table:
write.csv(store_Ps, "4.0_lm_genes_no_div_pres_abs.csv")


