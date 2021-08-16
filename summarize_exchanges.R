#### Summarize the metabolic exchnages across samples or biomes
### Input files from MetModels_merge_exchange_tables.py script, which processes MICOM results.

#### created: 18 May 2021
#### edited: 2 Aug 2021


rm(list=ls(all=TRUE))

### Input file
input_min_fluxes <- "0_input/minimal_fluxes_exchange_merged_test_dataset.csv"


########### PRODUCTION ###########

### Load MICOM results (needs to do thsi again for each table)
exchanges_samples <- read.csv(input_min_fluxes)

#### Calculate production (sum only positive values)
exchanges_samples[exchanges_samples < 0] <- 0 # keep only positive
sum_positive <- aggregate(exchanges_samples[,2:(ncol(exchanges_samples)-1)], list(exchanges_samples$sample), sum)
rownames(sum_positive) <- sum_positive[,1]
sum_positive <- sum_positive[-1]

write.csv(sum_positive, "1_Summary_exchanges/total_production.csv")



########### CONSUMPTION ###########
exchanges_samples <- read.csv(input_min_fluxes)

#### sum only negative values
cols <- sapply(exchanges_samples, is.numeric) # get only numeric columns
exchanges_samples[,cols][exchanges_samples[,cols] > 0] <- 0 # replace positive values with zero

sum_neg <- aggregate(exchanges_samples[,2:(ncol(exchanges_samples)-1)], list(exchanges_samples$sample), sum)
rownames(sum_neg) <- sum_neg[,1]
sum_neg <- sum_neg[-1]

write.csv(sum_neg, "1_Summary_exchanges/total_consumption.csv")


########### NET Production/Consumption ###########
exchanges_samples <- read.csv(input_min_fluxes)

net_prod <- aggregate(exchanges_samples[,2:(ncol(exchanges_samples)-1)], list(exchanges_samples$sample), sum)
rownames(net_prod) <- net_prod[,1]
net_prod <- net_prod[-1]

write.csv(net_prod, "1_Summary_exchanges/net_production_and_consumption.csv")




