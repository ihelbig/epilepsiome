# set directory
setwd("~/Downloads")

#activate libraries
library(tidyverse)
library(ggplot2)

# load files

## exac tibble includes mis_z and pLI for all genes in the human genome
exac <- read_tsv("ftp://ftp.broadinstitute.org/pub/ExAC_release/release1/manuscript_data/forweb_cleaned_exac_r03_march16_z_data_pLI.txt.gz")

## ccr tibble includes genes with constrained coding regions above the 99% percentile
ccr <- read_tsv("https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-018-0294-6/MediaObjects/41588_2018_294_MOESM3_ESM.txt")

#Show distribution of mis_z scores
ggplot(data=exac, aes(exac$mis_z)) + geom_histogram()
quantile(exac$mis_z, prob = c(0.90, 0.95, 0.99))

#Fix Excel contamination, 6-Mar instead of MARCH6
ccr[ccr$gene == "6-Mar",1] <- "MARCH6"

#Select relevant columns from exac file
exac %>% select(gene,mis_z,pLI) -> exac_lim

#filter most constrained genes
ccr_exac <- merge(ccr,exac_lim, by.x="gene",by.y="gene",all.x=T, all.y=T)
ccr_exac <- rename(ccr_exac, CCRs = number_of_99th_percentile_CCRs)
ccr_exac %>% filter(mis_z >3.688 | pLI >= 0.9) -> upper_exac

#Build "tornado plot", showing mis_z score per gene versus pLI
#genes with ccr > 99% in different color, size of point indicates number of CCRs
ggplot(ccr_exac, aes(x=mis_z, y=pLI)) + 
  geom_point(alpha = 0.7, aes(colour = factor(CCRs) , size = CCRs ) ) +
  geom_point(alpha = 0.1, aes(colour = factor(CCRs), size = 3 ) )

#chose lowest quantiles to find genes that may be missed by focussing on pLI and mis_z
quantile(exac$mis_z,c(0.01,0.1,0.05,0.1,0.25,0.5)) -> low_quant_mis_z
quantile(exac$pLI,c(0.01,0.1,0.05,0.1,0.25,0.5)) -> low_quant_pLI
ccr_exac_low <- merge(ccr,exac_lim, by.x="gene",by.y="gene",all.x=T, all.y=F)
ccr_exac_low %>% filter(
  (mis_z <= low_quant_mis_z[6] | is.na(mis_z) )
  & (pLI <= low_quant_pLI[6]) | is.na(pLI) ) -> lower_exac
