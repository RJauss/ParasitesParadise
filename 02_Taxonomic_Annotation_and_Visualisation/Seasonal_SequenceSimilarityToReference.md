Plot sequence similarity to reference database
================

# Sequence similarity to reference sequences

Let’s check how the sequence similarity of our OTUs is distributed.
Basically we plot the percentage of identity from our taxonomic
annotation against the number of sequences in a histogram.

## Load the data

First we need to load the data. The Taxonomy file is needed for the
percent identity, while the OTU table will be used to calculate the
total abundance of each OTU:

``` r
rm(list = ls())

library(ggplot2)
library(plyr)
library(ggpubr)

#setwd("02_Taxonomic_Annotation_and_Visualisation/")

TAX = read.csv("../00_Data/04_Oomycota_Seasonal_OTU_ContingencyTable_filtered_sequences_NCBI-nt_blasted_sorted_BestHit_TaxonomyTable_Oomycetes_sequences_vsearch-ITS1-BestHit_AnnotationRefined_noPipe_withFunction_noNA.tsv", 
                    header = F, 
                    sep = "\t", 
                    stringsAsFactors = T)
OTU_Table = as.data.frame(read.csv("../00_Data/05_Oomycota_Seasonal_OTU_Table_min-freq-1172_transposed_withMetadata.tsv", 
                     header = T, 
                     sep = "\t", 
                     stringsAsFactors = T))

SampleMetadata = OTU_Table[,1:17]
OTU_Table = OTU_Table[,18:ncol(OTU_Table)]

OTU_Table_Autumn2017 = 
  OTU_Table[SampleMetadata$SamplingDescription == "Autumn 2017",]
OTU_Table_Spring2018 = 
  OTU_Table[SampleMetadata$SamplingDescription == "Spring 2018",]
OTU_Table_Autumn2018 = 
  OTU_Table[SampleMetadata$SamplingDescription == "Autumn 2018",]
OTU_Table_Spring2019 = 
  OTU_Table[SampleMetadata$SamplingDescription == "Spring 2019",]

Abundances_Autumn2017 = colSums(OTU_Table_Autumn2017)
Abundances_Spring2018 = colSums(OTU_Table_Spring2018)
Abundances_Autumn2018 = colSums(OTU_Table_Autumn2018)
Abundances_Spring2019 = colSums(OTU_Table_Spring2019)

TAX_Autumn2017 = cbind(TAX, Abundances_Autumn2017)
TAX_Spring2018 = cbind(TAX, Abundances_Spring2018)
TAX_Autumn2018 = cbind(TAX, Abundances_Autumn2018)
TAX_Spring2019 = cbind(TAX, Abundances_Spring2019)

colnames(TAX_Autumn2017) = c("OTU_Number", "Order", "Family", "Genus", "Species", "ReferenceID", "PercentID", "Lifestyle", "Substrate", "Abundance")
colnames(TAX_Spring2018) = c("OTU_Number", "Order", "Family", "Genus", "Species", "ReferenceID", "PercentID", "Lifestyle", "Substrate", "Abundance")
colnames(TAX_Autumn2018) = c("OTU_Number", "Order", "Family", "Genus", "Species", "ReferenceID", "PercentID", "Lifestyle", "Substrate", "Abundance")
colnames(TAX_Spring2019) = c("OTU_Number", "Order", "Family", "Genus", "Species", "ReferenceID", "PercentID", "Lifestyle", "Substrate", "Abundance")
#TAX$OTU_ID = paste0("OTU", TAX$OTU_Number, "_", TAX$Species)

TAX_Autumn2017$SamplingDescription = "Autumn 2017"
TAX_Spring2018$SamplingDescription = "Spring 2018"
TAX_Autumn2018$SamplingDescription = "Autumn 2018"
TAX_Spring2019$SamplingDescription = "Spring 2019"

# Prepare the table for the Histogram

#Basically we now have two options: We can either round the identiy value and #concatenate the OTUs with the same rounded identity - which can then be plotted #in a simple bar chart, or we expand the table to make it readable for #`geom_histogram`.

#Here I stick with the second option, but nevertheless I round the percent #identity just in case you want to try the first option.


# This function rounds the value to the provided base
#mround <- function(x,base){
#    base*round(x/base)
#} 

# Here I divided by 100 to obtain values between 0 and 1
#TAX$PercentIDrounded = mround(TAX$PercentID, 1) / 100
# You can also round to the nearest e.g. 5% by putting "5" instead of "1"

# Aggregate Abundances with the same percent identity
PercentID_Autumn2017 = TAX_Autumn2017$PercentID
PercentID_Spring2018 = TAX_Spring2018$PercentID
PercentID_Autumn2018 = TAX_Autumn2018$PercentID
PercentID_Spring2019 = TAX_Spring2019$PercentID

AggregatedTAXpercent_Autumn2017 = 
  ddply(TAX_Autumn2017, "PercentID_Autumn2017", numcolwise(sum))
AggregatedTAXpercent_Autumn2017$PercentID_Autumn2017 =
  AggregatedTAXpercent_Autumn2017$PercentID_Autumn2017 / 100

AggregatedTAXpercent_Spring2018 = 
  ddply(TAX_Spring2018, "PercentID_Spring2018", numcolwise(sum))
AggregatedTAXpercent_Spring2018$PercentID_Spring2018 =
  AggregatedTAXpercent_Spring2018$PercentID_Spring2018 / 100

AggregatedTAXpercent_Autumn2018 = 
  ddply(TAX_Autumn2018, "PercentID_Autumn2018", numcolwise(sum))
AggregatedTAXpercent_Autumn2018$PercentID_Autumn2018 =
  AggregatedTAXpercent_Autumn2018$PercentID_Autumn2018 / 100

AggregatedTAXpercent_Spring2019 = 
  ddply(TAX_Spring2019, "PercentID_Spring2019", numcolwise(sum))
AggregatedTAXpercent_Spring2019$PercentID_Spring2019 =
  AggregatedTAXpercent_Spring2019$PercentID_Spring2019 / 100

# Expand the table. 
# This means that if e.g. 1000 sequences have the percent identity of 70.5 we print 70.5 1000 times
AggregatedTAXpercent_expanded_Autumn2017 = 
  AggregatedTAXpercent_Autumn2017[rep(row.names(AggregatedTAXpercent_Autumn2017), AggregatedTAXpercent_Autumn2017$Abundance), c(1,4)]
AggregatedTAXpercent_expanded_Autumn2017$SamplingDescription = "Autumn 2017"
colnames(AggregatedTAXpercent_expanded_Autumn2017) = 
  c("PercentID", "Abundance", "SamplingDescription")

AggregatedTAXpercent_expanded_Spring2018 = 
  AggregatedTAXpercent_Spring2018[rep(row.names(AggregatedTAXpercent_Spring2018), AggregatedTAXpercent_Spring2018$Abundance), c(1,4)]
AggregatedTAXpercent_expanded_Spring2018$SamplingDescription = "Spring 2018"
colnames(AggregatedTAXpercent_expanded_Spring2018) = 
  c("PercentID", "Abundance", "SamplingDescription")

AggregatedTAXpercent_expanded_Autumn2018 = 
  AggregatedTAXpercent_Autumn2018[rep(row.names(AggregatedTAXpercent_Autumn2018), AggregatedTAXpercent_Autumn2018$Abundance), c(1,4)]
AggregatedTAXpercent_expanded_Autumn2018$SamplingDescription = "Autumn 2018"
colnames(AggregatedTAXpercent_expanded_Autumn2018) = 
  c("PercentID", "Abundance", "SamplingDescription")

AggregatedTAXpercent_expanded_Spring2019 = 
  AggregatedTAXpercent_Spring2019[rep(row.names(AggregatedTAXpercent_Spring2019), AggregatedTAXpercent_Spring2019$Abundance), c(1,4)]
AggregatedTAXpercent_expanded_Spring2019$SamplingDescription = "Spring 2019"
colnames(AggregatedTAXpercent_expanded_Spring2019) = 
  c("PercentID", "Abundance", "SamplingDescription")

HistoData_reads = rbind(AggregatedTAXpercent_expanded_Autumn2017,
                  AggregatedTAXpercent_expanded_Spring2018, 
                  AggregatedTAXpercent_expanded_Autumn2018, 
                  AggregatedTAXpercent_expanded_Spring2019)


TAX_subset_Autumn2017 = subset(TAX_Autumn2017, select = "PercentID")
TAX_subset_Autumn2017$PercentID = TAX_subset_Autumn2017$PercentID / 100
TAX_subset_Autumn2017$SamplingDescription = "Autumn 2017"
TAX_subset_Autumn2017$Abundance = Abundances_Autumn2017

TAX_subset_Spring2018 = subset(TAX_Spring2018, select = "PercentID")
TAX_subset_Spring2018$PercentID = TAX_subset_Spring2018$PercentID / 100
TAX_subset_Spring2018$SamplingDescription = "Spring 2018"
TAX_subset_Spring2018$Abundance = Abundances_Spring2018

TAX_subset_Autumn2018 = subset(TAX_Autumn2018, select = "PercentID")
TAX_subset_Autumn2018$PercentID = TAX_subset_Autumn2018$PercentID / 100
TAX_subset_Autumn2018$SamplingDescription = "Autumn 2018"
TAX_subset_Autumn2018$Abundance = Abundances_Autumn2018

TAX_subset_Spring2019 = subset(TAX_Spring2019, select = "PercentID")
TAX_subset_Spring2019$PercentID = TAX_subset_Spring2019$PercentID / 100
TAX_subset_Spring2019$SamplingDescription = "Spring 2019"
TAX_subset_Spring2019$Abundance = Abundances_Spring2019

HistoData_OTUs = rbind(TAX_subset_Autumn2017, 
                       TAX_subset_Spring2018, 
                       TAX_subset_Autumn2018, 
                       TAX_subset_Spring2019)

HistoData_OTUs = HistoData_OTUs[HistoData_OTUs$Abundance > 0,]
HistoData_OTUs$Method = "No. of OTUs"

HistoData_reads$Method = "No. of Reads"

HistoData = rbind(HistoData_reads, HistoData_OTUs)

HistoData$SamplingDescription = factor(HistoData$SamplingDescription, levels=c('Autumn 2017','Spring 2018','Autumn 2018','Spring 2019'))
HistoData$Method = factor(HistoData$Method, levels=c('No. of Reads', 'No. of OTUs'))
```

## Plot the sequence similarity

Next we add the group (Oomycota and Cercozoa, respectively) which we
will use to fill the two different histograms. Then we combine the data
and load it into `ggplot`

``` r
# The data needs to be subsetted to exclude all abundances for a percent ID of 0
g = ggplot(HistoData[HistoData$PercentID > 0, ], aes(x = PercentID)) +
  geom_histogram(aes(y = ..count..), 
               alpha = 0.8, 
               color = NA, 
               bins = 25, fill = "darkslategrey") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) + 
  scale_y_continuous(labels = scales::comma) +
  theme_minimal() +
  labs(x = "Sequence similarity to reference database", 
       y = NULL, 
       title = NULL) +
  theme(axis.text=element_text(size=5), 
        axis.title=element_text(size=7, face = "bold"), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14, face = "bold"), 
        legend.position = "right", 
        legend.direction = "vertical", 
        strip.text = element_text(size=7, face = "bold"), 
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) + 
  facet_grid(cols = vars(SamplingDescription), 
             rows = vars(Method), scales = "free", switch = "y")
  #guides(fill = guide_legend(override.aes = list(alpha = 0.3)))

g
```

![](Seasonal_SequenceSimilarityToReference_files/figure-gfm/SequenceSimilarityPlotReads-1.png)<!-- -->

``` r
ggsave("SequenceSimilarityToReference.pdf", plot = g, 
       device = "pdf", dpi = 300, width = 8.5, height = 6, 
       units = "cm")
ggsave("SequenceSimilarityToReference.jpeg", plot = g, 
       device = "jpeg", dpi = 300, width = 8.5, height = 6, 
       units = "cm")
ggsave("SequenceSimilarityToReference.png", plot = g, 
       device = "png", dpi = 300, width = 8.5, height = 6, 
       units = "cm")
```
