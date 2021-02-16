Seasonal\_TaxonomyPerSeason
================

After plotting the total taxonomic composition, we now want to find out
which groups (orders or families or genera) dominate which habitat. Here
we will plot the taxonomic composition partitioned into the habitats in
a stacked bar chart.

## Prepare Data

``` r
rm(list = ls())
library(ggplot2)
library(plyr)
library(ggpubr)
library(reshape2)
library(viridis)

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
Abundances = colSums(OTU_Table)
TAX = cbind(TAX, Abundances)
colnames(TAX) = c("OTU_Number", "Order", "Family", "Genus", "Species", "ReferenceID", "PercentID", "Lifestyle", "Substrate", "Abundance")
TAX$OTU_ID = paste0("OTU", TAX$OTU_Number, "_", TAX$Species)
TAX$Lifestyle = as.factor(TAX$Lifestyle)
TAX$Class = "Oomycota"
TAX$Order = gsub("Order_NoHit", "Undetermined", 
                           TAX$Order)
TAX$Order = gsub("incertae_sedis", "Incertae sedis", 
                           TAX$Order)
```

## Aggregate Table

Next we aggregate the table by microhabitat and order. What we then have
is a 9x8 dataframe with 9 microhabitats and 8 orders

``` r
OrderDataFunction = function(MicrohabitatPerSampling, incidence){
OrderTable = OTU_Table
colnames(OrderTable) = TAX$Order
# First aggregate by Habitat
HabitatAggregatedOrderTable = 
  aggregate(OrderTable, 
            by = list(SampleMetadata$MicrohabitatPerSampling), 
            FUN = sum)
rownames(HabitatAggregatedOrderTable) = 
  HabitatAggregatedOrderTable$Group.1
HabitatAggregatedOrderTable = 
  HabitatAggregatedOrderTable[,-1]
# if you want to check the number of OTUs instead of the amount of reads:
# this is the point to convert it into a presence/absence matrix like this:
if(incidence == T){
HabitatAggregatedOrderTable[,] = 
  ifelse(HabitatAggregatedOrderTable[,] > 0, 1, 0)
}
#Then aggregate the (transposed) table by Order
HabitatAggregatedOrderTable = 
  aggregate(t(HabitatAggregatedOrderTable), 
            by = list(TAX$Order), 
            FUN = sum)
rownames(HabitatAggregatedOrderTable) = 
  HabitatAggregatedOrderTable$Group.1
HabitatAggregatedOrderTable = 
  HabitatAggregatedOrderTable[,-1]
HabitatAggregatedOrderTable = 
  as.data.frame(HabitatAggregatedOrderTable)
# melt the table so ggplot can interpret it as a histogram
data = HabitatAggregatedOrderTable
data$Order = rownames(data)
data2 = melt(data, id.vars = "Order")
# add the Stratum
data3 = data2
data3$Stratum = ifelse(grepl("Leaf Litter|^Soil", data2$variable), 
                       "Ground", "Canopy")
# for the plotting we need to change the column names
colnames(data3) = c("Order", "Microhabitat", "Abundance", "Stratum")

data4 = data3[grepl(MicrohabitatPerSampling, data3$Microhabitat),]

data4$Microhabitat = sub(paste0(" ", MicrohabitatPerSampling), "", data4$Microhabitat)


return(data4)
}
```

## Plot stacked bar chart

Now we plot it with `ggplot` and use the Stratum as a facet:

``` r
plot_OrderData = function(OrderData){
g = ggplot(OrderData, aes(x = Microhabitat, y = Abundance, fill = Order)) + 
  geom_bar(position = position_fill(), stat = "identity") + 
  facet_grid(cols = vars(Stratum), scales = "free", space = "free") + 
  scale_y_continuous(labels = scales::percent) +
  labs(y = "Relative\nAbundance") + 
  scale_fill_viridis(discrete = T, option = "B") +
  theme_minimal() +
  theme(axis.text=element_text(size=12), 
        axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title=element_text(size=14, face = "bold"), 
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14, face = "bold"), 
        strip.text = element_text(size=12, face = "bold"), 
        panel.grid = element_blank())
return(g)
}

OrderAutumn2017 = OrderDataFunction("Autumn 2017", incidence = F)
Plot_Autumn2017 = plot_OrderData(OrderAutumn2017)

OrderSpring2018 = OrderDataFunction("Spring 2018", incidence = F)
Plot_Spring2018 = plot_OrderData(OrderSpring2018)

OrderAutumn2018 = OrderDataFunction("Autumn 2018", incidence = F)
Plot_Autumn2018 = plot_OrderData(OrderAutumn2018)

OrderSpring2019 = OrderDataFunction("Spring 2019", incidence = F)
Plot_Spring2019 = plot_OrderData(OrderSpring2019)

combi = ggarrange(Plot_Autumn2017, Plot_Spring2018, 
                  Plot_Autumn2018, Plot_Spring2019, 
                  labels = c("AUTO"), 
                  ncol = 2, nrow = 2, 
                  common.legend = T, legend = "right", 
                  align = "hv", vjust = 2.5)

OrderAutumn2017$Sampling = "Autumn 2017"
OrderSpring2018$Sampling = "Spring 2018"
OrderAutumn2018$Sampling = "Autumn 2018"
OrderSpring2019$Sampling = "Spring 2019"

OrderDataAll = rbind(OrderAutumn2017, OrderSpring2018, 
                     OrderAutumn2018, OrderSpring2019)
OrderDataAll$Microhabitat = factor(OrderDataAll$Microhabitat, 
                                   levels = c("Arboreal Soil", "Bark", 
                                              "Deadwood", "Fresh Leaves", 
                                              "Hypnum", "Lichen", 
                                              "Orthotrichum", 
                                              "Leaf Litter", "Soil"))

g = ggplot(OrderDataAll, aes(x = Microhabitat, y = Abundance, fill = Order)) + 
  geom_bar(position = position_fill(), stat = "identity") + 
  facet_grid(cols = vars(Stratum), scales = "free", space = "free") + 
  scale_y_continuous(labels = scales::percent) +
  labs(y = "Relative\nAbundance") + 
  scale_fill_viridis(discrete = T, option = "B") +
  theme_minimal() +
  geom_vline(xintercept = 7.5) +
  theme(axis.text=element_text(size=12), 
        axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title=element_text(size=14, face = "bold"), 
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14, face = "bold"), 
        strip.text = element_text(size=12, face = "bold"), 
        panel.grid = element_blank()) + 
  facet_wrap(facets = Sampling~.)

ggsave("TaxonomyPerSeason.png", plot = g, 
       device = "png", dpi = 300, width = 18, height = 12, 
       units = "cm")
ggsave("TaxonomyPerSeason.jpeg", plot = g, 
       device = "jpeg", dpi = 300, width = 18, height = 12, 
       units = "cm")
ggsave("TaxonomyPerSeason.pdf", plot = g, 
       device = "pdf", dpi = 300, width = 18, height = 12, 
       units = "cm")

g
```

![](Seasonal_TaxonomyPerSeason_files/figure-gfm/Plot%20Barchart-1.png)<!-- -->
