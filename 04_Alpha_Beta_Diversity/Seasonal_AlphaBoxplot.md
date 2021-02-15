Plot alpha diversity indices in a boxplot
================

After checking the rarefaction and species richness, we can also check
other alpha diversity measurements. Here, we use the Simpson index,
which we calculate with the `vegan` package and visualise the diversity
in a boxplot.

## Load Data

``` r
rm(list = ls())

library(ggplot2)
library(vegan)
library(RColorBrewer)
library(ggpubr)
library(reshape2)
library(viridis)

#setwd("04_Alpha_Diversity/")

OTU_Table = as.data.frame(read.csv("../00_Data/05_Oomycota_Seasonal_OTU_Table_min-freq-1172_transposed_withMetadata.tsv", 
                     header = T, 
                     sep = "\t", 
                     stringsAsFactors = T))

SampleMetadata = OTU_Table[,1:17]
Microhabitat = SampleMetadata$Microhabitat
TreeSpecies = SampleMetadata$TreeSpecies
OTU_Table = OTU_Table[,18:ncol(OTU_Table)]
```

## Calculate Alpha Diversity

To run the diversity analyses, simply load the table and specify the
`index` - in this case: The Simpson index. Then convert it into a
dataframe and add the metadata and group.

``` r
#simpson = diversity(OTU_Table, index = "simpson")
#df = as.data.frame(simpson)
shannon = diversity(OTU_Table, index = "shannon")
df = as.data.frame(shannon)
df$richness = specnumber(OTU_Table)
df$evenness = df$shannon/log(df$richness)
df$Microhabitat = Microhabitat
rownames(df) = SampleMetadata$SampleID
df$Group = "Oomycota"
df$TreeSpecies = TreeSpecies
df$Season = SampleMetadata$Season
df$SamplingDescription = SampleMetadata$SamplingDescription
df$SamplingDescription = factor(df$SamplingDescription, levels=c('Autumn 2017','Spring 2018','Autumn 2018','Spring 2019'))

df_melted = melt(df)
```

    ## Using Microhabitat, Group, TreeSpecies, Season, SamplingDescription as id variables

``` r
df_melted$SamplingDescription = factor(df_melted$SamplingDescription, levels=c('Autumn 2017','Spring 2018','Autumn 2018','Spring 2019'))
```

## Plot the Figure

Now we put the diverity measurements into a habitat specific context. It
can be easiest visualised in a boxplot:

``` r
g = ggplot(df_melted, aes(x = Microhabitat, y = value, fill = Microhabitat)) + 
  stat_boxplot(geom = "errorbar", width = 0.2, show.legend = F) +
  geom_boxplot(show.legend = F) + 
  #geom_point(aes(shape = TreeSpecies), size = 4) +
  scale_fill_manual(values = c(viridis(7, direction = -1), 
                               "#8e8878", "#524640"), 
                    limits = c("Arboreal Soil", "Bark", "Deadwood", 
                                "Fresh Leaves", "Hypnum", "Lichen", 
                                "Orthotrichum", "Leaf Litter", "Soil")) + 
  scale_x_discrete(limits = c("Arboreal Soil", "Bark", "Deadwood", 
                                "Fresh Leaves", "Hypnum", "Lichen", 
                                "Orthotrichum", "Leaf Litter", "Soil")) + 
  #scale_y_continuous(limits = c(0, 1)) +
  theme_minimal() + 
  labs(y = "Alpha Diversity", 
       x = "Microhabitat") +
  #geom_dotplot(aes(x = Microhabitat, y = simpson, fill = TreeSpecies), 
  #             binaxis = "y", stackdir = "center", binwidth = 1, 
  #             binpositions = "all", dotsize = 0.01, 
  #             position = "stack") +
  theme(axis.text=element_text(size=8), 
        axis.title=element_text(size=10, face = "bold"), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14, face = "bold"), 
        legend.position = "right", 
        legend.direction = "vertical", 
        strip.text = element_text(size=9, face = "bold"), 
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5), 
        axis.text.x = element_text(angle = 45, hjust = 1)) + 
  #facet_grid(rows = vars(SamplingDescription), cols = vars(variable), scales = "free")
  facet_grid(variable ~ SamplingDescription, 
             scales = "free", switch = "y")

g
```

![](Seasonal_AlphaBoxplot_files/figure-gfm/OomycotaAlphaBoxPlot-1.png)<!-- -->

``` r
g$theme$axis.text.x$angle = 45
g$theme$axis.text.x$hjust = 1
g$theme$axis.text.x$vjust = 1

g$theme$plot.subtitle$hjust = 0


ggsave("AlphaBoxplot.png", plot = g, 
       device = "png", dpi = 300, width = 18, height = 12, 
       units = "cm")
ggsave("AlphaBoxplot.jpeg", plot = g, 
       device = "jpeg", dpi = 300, width = 18, height = 12, 
       units = "cm")
ggsave("AlphaBoxplot.pdf", plot = g, 
       device = "pdf", dpi = 300, width = 18, height = 12, 
       units = "cm")
```
