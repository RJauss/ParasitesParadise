Seasonal\_TernaryPlot
================

## Load data

Ternary plots are a good way to partition OTU abundances on three
habitats. In this case, we focus on the concatenated canopy habitats
vs. soil vs. leaf litter. First we load the data:

``` r
rm(list = ls())

library(ggplot2)
library(ggtern)
library(lemon)
library(viridis)
library(funrar)

OTU_Table = read.csv("../00_Data/05_Oomycota_Seasonal_OTU_Table_min-freq-1172_transposed_withMetadata.tsv", 
                     header=T, stringsAsFactors=TRUE, sep="\t")
TaxonomyTable = read.csv("../00_Data/04_Oomycota_Seasonal_OTU_ContingencyTable_filtered_sequences_NCBI-nt_blasted_sorted_BestHit_TaxonomyTable_Oomycetes_sequences_vsearch-ITS1-BestHit_AnnotationRefined_noPipe_withFunction_noNA.tsv", 
                         header = F, stringsAsFactors = T, sep = "\t", row.names = 1)
TaxonomyTable = TaxonomyTable[,c(1:4, 7:8)]
colnames(TaxonomyTable) = c("Order", "Family", "Genus", 
                            "Species", "Lifestyle", "Substrate")
TaxonomyTable$Species = gsub("Genus_NoHit_Species_NoHit", "undetermined", 
                             TaxonomyTable$Species)
TaxonomyTable$Species = gsub("Species_NoHit", "sp.", 
                             TaxonomyTable$Species)
TaxonomyTable$Species = paste0("OTU", rownames(TaxonomyTable),
                               "_", TaxonomyTable$Species)
TaxonomyTable$Order = gsub("Order_NoHit", "Undetermined", 
                           TaxonomyTable$Order)
TaxonomyTable$Order = gsub("incertae_sedis", "Incertae sedis", 
                           TaxonomyTable$Order)

species = OTU_Table[,18:ncol(OTU_Table)]
species_mat = as.matrix(species)
SampleMetadata = as.data.frame(OTU_Table[,1:17])
```

## Aggregate

Now we aggregate by soil, leaf litter and the concatenated canopy
habitats. To do so we generate a vector where all canopy habitats are
labeled as “canopy”, so we can easily aggregate afterwards:

``` r
TestVector = vector()
for(i in as.vector(SampleMetadata$Microhabitat)){
  if(i %in% c("Arboreal Soil", "Bark", "Deadwood", "Fresh Leaves", 
              "Hypnum", "Orthotrichum", "Lichen")){
  TestVector = c(TestVector, "Canopy")
} else if(i == "Soil"){
  TestVector = c(TestVector, "Soil")
} else {
  TestVector = c(TestVector, "Leaf Litter")
}
}

SampleMetadata$Ternary = TestVector

Terndata = aggregate(species_mat, 
            by = list(SampleMetadata$Ternary), 
            FUN = mean)
rownames(Terndata) = Terndata$Group.1
Terndata = Terndata[,-1]
Terndata = as.data.frame(t(Terndata))
colnames(Terndata) = c("Canopy", "LeafLitter", "Soil")
Terndata = cbind(Terndata, TaxonomyTable)
```

## Plot

Now we plot the dataframe with the ggtern package:

``` r
t = ggtern(data = Terndata, aes(x = LeafLitter, y = Canopy, z = Soil)) + 
  geom_point(aes(color = Lifestyle), size = 2.5) +
  scale_color_manual(values = c("lightblue4", "indianred4", 
                               "khaki", "grey"), 
                    limits = c("hemibiotroph", "obligate biotroph", 
                               "saprotroph", "undetermined")) +
  theme_arrowdefault() + 
  theme(legend.position = c(0.92, 0.12),
        legend.justification = c(0.9, 0.1), 
        axis.text=element_text(size=8), 
        axis.title=element_text(size=6, face = "bold", hjust = 0.5), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 10, face = "bold"), 
        #legend.position = "right", 
        #legend.direction = "vertical", 
        strip.text = element_text(size=10, face = "bold"), 
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) + 
  labs(x="Leaf\nLitter",y="Canopy",z="Soil") +
  guides(color = guide_legend(override.aes = list(size = 2.5))) +
  #facet_grid(cols = vars(Order))
  facet_wrap(vars(Order), ncol = 3, nrow = 3)


ggsave("TernaryPlot.pdf", plot = t, 
       device = "pdf", dpi = 300, width = 18, height = 18, 
       units = "cm")
ggsave("TernaryPlot.png", plot = t, 
       device = "png", dpi = 300, width = 18, height = 18, 
       units = "cm")
ggsave("TernaryPlot.jpeg", plot = t, 
       device = "jpeg", dpi = 300, width = 18, height = 18, 
       units = "cm")

t
```

![](Seasonal_Ternary_files/figure-gfm/Plot%20ggtern-1.png)<!-- -->

## Peronosporales & Pythiales

Additionally we can focus on the two biggest oomycete orders and
partition their abundances on the sampling events as well:

``` r
SampleMetadata$TernarySeason = 
  ifelse(grepl("Canopy", SampleMetadata$StratumPerSampling),
         as.character(SampleMetadata$StratumPerSampling), 
         as.character(SampleMetadata$MicrohabitatPerSampling))

TerndataSeason = aggregate(species_mat, 
            by = list(SampleMetadata$TernarySeason), 
            FUN = mean)
rownames(TerndataSeason) = TerndataSeason$Group.1
TerndataSeason = TerndataSeason[,-1]
TerndataSeason = as.data.frame(t(TerndataSeason))
#colnames(Terndata) = c("Canopy", "LeafLitter", "Soil")
#Terndata = cbind(Terndata, TaxonomyTable)
#reshape2::melt(TerndataSeason, measure.vars=grep("^Canopy", colnames(TerndataSeason)), 
#     id.vars = c("Canopy Autumn 2017", "Canopy Spring 2018", 
#                 "Canopy Autumn 2018", "Canopy Spring 2019"))

TerndataSeason_long = data.frame(Canopy =c(TerndataSeason$`Canopy Autumn 2017`, 
                                          TerndataSeason$`Canopy Spring 2018`, 
                                          TerndataSeason$`Canopy Autumn 2018`, 
                                          TerndataSeason$`Canopy Spring 2019`))
TerndataSeason_long$LeafLitter = c(TerndataSeason$`Leaf Litter Autumn 2017`, 
                                   TerndataSeason$`Leaf Litter Spring 2018`, 
                                   TerndataSeason$`Leaf Litter Autumn 2018`, 
                                   TerndataSeason$`Leaf Litter Spring 2019`)
TerndataSeason_long$Soil = c(TerndataSeason$`Soil Autumn 2017`, 
                                   TerndataSeason$`Soil Spring 2018`, 
                                   TerndataSeason$`Soil Autumn 2018`, 
                                   TerndataSeason$`Soil Spring 2019`)
TerndataSeason_long$Sampling = rep(c("Autumn 2017", "Spring 2018", 
                                     "Autumn 2018", "Spring 2019"), each = 375)
TaxonomyTable_long = rbind(TaxonomyTable, TaxonomyTable, 
                           TaxonomyTable, TaxonomyTable)

TerndataSeason_long = cbind(TerndataSeason_long, TaxonomyTable_long)

ColSumsPerSampling = function(a){as.vector(colSums(species_mat[SampleMetadata$SamplingDescription == a,]))}

TerndataSeason_long$Abundance = 
  c(ColSumsPerSampling("Autumn 2017"), 
    ColSumsPerSampling("Spring 2018"),
    ColSumsPerSampling("Autumn 2018"),
    ColSumsPerSampling("Spring 2019"))

TerndataSeason_long_PeronoPythiales = TerndataSeason_long[TerndataSeason_long$Order %in% c("Peronosporales", "Pythiales"),]

TerndataSeason_long_PeronoPythiales$Sampling =
  factor(as.character(TerndataSeason_long_PeronoPythiales$Sampling), 
         levels = c("Autumn 2017", "Spring 2018", 
                    "Autumn 2018", "Spring 2019"))

# Skip all rows which are exclusively 0
TerndataSeason_long_PeronoPythiales = 
  TerndataSeason_long_PeronoPythiales[rowSums(TerndataSeason_long_PeronoPythiales[,1:3]) != 0,]
```

## Plot Peronosporales & Pythiales

``` r
p = ggtern(data = TerndataSeason_long_PeronoPythiales, aes(x = LeafLitter, y = Canopy, z = Soil)) + 
  geom_point(aes(color = Lifestyle), size = 2.5) +
  scale_color_manual(values = c("lightblue4", "indianred4", 
                               "khaki", "grey"), 
                    limits = c("hemibiotroph", "obligate biotroph", 
                               "saprotroph", "undetermined")) +
  #scale_color_viridis_d(option = "magma") +
  #geom_tri_tern(stat = ..count..)
  theme_arrowdefault() + 
  theme(#legend.position = c(0.92, 0.12),
        #legend.justification = c(0.9, 0.1), 
        axis.text=element_text(size=5), 
        axis.title=element_text(size=6, face = "bold", hjust = 0.5), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 10, face = "bold"), 
        legend.position = "bottom", 
        legend.direction = "vertical", 
        strip.text = element_text(size=8, face = "bold"), 
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)#, 
        #panel.spacing = unit(1.5, "lines")
        ) + 
  labs(x="Leaf\nLitter",y="Canopy",z="Soil") +
  guides(color = guide_legend(override.aes = list(size = 2.5))) +
  #facet_grid(cols = vars(Order))
  facet_grid(rows = vars(Sampling), cols = vars(Order), switch = "y")

ggsave("TernaryPlotSeason.pdf", plot = p, 
       device = "pdf", dpi = 300, width = 8.5, height = 16, 
       units = "cm")
ggsave("TernaryPlotSeason.png", plot = p, 
       device = "png", dpi = 300, width = 8.5, height = 16, 
       units = "cm")
ggsave("TernaryPlotSeason.jpeg", plot = p, 
       device = "jpeg", dpi = 300, width = 8.5, height = 16, 
       units = "cm")

p
```

![](Seasonal_Ternary_files/figure-gfm/Plot%20PeronosporalesPythiales-1.png)<!-- -->
