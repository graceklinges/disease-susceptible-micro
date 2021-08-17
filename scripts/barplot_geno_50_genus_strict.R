library("phyloseq")
library("data.table")
library("plyr")
library("dplyr")
library("tidyr")
library("ggplot2")
library("reshape2")
library("indicspecies")
library("ggnetwork")
library("ape")
library("microbiome")
library("ggthemes")
library("cowplot")
library("ggsignif")

setwd("~/Mote_nutrient_experiment/data")

# load the rarefied data table with all species
load(file = "ps_rare_g50.RData") #renamed, rarefied ps object
ps <- ps_rarefied
rm(ps_rarefied)
sample_sums(ps) #should be 5098
ps = subset_samples(ps, SampleID != "Negative")
mapfile = "Mapping-file-full-renamed.txt"
map = import_qiime_sample_data(mapfile)
sample_data(ps) <- map

# transform to relative abundance
ps_rel <- transform(ps, "compositional")
sample_sums(ps_rel) #1
# melt the data at the Genus level
ps_rel_genus_melt <- ps_rel %>%
  tax_glom(taxrank = "Genus") %>%
  psmelt()
head(ps_rel_genus_melt)
write.csv(ps_rel_genus_melt, file = "relative_abundance.csv")


# Get genera with mean relative abundance >0.01 across all samples 
gen_sum <- ps_rel_genus_melt %>% group_by(Genus) %>% dplyr::summarise(Aver = mean(Abundance))
gen_sub <- gen_sum[which(gen_sum$Aver > 0.01),]
names <- gen_sub$Genus
names

# Replace genus with <0.01 abundance with "NA"
ps_rel_genus_melt$genus <- ps_rel_genus_melt$Genus

ps_rel_genus_melt$genus[ps_rel_genus_melt$genus != "Aquarickettsia" &
                           ps_rel_genus_melt$genus != "Unclassified_ASV_2" &
                           ps_rel_genus_melt$genus != "Unclassified_ASV_6" &
                           ps_rel_genus_melt$genus != "Family_Francisellaceae" &
                           ps_rel_genus_melt$genus != "Spirochaeta_2"] <- NA

#levels(ps_rel_fam_melt_50$Nutrient) <- new_order

# plot
bar_genus = ggplot(ps_rel_genus_melt, aes(x = reorder(NW_T0_merge, Exposure_weeks), y=Abundance)) + 
  geom_bar(stat="identity", position="fill", aes(fill = reorder(genus, Abundance))) +
  scale_fill_manual(values=c("#0072B2", "#CC3300", "#ffff00","#660066", "#56B4E9","#FFFFFF"), 
                    breaks = c("Aquarickettsia", "Unclassified_ASV_2", "Unclassified_ASV_6","Family_Francisellaceae", "Spirochaeta_2"),
                    labels = c("Aquarickettsia", "Unclassified ASV 2", "Unclassified ASV 6","Family Francisellaceae", "Spirochaeta"),
                    na.value = "transparent") +
  facet_grid(~reorder(Nutrient, Nutrient_order), scales = "free_x", space = "free_x") +
  theme_bw() +
  #theme(legend.text = element_text(face = "italic")) +
  ylab("Relative Abundance") +
  xlab("Treatment Weeks") +
  theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
  labs(fill = "Genus", size = "Genus")
bar_genus
#exporting as svg so I can italicize Aquarickettsia in Illustrator
ggsave(filename="~/Mote_nutrient_experiment/plots/Geno_50_NW_genus.svg", plot=bar_genus, device="svg", width = 7, height = 4, dpi=500)
