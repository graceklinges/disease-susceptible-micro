library("phyloseq")
library("plyr")
library("dplyr")
library("tidyr")
library("ggplot2")
library("reshape2")
library("ggnetwork")
library("ggthemes")
library("cowplot")

setwd("~/Mote_nutrient_experiment/data/")
load(file = "ps_rare_g50.RData")
ps <- ps_rarefied
rm(ps_rarefied)
sample_sums(ps) #should be 5098
mapfile = "~/Mote_nutrient_experiment/data/Mapping-file-full-renamed_50.txt" #confirm correct mapping file associated
map = import_qiime_sample_data(mapfile)
sample_data(ps) <- map


nutrient <- merge_samples(ps, "Nutrient_weeks", fun = mean)
nutrient_rel <- transform(nutrient, "compositional")
nsamples(nutrient_rel)

ps_genus <- tax_glom(nutrient_rel, "Genus")
ps_melt <- psmelt(ps_genus)

# summarize by Genus, sort, get top taxa
genus <- ps_melt %>% group_by(Genus) %>% summarise(average = mean(Abundance))
genus_sort <- arrange(genus, desc(average))
top <- genus_sort[which(genus_sort$average > 0.01),]
top <- top$Genus
ps_melt_top <- ps_melt[ps_melt$Genus %in% top,]
top

write.csv(ps_melt_top, file = "ps_melt_top.csv") #had to manually fix broken variables in excel...
rm(ps_melt_top)
ps_melt_top <- read.csv(file = "ps_melt_top_fixed.csv")

# Get genera with mean relative abundance >0.01 across all samples 

ps_melt_top$Genus <- factor(ps_melt_top$Genus, levels = c("Unclassified_ASV_2", "Aquarickettsia", "Family_Francisellaceae", "Unclassified_ASV_6",
                                                          "Spirochaeta_2"))

colorblind_pallette = c("#0072B2", "#CC3300", "#F4B90C","#660066", "#56B4E9","#FFFFFF")

nozero<- subset(ps_melt_top, Abundance > 0)
p <- ggplot(nozero, aes(x = Sample, y = reorder(OTU, Abundance))) +
  geom_point(aes(size=Abundance, fill = Genus), shape = 21) + #shape = 21 tells it to do colored circles
  facet_grid(~reorder(Nutrient, Nutrient_order), scales = "free", space = "free", switch = "y") + #splits into resistant v susceptible categories
  theme_facet() +
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 1), 
        axis.text.y = element_blank(), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 0.8)) +  
  scale_fill_manual(values = colorblind_pallette, 
                    breaks = c("Aquarickettsia", "Unclassified_ASV_2", "Unclassified_ASV_6", "Family_Francisellaceae", "Spirochaeta_2"),
                    labels = c("Aquarickettsia", "Unclassified ASV 2", "Unclassified ASV 6", "Family Francisellaceae", "Spirochaeta")) +  
  guides(fill = guide_legend(override.aes = list(size=4))) + #makes legend circles bigger
  labs(x= "", y = "", fill = "Genus", size = "Relative \nAbundance")
p
ggsave(filename="~/Mote_nutrient_experiment/plots/g50_bubble.png", plot=p, device="png",  width = 12, height = 3.5, dpi=500)
ggsave(filename="~/Mote_nutrient_experiment/plots/g50_bubble.svg", plot=p, device="svg",  width = 12, height = 3.5, dpi=500)

