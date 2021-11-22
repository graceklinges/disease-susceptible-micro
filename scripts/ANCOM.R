library(exactRankTests)
library(dplyr)
library(ggplot2)
library(compositions)
library(phyloseq)
library(robustbase)
library(ggpubr)
library("gridExtra")
library("cowplot")

source("~/ancom_v2.1.R") #sourced from https://github.com/FrederickHuangLin/ANCOM/tree/master/scripts
setwd("~/Mote_nutrient_experiment/data")

# load unrarefied,pruned, renamed data for differential abundance analysis
load(file = "ps_rename_g50.RData.RData")
sample_sums(ps)
#remove negatives
ps = subset_samples(ps, SampleID != "Negative")
#remove outlier (discovered from dispersion plots)
ps <- subset_samples(ps, SampleID != "LP-50-2-T0")
mapfile = "~/Mote_nutrient_experiment/data/Mapping-file-full-renamed.txt"
map = import_qiime_sample_data(mapfile)
sample_data(ps) <- map
ps@sam_data
ps@tax_table

# agglomerate to genus
ps = filter_taxa(ps, function(x) sum(x > 10) > (0.2*length(x)), TRUE)
ps <- tax_glom(ps, "Genus")
ntaxa(ps) #29

ps <- subset_samples(ps, No_level_no_T0 == "Ammonium")

OTUdf <- as.data.frame(t(otu_table(ps)))
metadf <- read.delim(file = "~/Mote_nutrient_experiment/data/Mapping-file-full-renamed.txt")
colnames(metadf)[1] <- "Sample.ID"
metadf$Exposure_weeks <- as.factor(metadf$Exposure_weeks)
levels(metadf$Exposure_weeks)

#for reciprocal comparison
metadf$Exposure_weeks <- factor(metadf$Exposure_weeks, levels = c("3", "0", "6")) 
levels(metadf$Exposure_weeks)                    

# Step 1: Data preprocessing
feature_table = OTUdf; meta_data = metadf; sample_var = "Sample.ID"; group_var = NULL
out_cut = 0.05; zero_cut = 0.90; lib_cut = 1000; neg_lb = FALSE
prepro = feature_table_pre_process(feature_table, meta_data, sample_var, group_var, 
                                   out_cut, zero_cut, lib_cut, neg_lb)

#zero_cut: Taxa with proportion of zeroes greater than zero_cut are not included in the analysis.
#out_cut: observations with proportion of mixture distribution less than out_cut (or greater than 1- out_cut) 
#will be detected as outliers

feature_table = prepro$feature_table # Preprocessed feature table
meta_data = prepro$meta_data # Preprocessed metadata
struc_zero = prepro$structure_zeros # Structural zero info

# Step 2: ANCOM
main_var = "Exposure_weeks"; p_adj_method = "fdr"; alpha = 0.05
adj_formula = NULL ; rand_formula = "~ 1 | Tank"

res = ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method, 
            alpha, adj_formula, rand_formula)

resdf <- as.data.frame(res$out)

# add taxonomy
tax<-as(tax_table(ps),"matrix")
tax_cols <- c("Kingdom", "Phylum", "Class", "Order","Family","Genus", "Species")
tax<-as.data.frame(tax)
colnames(tax)
tax$taxa_id <- rownames(tax)
rownames(tax) <- NULL
dim(tax)
tax <- tax[,c(8,1:7)] #taxID
keep <- tax[tax$taxa_id %in% rownames(feature_table), ]
tax <- keep
rm(keep)

# Number of taxa except structural zeros
n_taxa = ifelse(is.null(struc_zero), nrow(feature_table), sum(apply(struc_zero, 1, sum) == 0))
# Cutoff values for declaring differentially abundant taxa
cut_off = c(0.9 * (n_taxa -1), 0.8 * (n_taxa -1), 0.7 * (n_taxa -1), 0.6 * (n_taxa -1))
names(cut_off) = c("detected_0.9", "detected_0.8", "detected_0.7", "detected_0.6")

# Annotation data
dat_ann = data.frame(x = min(res$fig$data$x), y = cut_off["detected_0.8"], label = "W[0.8]")

#test figure
fig = res$fig +  
   geom_hline(yintercept = cut_off["detected_0.8"], linetype = "dashed") + 
  geom_text(data = dat_ann, aes(x = x, y = y, label = label), 
              size = 4, vjust = -0.5, hjust = 0, color = "orange", parse = TRUE)
fig  
 
# # specialized plot
figdf <- as.data.frame(fig$data)
figdf <- cbind(figdf,tax, by = "taxa.id")
figdf <- figdf[,-6]
figdf$group <- as.factor(figdf$group)
levels(figdf$group) <- c("3 vs. 0 weeks", "6 vs. 0 weeks")

#reciprocal contrast
levels(figdf$group) <- c("0 vs. 3 weeks", "6 vs. 3 weeks")

# order genus
x = tapply(figdf$y, figdf$Genus, function(x) max(x))
x = sort(x, TRUE)
figdf$Genus = factor(as.character(figdf$Genus), levels=names(x))
figdf$col_genus <- figdf$Genus

#nitrate only
figdf$col_genus[figdf$col_genus != "Family_Helicobacteraceae" &
                  figdf$col_genus != "Family_Francisellaceae" &
                  figdf$col_genus != "Unclassified_ASV_6"] <- NA

#combined only, kept taxa sig at 0.7
figdf$col_genus[figdf$col_genus != "Aquarickettsia" &
                  figdf$col_genus != "Family_Francisellaceae" &
                  figdf$col_genus != "Spirochaeta_2" &
                  figdf$col_genus != "Unclassified_ASV_6"] <- NA


#phosphate only
figdf$col_genus[figdf$col_genus != "Aquarickettsia"] <- NA

#ammonium and no treatment
figdf$col_genus <- NA

#geno 50 nutrient only
figdf$col_genus[figdf$col_genus != "Unclassified_ASV_6" &
                  figdf$col_genus != "Aquarickettsia" &
                  figdf$col_genus != "Family_Francisellaceae" &
                  figdf$col_genus != "Spirochaeta_2" &
                  figdf$col_genus != "Family_Helicobacteraceae"] <- NA

levels(figdf$col_genus)
#add new factor
figdf$col_genus <- factor(figdf$col_genus, levels = c(levels(figdf$col_genus), "Other"))
# convert NAs to other
figdf$col_genus[is.na(figdf$col_genus)] = "Other"

g50_nut_colors = c('#BE0032', '#F3C300', '#831DCB', '#37A337','#0067A5', '#848482')
g50_nitrate_colors = c('#0067A5', '#37A337','#F3C300', '#848482')
g50_phosphate_colors = c('#BE0032', '#848482')
g50_ammonium_colors = c('#848482')
g50_combined_colors = c('#BE0032','#37A337', '#831DCB', '#F3C300', '#848482')

# Taxa above the dashed line are significant with the null‐hypothesis rejected 80% of the time.
ggfig <- ggplot(figdf, aes(x = x, y = y, color = col_genus)) +
  geom_vline(xintercept = 0, color = "grey") +
  geom_point(size = 2) +
  facet_grid(~group) +
  theme_bw() +
  ylab("W statistic") +
  xlim(-4, 2) +
  theme(axis.title.x=element_blank()) +
  theme(legend.position = "none") +
  scale_color_manual(name = "Genus", values = g50_ammonium_colors) +
  geom_hline(yintercept = cut_off["detected_0.8"], linetype = "dashed") +
  annotate("text", label = "W = 0.8", size = 3.5, x = -2, y = 23.5, color = "orange")
ggfig

ggfig2 <- ggplot(subset(figdf, group == "6 vs. 3 weeks"), aes(x = x, y = y, color = col_genus)) +
  geom_vline(xintercept = 0, color = "grey") +
  geom_point(size = 2) +
  theme_bw() +
  facet_grid(~group) +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  theme(axis.title.x=element_blank()) +
  xlim(-4, 2) +
  theme(legend.position = "none") +
  scale_color_manual(name = "Genus", values = g50_ammonium_colors) +
  geom_hline(yintercept = cut_off["detected_0.8"], linetype = "dashed") + 
  annotate("text", label = "W = 0.8", size = 3.5, x = -2, y = 23.5, color = "orange")
ggfig2

#ammonium panel
ammonium <- ggarrange(ggfig + theme(plot.margin = margin(r = 1)), ggfig2 + theme(plot.margin = margin(l = 1)), nrow = 1, align = "h", widths = c(2.1,1))
ammonium <- annotate_figure(ammonium,
                          bottom = text_grob("CLR Mean Difference", hjust = 1, color = "black", size = 10))
ammonium

#nitrate panel
nitrate <- ggarrange(ggfig + theme(plot.margin = margin(r = 1)), ggfig2 + theme(plot.margin = margin(l = 1)), nrow = 1, align = "h", widths = c(2.1,1))
nitrate <- annotate_figure(nitrate,
                            bottom = text_grob("CLR Mean Difference", hjust = 1, color = "black", size = 10))
nitrate

#phosphate panel
phosphate <- ggarrange(ggfig + theme(plot.margin = margin(r = 1)), ggfig2 + theme(plot.margin = margin(l = 1)), nrow = 1, align = "h", widths = c(2.1,1))
phosphate <- annotate_figure(phosphate,
                           bottom = text_grob("CLR Mean Difference", hjust = 1, color = "black", size = 10))
phosphate

#combined panel
combined <- ggarrange(ggfig + theme(plot.margin = margin(r = 1)), ggfig2 + theme(plot.margin = margin(l = 1)), nrow = 1, align = "h", widths = c(2.1,1))
combined <- annotate_figure(combined,
                             bottom = text_grob("CLR Mean Difference", hjust = 1, color = "black", size = 10))
combined

#merge em all
p <- ggarrange(ammonium, combined, nitrate, phosphate,
          common.legend = TRUE, legend = "bottom", labels = c("A", "B", "C", "D"))

ggsave(filename="~/Mote_nutrient_experiment/plots/ANCOM_50_nitrate.png", plot=nitrate, device="png",  width = 6.5, height = 3, dpi=500)
ggsave(filename="~/Mote_nutrient_experiment/plots/ANCOM_50_ammonium.svg", plot=merged2, device="svg", dpi=500)
ggsave(filename="~/Mote_nutrient_experiment/plots/ANCOM_50_phosphate.png", plot=phosphate, device="png",  width = 6.5, height = 3, dpi=500)


ggsave(filename="~/Mote_nutrient_experiment/plots/merged.svg", plot=p, device="svg",  width = 12, height = 6, dpi=500)
