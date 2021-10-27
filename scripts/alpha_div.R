library(phyloseq)
library(multcomp)
library(picante)
library(ggplot2 )
library(tidyverse)
library(plyr)
library(dplyr)

setwd("~/Mote_nutrient_experiment/data")

#load for figure 2A
load(file = "ps_rare_g50.RData") #renamed, unfiltered, rarefied ps object

#load for supp fig. 2
#load(file = "ps_no_rick_rare.RData") #no aquarickettsia dataset, rarefied to 351

ps <- ps_rarefied
rm(ps_rarefied)
sample_sums(ps) #should be 5098, or 351 for no rick
ps = subset_samples(ps, SampleID != "Negative")
mapfile = "Mapping-file-full-renamed.txt"
map = import_qiime_sample_data(mapfile)
sample_data(ps) <- map

#functions
asinTransform <- function(p) { asin(sqrt(p)) }

# Estimate faith's phylogenetic diversity 
estimate_pd <- function(phylo){
  
  # Error if input is not of class phylo
  if(class(phylo) != "phyloseq"){
    stop("Input file is not of class 'phyloseq'.")
  }
  
  # Error if no class phy_tree
  if(!(.hasSlot(phylo, "phy_tree"))){
    stop("Could not find tree slot in phylo object.")
  }
  
  # Transpose if needed
  # Adapted from phyloseq/vegan import
  OTU <- phyloseq::otu_table(phylo)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  
  # Get matrix version of OTU table
  otutable <- as(OTU, "matrix")
  
  # Get phylogenetic tree from pyloseq object
  tree <- phyloseq::phy_tree(phylo)
  
  # Print status message
  message("Calculating Faiths PD-index...")
  
  # If object is greater than 10mb, then print status message
  if(object.size(otutable) > 10000000){
    message("This is a large object, it may take awhile...")
  }
  
  # Calculate Faith's PD-index
  pdtable <- picante::pd(otutable, tree, include.root = F)
  
  # Return data frame of results
  return(pdtable)
}

# Calculate standard error
sderr <- function(x) {sd(x)/sqrt(length(x))}
#+++++++++++++++++++++++++
# Function to calculate the mean and the standard deviation
# for each group
#+++++++++++++++++++++++++
# data : a data frame
# varname : the name of a column containing the variable
#to be summariezed
# groupnames : vector of column names to be used as
# grouping variables
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sderr(x[[col]]), na.rm=TRUE)
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

normality.plots <- function(x) {
  par(mfrow=c(2,2))
  hist(residuals(x), main = "Histogram", xlab = "Values")
  boxplot(residuals(x), main = "Boxplot", ylab = "Values")
  qqnorm(residuals(x))
  qqline(residuals(x))
  plot(density(residuals(x)), main = "Kernel Density Estimate")
}

#Initialize matrices to store alpha diversity estimates
nsamp = nsamples(ps)
observed <- matrix(nrow = nsamp)
row.names(observed) <- sample_names(ps)
simpson <- matrix(nrow =nsamp)
row.names(simpson) <- sample_names(ps)
shannon <- matrix(nrow =nsamp)
row.names(shannon) <- sample_names(ps)

#Calculate statistics
# Options for measures = ("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")
# Calculate observed
obs <- as.numeric(as.matrix(estimate_richness(ps, measures = "Observed")))
observed[ ,] <- obs
colnames(observed) [1] <- "observed"
# Calculate simpson
simp <- as.numeric(as.matrix(estimate_richness(ps, measures = "Simpson")))
simpson[ ,] <- simp
colnames(simpson) [1] <- "simpson"
# Calculate shannon
shan <- as.numeric(as.matrix(estimate_richness(ps, measures = "Shannon")))
shannon[ ,] <- shan
colnames(shannon) [1] <- "shannon"

#Combine our estimates for observed, simpson, and Shannon indices into one dataframe
alpha <- cbind(observed,simpson,shannon)
head(alpha)
# Add the sample metadata into this dataframe
s <- data.frame(sample_data(ps))
s$Exposure_weeks <- as.factor(s$Exposure_weeks)
levels(s$Exposure_weeks) <- c("0","3", "6")
alphadiv <- cbind(alpha, s)
head(alphadiv)
alphadiv <- alphadiv[,-5]
head(alphadiv)

write.csv(alphadiv, file = "./alphadiv_5098_rare_geno50.csv")
#write.csv(alphadiv, file = "./alphadiv_351_rare_geno50.csv")

# # Observed
observed1 <- data_summary(alphadiv, varname = "observed", groupnames = c("Exposure_weeks"))
observed2 <- data_summary(alphadiv, varname = "observed", groupnames = c("NW_T0_merge"))
observed3 <- data_summary(alphadiv, varname = "observed", groupnames = c("Nutrient_no_level"))
observed4 <- data_summary(alphadiv, varname = "observed", groupnames = c("No_level_weeks_noT0"))

# # Simpson
simp1 <- data_summary(alphadiv, varname = "simpson", groupnames = c("Exposure_weeks"))
simp2 <- data_summary(alphadiv, varname = "simpson", groupnames = c("NW_T0_merge"))
simp3 <- data_summary(alphadiv, varname = "simpson", groupnames = c("Nutrient_no_level"))
simp4 <- data_summary(alphadiv, varname = "simpson", groupnames = c("No_level_weeks_noT0"))

# # Shannon
shan1 <- data_summary(alphadiv, varname = "shannon", groupnames = c("Exposure_weeks"))
shan2 <- data_summary(alphadiv, varname = "shannon", groupnames = c("NW_T0_merge"))
shan3 <- data_summary(alphadiv, varname = "shannon", groupnames = c("Nutrient_no_level"))
shan4 <- data_summary(alphadiv, varname = "shannon", groupnames = c("No_level_weeks_noT0"))

# #Kruskal Wallis -- observed
# kruskal.test(observed ~ Nutrient_geno, data = alphadiv)
# kruskal.test(observed ~ Genotype, data = alphadiv)
# kruskal.test(observed ~ Nutrient, data = alphadiv)
# pairwise.wilcox.test(alphadiv$observed, alphadiv$Nutrient_geno, p.adjust.method = 'fdr')
# 

# #Kruskal Wallis -- Simpson's 
kruskal.test(simpson ~ Exposure_weeks, data = alphadiv) #sig with and without rick
kruskal.test(simpson ~ NW_T0_merge, data = alphadiv) #sig
kruskal.test(simpson ~ Nutrient_no_level, data = alphadiv) #sig
kruskal.test(simpson ~ No_level_weeks_noT0, data = alphadiv) #sig
kruskal.test(simpson ~ Nutrient, data = alphadiv) #NS
kruskal.test(simpson ~ No_level_no_T0, data = alphadiv) #NS

# 0 vs 3 and 0 vs 6 sig, 3 vs 6 NS
pairwise.wilcox.test(alphadiv$simpson, alphadiv$Exposure_weeks, p.adjust.method = 'fdr')
pairwise.wilcox.test(alphadiv$shannon, alphadiv$Exposure_weeks, p.adjust.method = 'fdr')

#All treatments significantly different from T0, supports NMDS plots
KW_simpsons_Nutrient_no_level <- pairwise.wilcox.test(alphadiv$simpson, alphadiv$Nutrient_no_level, p.adjust.method = 'fdr')

#most NS except comparisons to T0
KW_simpsons_NW_T0_merge <- pairwise.wilcox.test(alphadiv$simpson, alphadiv$NW_T0_merge, p.adjust.method = 'fdr')

#for plot
KW_simpsons_No_level_weeks_noT0 <- pairwise.wilcox.test(alphadiv$simpson, alphadiv$No_level_weeks_noT0, p.adjust.method = 'fdr')

#All treatments significantly different from T0, supports NMDS plots
KW_simpsons_T0_merge <- pairwise.wilcox.test(alphadiv$simpson, alphadiv$T0_merge, p.adjust.method = 'fdr')

write.csv(KW_simpsons_Nutrient_no_level[["p.value"]], file = "KW_simpsons_Nutrient_no_level.csv")
write.csv(KW_simpsons_NW_T0_merge[["p.value"]], file = "KW_simpsons_NW_T0_merge.csv")
write.csv(KW_simpsons_T0_merge[["p.value"]], file = "KW_simpsons_T0_merge.csv")
write.csv(KW_simpsons_No_level_weeks_noT0[["p.value"]], file = "KW_simpsons_No_level_weeks_noT0.csv")


write.csv(KW_simpsons_Nutrient_no_level[["p.value"]], file = "KW_simpsons_Nutrient_no_level_rick_rem.csv")
write.csv(KW_simpsons_NW_T0_merge[["p.value"]], file = "KW_simpsons_NW_T0_merge_rick_rem.csv")
write.csv(KW_simpsons_T0_merge[["p.value"]], file = "KW_simpsons_T0_merge_rick_rem.csv")
write.csv(KW_simpsons_No_level_weeks_noT0[["p.value"]], file = "KW_simpsons_No_level_weeks_noT0_rick_rem.csv")
# 
# #Kruskal Wallis -- Shannon's
kruskal.test(shannon ~ Exposure_weeks, data = alphadiv) #sig
kruskal.test(shannon ~ NW_T0_merge, data = alphadiv) #sig
kruskal.test(shannon ~ Nutrient_no_level, data = alphadiv) #sig
kruskal.test(shannon ~ No_level_weeks_noT0, data = alphadiv) #sig
kruskal.test(shannon ~ Nutrient, data = alphadiv) #NS
kruskal.test(shannon ~ No_level_no_T0, data = alphadiv) #NS

#plots
myCols <- c("#EF476F", "#06D6A0", "#FFC107")

#goes with stats from KW_simpsons_No_level_weeks_noT0, use alpha_div_sig_letters to get letters
A <- ggplot(alphadiv, aes(x=Exposure_weeks, y=simpson)) +
  facet_wrap(~No_level_no_T0)+
  geom_boxplot(outlier.shape = NA, color = "gray35") +
  geom_point(aes(color = Exposure_weeks), 
             position = position_jitter(width = .25, height = 0)) +
  ylab("Simpson's Diversity Index") +
  theme_bw() +
  theme(aspect.ratio = 1.8) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = c(0.9,0.3),
        legend.title = element_text(color = "black", size = 12),
        legend.text = element_text(color = "black", size = 12)) +
  scale_colour_manual(values = myCols) +
  stat_summary(geom = 'text', label = c("abc","adef","abde","c", "f", "df", "c", "de", "abe", "bc", "df", "abdef","c", "f", "f"), fun = max, vjust = -1, size = 3.5) + #stats with rick
    ylim(0,1) +
  ggtitle("Simpson's Diversity by Nutrient") +
theme(plot.title = element_text(hjust = 0.5))
A
#merged with PCA in Illustrator
ggsave(filename="~/Mote_nutrient_experiment/plots/Nutrient_No_level_weeks_noT0_simpsons.svg", plot=A, device="svg", dpi=500)

#load 
B <- ggplot(alphadiv, aes(x=Exposure_weeks, y=simpson)) +
  facet_wrap(~No_level_no_T0)+
  geom_boxplot(outlier.shape = NA, color = "gray35") +
  geom_point(aes(color = Exposure_weeks), 
             position = position_jitter(width = .25, height = 0)) +
  ylab("Simpson's Diversity Index") +
  theme_bw() +
  theme(aspect.ratio = 1.8) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = c(0.9,0.3),
        legend.title = element_text(color = "black", size = 12),
        legend.text = element_text(color = "black", size = 12)) +
  scale_colour_manual(values = myCols) +
  stat_summary(geom = 'text', label = c("abc","abcd","abc","a", "d", "bd", "abc", "bcd", "abcd", "a", "abcd", "abcd","ac", "bd", "abcd"), fun = max, vjust = -1, size = 3.5) + #stats without rick
  ylim(0,1) +
  ggtitle("Simpson's Diversity by Nutrient, Aquarickettsia Removed") +
  theme(plot.title = element_text(hjust = 0.5))
B
ggsave(filename="~/Mote_nutrient_experiment/plots/Nutrient_No_level_weeks_noT0_simpsons_rick_rem.png", plot=B, device="png", dpi=500)


