library('phyloseq'); packageVersion('phyloseq')
library('vegan'); packageVersion('vegan')
library("usedist")
library("ggplot2")
library("microbiome")
library("multcompView")

setwd("~/Mote_nutrient_experiment/data")

pairwise.adonis.dm <- function(x,factors,stratum=NULL,p.adjust.m="fdr",perm=999){
  
  library(vegan)
  if(class(x) != "dist") stop("x must be a dissimilarity matrix (dist object)")
  co = as.matrix(combn(unique(factors),2))
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  
  
  
  for(elem in 1:ncol(co)){
    sub_inds <- factors %in% c(as.character(co[1,elem]),as.character(co[2,elem]))
    resp <- as.matrix(x)[sub_inds,sub_inds]
    ad = adonis(as.dist(resp) ~
                  
                  factors[sub_inds], strata=stratum[sub_inds], permutations=perm);
    
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
    
  }
  
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted)
  return(pairw.res)
}


load(file = "ps_clr_g50.RData") #renamed, clr transformed data
load(file = "clr_no_rick.RData") #renamed, clr transformed data with no aquarickettsia
sample_sums(clr)

mapfile = "Mapping-file-full-renamed.txt" #confirm correct mapping file associated
map = import_qiime_sample_data(mapfile)
sample_data(clr) <- map
sample_data(clr_no_rick) <- map
levels(clr@sam_data$Nutrient_no_level) <- c("Ammonium", "Nitrate", "Phosphate", "Combined", "No Treatment", "T0")


### Calculate distance matrices
# only calculating euclidean distance for clr transformed data because we have negative values, and with the clr transform we've moved our data into "real space"
clr_euc <- phyloseq::distance(clr, method = "euclidean")

### PERMANOVAs

# PERMANOVA's with Adonis -  Permutational Multivariate Analasis of Variance
# Make a data frame from the sample_data
sampledf <- data.frame(sample_data(clr))
# Adonis 
adonis(clr_euc ~ Nutrient_no_level*Level, data = sampledf) #sig, p = 0.001

pairwise.adonis.dm(clr_euc, sample_data(clr)$Nutrient_no_level, p.adjust.m = "fdr")
#comparisons against T0 significant, explore further. 

pairwise.adonis.dm(clr_euc, sample_data(clr)$Exposure_weeks, p.adjust.m = "fdr")
#all significant, interesting because non clr transformed data was NS for 3 vs 6

pairwise.adonis.dm(clr_euc, sample_data(clr)$T0_merge, p.adjust.m = "fdr")
#comparisons against T0 significant...

pairwise.adonis.dm(clr_euc, sample_data(clr)$Nutrient, p.adjust.m = "fdr")
#no comparison significant

pairwise.adonis.dm(clr_euc, sample_data(clr)$No_level_weeks, p.adjust.m = "fdr")
#T0 vs CTRL3 NS at p <0.01

pairwise.adonis.dm(clr_euc, sample_data(clr)$Level, p.adjust.m = "fdr")
#0 vs H and vs L sig, NS H vs L

### PERMDISPR
anova(betadisper(clr_euc, sampledf$Nutrient, bias.adjust = TRUE)) #NS
anova(betadisper(clr_euc, sampledf$Nutrient_no_level, bias.adjust = TRUE)) #significant
# Response: Distances
# Df Sum Sq Mean Sq F value    Pr(>F)    
# Groups      5 1779.7  355.95  11.548 2.225e-09 ***
#   Residuals 143 4407.7   30.82  

anova(betadisper(clr_euc, sampledf$Exposure_weeks, bias.adjust = TRUE)) #sig
# Response: Distances
# Df Sum Sq Mean Sq F value    Pr(>F)    
# Groups      2 1269.1  634.55  20.074 1.986e-08 ***
#   Residuals 146 4615.0   31.61                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

anova(betadisper(clr_euc, sampledf$No_level_weeks, bias.adjust = TRUE))
# Response: Distances
# Df Sum Sq Mean Sq F value    Pr(>F)    
# Groups     10 1937.9 193.789  5.7772 3.153e-07 ***
#   Residuals 138 4629.1  33.544                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

p.adjust(permutest(betadisper(clr_euc, sampledf$Exposure_weeks, 
                              bias.adjust = TRUE), 
                   pairwise=TRUE)$pairwise$permuted, method = 'fdr')
#3 to 6 NS

anova(betadisper(clr_euc, sampledf$Level, bias.adjust = TRUE))
p.adjust(permutest(betadisper(clr_euc, sampledf$Level, 
                              bias.adjust = TRUE), 
                   pairwise=TRUE)$pairwise$permuted, method = 'fdr')
#0 to high and low sig, high to low NS

p.adjust(permutest(betadisper(clr_euc, sampledf$No_level_weeks, 
                              bias.adjust = TRUE), 
                   pairwise=TRUE)$pairwise$permuted, method = 'fdr')

#PCA ordination
ord_clr <- phyloseq::ordinate(clr, "RDA") #RDA without constraints = PCA
phyloseq::plot_scree(ord_clr) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n") #screen plot plateaus quickly
head(ord_clr$CA$eig)
sapply(ord_clr$CA$eig[1:5], function(x) x / sum(ord_clr$CA$eig))    

clr@sam_data$Exposure_weeks <- as.factor(clr@sam_data$Exposure_weeks)
#plot
clr1 <- ord_clr$CA$eig[1] / sum(ord_clr$CA$eig)
clr2 <- ord_clr$CA$eig[2] / sum(ord_clr$CA$eig)
PCA <- phyloseq::plot_ordination(clr, ord_clr, type="samples", color="Nutrient_no_level", shape="Exposure_weeks") + 
  geom_point(size = 2) +
  coord_fixed(((clr2 / clr1))*4) +
  stat_ellipse(aes(group = Nutrient_no_level), linetype = 1) +
  ggtitle("Principal Component Analysis, All Treatments and Timepoints") +
  labs(color = "Treatment") +
  labs(shape = "Exposure Weeks") +
  scale_color_manual(values = c("#BE0032", "#009E73","#F3C300","#848482","#0067A5", "#000000"))
PCA
ggsave(filename="~/Mote_nutrient_experiment/plots/PCA_clr.png", plot=PCA, device="png", dpi=700)

adonis(clr_euc ~ Nutrient_no_level*Level, data = sampledf) 
dispr <- betadisper(clr_euc, sampledf$Nutrient_no_level, bias.adjust = TRUE)
dispr
#Average distance to median:
#  Ammonium  Combined      Ctrl   Nitrate Phosphate        T0 
#17.09     12.88     15.19     18.40     13.56     22.07 
plot(dispr, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")
boxplot(dispr, main = "", xlab = "") #quick dispersion plot
permutest(dispr)
#dispersion highest in T0s, way more samples

# extract distance to centroid
clr_euc <- phyloseq::distance(clr, method = "euclidean")
sampledf <- data.frame(sample_data(clr))
disp <- betadisper(clr_euc, sampledf$No_level_weeks, bias.adjust = TRUE)
dispd <- as.data.frame(disp$distances)
dispd <- cbind(dispd, sample_data(clr))
colnames(dispd)[1] <- "distance"

anova(betadisper(clr_euc, sampledf$No_level_weeks, bias.adjust = TRUE)) #sig
anova(betadisper(clr_euc, sampledf$Exposure_weeks, bias.adjust = TRUE)) #also sig
dispersion_stats <- p.adjust(permutest(betadisper(clr_euc, sampledf$No_level_weeks_noT0, 
                                                  bias.adjust = TRUE), 
                                       pairwise=TRUE)$pairwise$permuted, method = 'fdr')
write.table(dispersion_stats, file = "dispersion_stats_clr.txt", sep = "\t")
#manually turned this into a triangular matrix
dispersion_table <- read.table("dispersion_matrix_clr.txt", header = TRUE, sep = "\t", row.names = 1)
myletters<-multcompLetters(dispersion_table,compare="<=",threshold=0.05,Letters=letters)

dispd$Exposure_weeks <- as.factor(dispd$Exposure_weeks)
levels(dispd$Exposure_weeks) <- c("0","3", "6")
levels(dispd$Nutrient_no_level) <- c("Ammonium", "Nitrate", "Phosphate", "Combined", "No Treatment", "T0")
myCols <- c("#EF476F", "#06D6A0", "#FFC107")

dispersion_plot <- ggplot(dispd, aes(x=Exposure_weeks, y=distance)) +
  facet_wrap(~No_level_no_T0)+
  geom_boxplot(outlier.shape = NA, color = "gray35") +
  geom_point(aes(color = Exposure_weeks), 
             position = position_jitter(width = .25, height = 0)) +
  ylab("Distance-to-centroid") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = c(0.85, 0.3)) +
  labs(color = "Exposure Weeks") +
  scale_colour_manual(values = myCols) +
  stat_summary(geom = 'text', label = c("a","bcd","abc", "ab", "d", "c", "a", "abcd", "abc",  "ab", "abcd", "abcd","ab", "bcd", "cd"), fun.y = max, vjust = -0.5, size = 4) +
  ggtitle("Dispersion Over Time, All Treatments") +
  scale_y_continuous(expand = expand_scale(mult = c(.1)))
dispersion_plot
ggsave(filename="~/Mote_nutrient_experiment/plots/Dispersion_all_clr.png", plot=dispersion_plot, height = 6, width = 6, device="png", dpi=500)

#Dispersion with Aquarickettsia removed
load(file = "clr_no_rick.RData")
sample_data(clr_no_rick) <- map
clr_no_rick <- subset_samples(clr_no_rick, SampleID != "LP-50-2-T0")
clr_euc <- phyloseq::distance(clr_no_rick, method = "euclidean")
sampledf <- data.frame(sample_data(clr_no_rick))
disp <- betadisper(clr_euc, sampledf$No_level_weeks, bias.adjust = TRUE)
dispd <- as.data.frame(disp$distances)
dispd <- cbind(dispd, sample_data(clr_no_rick))
colnames(dispd)[1] <- "distance"
boxplot(disp, main = "", xlab = "") 

anova(betadisper(clr_euc, sampledf$No_level_weeks, bias.adjust = TRUE))
dispersion_stats <- p.adjust(permutest(betadisper(clr_euc, sampledf$No_level_weeks_noT0, 
                                                  bias.adjust = TRUE), 
                                       pairwise=TRUE)$pairwise$permuted, method = 'fdr')
write.table(dispersion_stats, file = "dispersion_stats_no_rick_clr.txt", sep = "\t")
#manually turned this into a triangular matrix
dispersion_table <- read.table("dispersion_matrix_no_rick_clr.txt", header = TRUE, sep = "\t", row.names = 1)
myletters<-multcompLetters(dispersion_table,compare="<=",threshold=0.05,Letters=letters)

dispd$Exposure_weeks <- as.factor(dispd$Exposure_weeks)
levels(dispd$Exposure_weeks) <- c("0","3", "6")
myCols <- c("#EF476F", "#06D6A0", "#FFC107")
no_rick_plot <- ggplot(dispd, aes(x=Exposure_weeks, y=distance)) +
  facet_wrap(~No_level_no_T0)+
  geom_boxplot(outlier.shape = NA, color = "gray35") +
  geom_point(aes(color = Exposure_weeks), 
             position = position_jitter(width = .25, height = 0)) +
  ylab("Distance-to-centroid") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = c(0.85, 0.3)) +
  labs(color = "Exposure Weeks") +
  scale_colour_manual(values = myCols) +
  stat_summary(geom = 'text', label = c("ab","ab","ab","ab","b","ab","ab","ab","ab","a","ab","ab","ab","ab","b"), fun.y = max, vjust = -0.5, size = 4) +
  ggtitle("Dispersion Over Time, Aquarickettsia Removed") +
  scale_y_continuous(expand = expand_scale(mult = c(.1)))
no_rick_plot
ggsave(filename="~/Mote_nutrient_experiment/plots/Dispersion_rick_removed.png", plot=no_rick_plot, height = 6, width = 6, device="png", dpi=500)
