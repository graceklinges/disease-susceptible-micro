library('phyloseq'); packageVersion('phyloseq')
library('vegan'); packageVersion('vegan')
library("usedist")
library("ggplot2")
library("microbiome")
library("multcompView")

setwd("~/Mote_nutrient_experiment/data")

pairwise.adonis.dm <- function(x,factors,stratum=NULL,p.adjust.m="bonferroni",perm=999){
  
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


load(file = "ps_rare_g50.RData") #renamed, rarefied data
ps <- ps_rarefied
rm(ps_rarefied)
sample_sums(ps)

mapfile = "Mapping-file-full-renamed.txt" #confirm correct mapping file associated
map = import_qiime_sample_data(mapfile)
sample_data(ps) <- map

ps <- subset_samples(ps, SampleID != "LP-50-2-T0") #outlier

### Calculate distance matrices
# Calculate all four distance matrices: weighted unifrac, unweighted unifrac, bray curtis, binary jaccard
ps_bc <- phyloseq::distance(ps, method = "bray")
ps_bj <- distance(ps, method = "jaccard", binary =TRUE)

### PERMANOVAs

# PERMANOVA's with Adonis -  Permutational Multivariate Analysis of Variance
# Make a data frame from the sample_data
sampledf <- data.frame(sample_data(ps))
# Adonis 
adonis(ps_bc ~ Nutrient_no_level*Level, data = sampledf)

pairwise.adonis.dm(ps_bc, sample_data(ps)$NW_T0_merge, p.adjust.m = "fdr")
#too many comparisons
pairwise.adonis.dm(ps_bc, sample_data(ps)$Nutrient_no_level, p.adjust.m = "fdr")
#some sig comparisons, explore further
pairwise.adonis.dm(ps_bj, sample_data(ps)$Nutrient_no_level, p.adjust.m = "fdr")
#fewer sig comparisons than BC

pairwise.adonis.dm(ps_bc, sample_data(ps)$Exposure_weeks, p.adjust.m = "fdr")
#0 vs 3 and 0 vs 6 significant, 3 vs 6 ns

pairwise.adonis.dm(ps_bc, sample_data(ps)$T0_merge, p.adjust.m = "fdr")
#comparisons against T0 significant...

pairwise.adonis.dm(ps_bc, sample_data(ps)$Nutrient, p.adjust.m = "fdr")
#no comparison significant

### PERMDISPs
anova(betadisper(ps_bc, sampledf$Nutrient, bias.adjust = TRUE)) #NS
anova(betadisper(ps_bc, sampledf$Nutrient_no_level, bias.adjust = TRUE)) #significant
#Response: Distances
#Df  Sum Sq  Mean Sq F value    Pr(>F)    
#Groups      5 0.59961 0.119922    13.5 8.473e-11 ***
#  Residuals 144 1.27915 0.008883                      
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
anova(betadisper(ps_bc, sampledf$Exposure_weeks, bias.adjust = TRUE)) #sig
#Response: Distances
#Df Sum Sq  Mean Sq F value    Pr(>F)    
#Groups      2 0.4504 0.225202  23.406 1.498e-09 ***
#  Residuals 147 1.4144 0.009622                      
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

p.adjust(permutest(betadisper(ps_bc, sampledf$Exposure_weeks, 
                              bias.adjust = TRUE), 
                   pairwise=TRUE)$pairwise$permuted, method = 'fdr')
#same thing as before, 3 to 6 NS

anova(betadisper(ps_bc, sampledf$Level, bias.adjust = TRUE))
p.adjust(permutest(betadisper(ps_bc, sampledf$Level, 
                              bias.adjust = TRUE), 
                   pairwise=TRUE)$pairwise$permuted, method = 'fdr')
#0 to high and low sig, high to low NS

#PCoA plot
ps_T0 <- subset_samples(ps, Nutrient_no_level == "T0")
ps_not_t0 <- subset_samples(ps, Nutrient_no_level != "T0")
ps_bc <- phyloseq::distance(ps, method = "bray")
ps_bc <- phyloseq::distance(ps_T0, method = "bray")
ps_bc <- phyloseq::distance(ps_not_t0, method = "bray")

logt  = transform_sample_counts(ps, function(x) log(1 + x) )
out.pcoa.logt <- ordinate(logt, method = "PCoA", distance = "bray")
evals <- out.pcoa.logt$values$Eigenvalues

plot_ordination(logt, out.pcoa.logt, type = "samples", 
                color = "Nutrient_no_level") + labs(col = "Nutrient") +
  coord_fixed(sqrt(evals[2] / evals[1]))

# Prepare data for plotting
meta <- as.data.frame(sample_data(ps_not_t0))
table <- as.data.frame(otu_table(ps_not_t0))
head(meta)

ctrl <- rownames(meta[which(meta[,14] =="Ctrl"),])
T0 <- rownames(meta[which(meta[,14] =="T0"),])
nut <- rownames(meta[which(meta[,14] !="Ctrl" & meta[,14] !="T0"),])
ammonium <- rownames(meta[which(meta[,14] =="Ammonium"),])
nitrate <- rownames(meta[which(meta[,14] =="Nitrate"),])
phosphate <- rownames(meta[which(meta[,14] =="Phosphate"),])
combined <- rownames(meta[which(meta[,14] =="Combined"),])
three <- rownames(meta[which(meta[,6] =="3"),])
six <- rownames(meta[which(meta[,6] =="6"),])
Nutrient <- meta$Nutrient_no_level
Weeks <- meta$Exposure_weeks


# params for plotting
dims <- c(1,2)
ellp.kind <- "ehull"

# ordinate with Bray Curtis
object <- metaMDSiter(ps_bc, k=2, trymax = 1000, maxit = 1000, autotransform=FALSE) #wow finally reached solution at run 977 stress 0.1109059 
save(object, file = "ord_bc_50.RData")

mds.fig <- ordiplot(object, display = "sites", type = "none", choices = dims)
plot(1, type="n", main= "NMDS All Timepoints", xlab="MDS1", ylab="MDS2", xlim=c(-4.5,1.4), ylim=c(-3,3))
points(mds.fig, "sites", pch = 18, col = "#1E88E5", select = T0)
points(mds.fig, "sites", pch = 17, col = "#D81B60", select = nut)
points(mds.fig, "sites", pch = 19, col = "#009E73", select = ctrl)
ordiellipse(object, Nutrient, conf = 0.95, label = FALSE, choices = dims, kind = ellp.kind, col = "#009E73", lwd = 2, show.groups = "Ctrl")
ordiellipse(object, Nutrient, conf = 0.95, label = FALSE, choices = dims, kind = ellp.kind, col = "#1E88E5", lwd = 2, show.groups = "T0")
ordiellipse(object, Nutrient, conf = 0.95, label = FALSE, choices = dims, kind = ellp.kind, col = "#D81B60", lwd = 2, show.groups = "Phosphate")
ordiellipse(object, Nutrient, conf = 0.95, label = FALSE, choices = dims, kind = ellp.kind, col = "#D81B60", lwd = 2, show.groups = "Ammonium")
ordiellipse(object, Nutrient, conf = 0.95, label = FALSE, choices = dims, kind = ellp.kind, col = "#D81B60", lwd = 2, show.groups = "Nitrate")
ordiellipse(object, Nutrient, conf = 0.95, label = FALSE, choices = dims, kind = ellp.kind, col = "#D81B60", lwd = 2, show.groups = "Combined")
legend("bottomleft", legend = c("T0", "Control", "Nutrient"), pch = c(18, 19, 17), col = c("#1E88E5","#009E73","#D81B60"))
text(0.7,2.8, labels = c("stress=0.111"))

mds.fig <- ordiplot(object, display = "sites", type = "points", choices = dims)
plot(1, type="n", main= "NMDS 3 and 6 Weeks", xlab="MDS1", ylab="MDS2", xlim=c(-3.75,3), ylim=c(-3,2.2))
points(mds.fig, "sites", pch = 20, col = "#BE0032", select = ammonium)
points(mds.fig, "sites", pch = 19, col = "#F3C300", select = nitrate)
points(mds.fig, "sites", pch = 18, col = "#0067A5", select = phosphate)
points(mds.fig, "sites", pch = 15, col = "#009E73", select = combined) #square
points(mds.fig, "sites", pch = 17, col = "#848482", select = ctrl) #triangle
ordiellipse(object, Nutrient, conf = 0.95, label = FALSE, choices = dims, kind = ellp.kind, col = "#848482", lwd = 2, show.groups = "Ctrl")
ordiellipse(object, Nutrient, conf = 0.95, label = FALSE, choices = dims, kind = ellp.kind, col = "#0067A5", lwd = 2, show.groups = "Phosphate")
ordiellipse(object, Nutrient, conf = 0.95, label = FALSE, choices = dims, kind = ellp.kind, col = "#BE0032", lwd = 2, show.groups = "Ammonium")
ordiellipse(object, Nutrient, conf = 0.95, label = FALSE, choices = dims, kind = ellp.kind, col = "#F3C300", lwd = 2, show.groups = "Nitrate")
ordiellipse(object, Nutrient, conf = 0.95, label = FALSE, choices = dims, kind = ellp.kind, col = "#009E73", lwd = 2, show.groups = "Combined")
legend("bottomleft", legend = c("Control", "Ammonium", "Nitrate", "Phosphate", "Combined"), pch = c(17, 20, 19, 18, 15), col = c("#848482","#BE0032","#F3C300", "#0067A5", "#009E73"))
text(2.3,2, labels = c("stress=0.08"))

mds.fig <- ordiplot(object, display = "sites", type = "points", choices = dims)
plot(1, type="n", main= "NMDS 3 and 6 Weeks", xlab="MDS1", ylab="MDS2", xlim=c(-3.75,3), ylim=c(-3,2.2))
points(mds.fig, "sites", pch = 20, col = "#BE0032", select = three)
points(mds.fig, "sites", pch = 15, col = "#009E73", select = six)
ordiellipse(object, Weeks, conf = 0.95, label = FALSE, choices = dims, kind = ellp.kind, col = "#BE0032", lwd = 2, show.groups = "3")
ordiellipse(object, Weeks, conf = 0.95, label = FALSE, choices = dims, kind = ellp.kind, col = "#009E73", lwd = 2, show.groups = "6")
legend("bottomleft", legend = c("Three Weeks", "Six Weeks"), pch = c(20, 15), col = c("#BE0032", "#009E73"))
text(2.3,2, labels = c("stress=0.08"))

# extract distance to centroid
load(file = "ps_rare_g50.RData") #renamed, rarefied data
ps <- ps_rarefied
ps <- subset_samples(ps, SampleID != "LP-50-2-T0")
ps_ctrl <- subset_samples(ps, Nutrient == "CTRL")
ps_nutrient <- subset_samples(ps, Nutrient != "CTRL")
ps_bc <- phyloseq::distance(ps, method = "bray")
ps_ctrl_bc <- phyloseq::distance(ps_ctrl, method = "bray")
ps_nutrient_bc <- phyloseq::distance(ps_nutrient, method = "bray")

sampledf <- data.frame(sample_data(ps))
disp <- betadisper(ps_bc, sampledf$Exposure_weeks, bias.adjust = TRUE)
dispd <- as.data.frame(disp$distances)
dispd <- cbind(dispd, sample_data(ps))
colnames(dispd)[1] <- "distance"

anova(betadisper(ps_bc, sampledf$No_level_weeks, bias.adjust = TRUE))
dispersion_stats <- p.adjust(permutest(betadisper(ps_bc, sampledf$No_level_weeks, 
                              bias.adjust = TRUE), 
                   pairwise=TRUE)$pairwise$permuted, method = 'fdr')
write.table(dispersion_stats, file = "dispersion_stats.txt", sep = "\t")
#manually turned this into a triangular matrix
dispersion_table <- read.table("dispersion_matrix.txt", header = TRUE, sep = "\t", row.names = 1)
myletters<-multcompLetters(dispersion_table,compare="<=",threshold=0.05,Letters=letters)

dispd$Exposure_weeks <- as.factor(dispd$Exposure_weeks)
levels(dispd$Exposure_weeks) <- c("0","3", "6")
myCols <- c("#EF476F", "#06D6A0", "#FFC107")
a <- ggplot(dispd, aes(x=Exposure_weeks, y=distance)) +
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
  scale_colour_manual(values = myCols) +
  stat_summary(geom = 'text', label = c("abcd","abe","a", "cd", "ef", "f", "bcd", "abcdef", "abcdef", "abcd", "abce", "abce", "d", "f", "f"), fun.y = max, vjust = -0.5, size = 4) +
  ggtitle("Dispersion Over Time, All Treatments") +
  scale_y_continuous(expand = expand_scale(mult = c(.1)))
a
ggsave(filename="Dispersion_all.png", plot=a, device="png", dpi=500)

#Dispersion with Aquarickettsia removed
subset <- load(file = "ps_rare_g50.RData")
ps_bc <- phyloseq::distance(subset, method = "bray")
sampledf <- data.frame(sample_data(subset))
disp <- betadisper(ps_bc, sampledf$Exposure_weeks, bias.adjust = TRUE)
dispd <- as.data.frame(disp$distances)
dispd <- cbind(dispd, sample_data(subset))
colnames(dispd)[1] <- "distance"

anova(betadisper(ps_bc, sampledf$No_level_weeks, bias.adjust = TRUE))
dispersion_stats <- p.adjust(permutest(betadisper(ps_bc, sampledf$No_level_weeks, 
                                                  bias.adjust = TRUE), 
                                       pairwise=TRUE)$pairwise$permuted, method = 'fdr')
write.table(dispersion_stats, file = "dispersion_stats_no_rick.txt", sep = "\t")
#there were no significant comparisons

dispd$Exposure_weeks <- as.factor(dispd$Exposure_weeks)
levels(dispd$Exposure_weeks) <- c("0","3", "6")
myCols <- c("#EF476F", "#06D6A0", "#FFC107")
a <- ggplot(dispd, aes(x=Exposure_weeks, y=distance)) +
  facet_wrap(~No_level_no_T0)+
  geom_boxplot(outlier.shape = NA, color = "gray35") +
  geom_point(aes(color = Exposure_weeks), 
             position = position_jitter(width = .25, height = 0)) +
  ylab("Distance-to-centroid") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none") +
  scale_colour_manual(values = myCols) +
  stat_summary(geom = 'text', label = c("a","a","a","a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a"), fun.y = max, vjust = -1, size = 4) +
  ggtitle("Dispersion Over Time, Aquarickettsia Removed") +
  scale_y_continuous(expand = expand_scale(mult = c(.1)))
a
ggsave(filename="Dispersion_rick_removed.png", plot=a, device="png", dpi=500)