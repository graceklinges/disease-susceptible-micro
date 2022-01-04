library(corncob)
library(phyloseq)
library(magrittr)
library(ggplot2)

setwd("~/Mote_nutrient_experiment/data")

# load unrarefied, renamed data for differential abundance analysis
load(file = "ps_rename_g50.RData")
sample_sums(ps)
#remove negatives
ps = subset_samples(ps, SampleID != "Negative")

mapfile = "Mapping-file-full-renamed.txt"
map = import_qiime_sample_data(mapfile)
sample_data(ps) <- map

#reordering levels
ps@sam_data$Nutrient_no_level <- as.factor(ps@sam_data$Nutrient_no_level)
ps@sam_data$Nutrient_no_level <- factor(ps@sam_data$Nutrient_no_level, levels = c("T0", "No Treatment", "Ammonium", "Combined", "Nitrate", "Phosphate")) 
ps@sam_data$No_level_weeks <- as.factor(ps@sam_data$No_level_weeks)
ps@sam_data$No_level_weeks <- factor(ps@sam_data$No_level_weeks, levels = c("T0", "NT3", "NT6", "A3", "A6", "N3", "N6", "C3", "C6", "P3", "P6")) 
ps@sam_data$No_level_weeks

#turn exposure weeks into categorical variable
ps@sam_data$Exposure_weeks <- as.character(ps@sam_data$Exposure_weeks)

otu_table(ps)[1:3, 1:3]
sample_data(ps)[1:3, ]

ps_corn <- tax_glom(ps, "Genus")
ntaxa(ps_corn) #576

#fitting a simple model. ASV1 is Aquarickettsia
corncob <- bbdml(formula = ASV_1 ~ 1,
                  phi.formula = ~ 1,
                  data = ps_corn)

# #plot data with model fit on rel abundance scale
plot(corncob, B = 50)
#The points represent the relative abundances. The bars represent the 95% prediction intervals for the
#observed relative abundance by sample. The parameter B determines the number of bootstrap simulations
#used to approximate the prediction intervals. For purposes of this tutorial, we use a small value B = 50 for
#computational purposes, but recommend a higher setting for more accurate intervals, such as the default B =1000.

#Now let’s look at the same plot, but on the counts scale with 95% prediction intervals 
#(since counts is not a parameter). To do this, we add the option total = TRUE to our plotting code.
# plot(corncob, total = TRUE, B = 50)
# plot(corncob, total = TRUE, color = "Nutrient_no_level", B = 50)
# plot(corncob, color = "Nutrient_no_level", B = 50)

#adding covariates. Modeling the expected relative abundance and the variability of the counts with nutrient as a covariate.
#modeling Aquarickettsia
corncob_nl <- bbdml(formula = ASV_1 ~ Nutrient_no_level,
                    phi.formula = ~ Nutrient_no_level,
                    data = ps_corn)

#test with fewer bootstraps
plot(corncob_nl, color = "Nutrient_no_level", total = TRUE, B = 50)
plot <- plot(corncob_nl, color = "Nutrient_no_level", B = 50, sample_names = FALSE)
plot

#figure 5
plot <- plot(corncob_nl, color = "Nutrient_no_level", B = 1000, sample_names = FALSE)
plot
#this model seems to fit the data much better

#I manually adjusted the color in Illustrator and added bars for mean and std dev
ggsave(filename="~/Mote_nutrient_experiment/plots/corncob.svg", plot=plot, width=10, height=8)

#now let's compare the models with a likelihood ratio test
#test null hypothesis that likelihood of model with covariates equals likelihood of model without covariates
lrtest(mod_null = corncob, mod = corncob_nl)
# p value was 1.49774e-21, super low p value means that these is a statistically significant difference in the likelihood of the models

summary(corncob_nl)
corncob_stats <- summary(corncob_nl)
write.table(corncob_stats$coefficients, file = "corncob_stats.txt", sep = "\t")

#if we run this without changing the levels (compare against Ammonium) we can see that the Nutrient_no_levelT0 abundance coefficient is strongly negative and statistically significant
#this suggests that ASV1 is differentially abundant across nutrient treatment and samples from T0 have a lower relative abundance

#with levels changed all treatments have a positive coefficient compared with T0, but phosphate is the highest
#all treatments had lower dispersion than T0, only phosphate and combined significantly decreased dispersion

#now let's look at other taxa
#keeping same model covariates as before, but removing response term
set.seed(1) #here controlling for the effect of Nutrient on dispersion
da_analysis <- differentialTest(formula = ~ Nutrient_no_level,
                                phi.formula = ~ Nutrient_no_level,
                                formula_null = ~ 1,
                                phi.formula_null = ~ Nutrient_no_level,
                                test = "Wald", boot = FALSE,
                                data = ps_corn,
                                fdr_cutoff = 0.05)
da_analysis$significant_taxa
#"ASV_1"  "ASV_3"  "ASV_4"  "ASV_67"
otu_to_taxonomy(OTU = da_analysis$significant_taxa, data = ps_corn)
summary(da_analysis)
corncob_stats <- summary(da_analysis)
plot <- plot(da_analysis, level = c("Genus"))
plot
ggsave(filename="~/Mote_nutrient_experiment/plots/da_nutrient.svg", plot=plot, device="svg",  width = 20, height = 6, dpi=500)

set.seed(1) #here controlling for the effect of No_level_weeks on dispersion
da_analysis <- differentialTest(formula = ~ No_level_weeks,
                                phi.formula = ~ No_level_weeks,
                                formula_null = ~ 1,
                                phi.formula_null = ~ No_level_weeks,
                                test = "Wald", boot = FALSE,
                                data = ps_corn,
                                fdr_cutoff = 0.05)
da_analysis$significant_taxa
#"ASV_1" "ASV_2" "ASV_4"
otu_to_taxonomy(OTU = da_analysis$significant_taxa, data = ps_corn)
plot(da_analysis, level = c("Genus"))
#not a great plot, too many panels

set.seed(1) #here controlling for the effect of exposure weeks on dispersion
da_analysis <- differentialTest(formula = ~ Exposure_weeks,
                                phi.formula = ~ Exposure_weeks,
                                formula_null = ~ 1,
                                phi.formula_null = ~ Exposure_weeks,
                                test = "Wald", boot = FALSE,
                                data = ps_corn,
                                fdr_cutoff = 0.05)
da_analysis$significant_taxa

#now let's look for taxa that are differentially variable (dispersion)
set.seed(1)
dv_analysis <- differentialTest(formula = ~ Nutrient_no_level,
                                phi.formula = ~ Nutrient_no_level,
                                formula_null = ~ Nutrient_no_level,
                                phi.formula_null = ~ 1,
                                test = "LRT", boot = FALSE,
                                data = ps_corn,
                                fdr_cutoff = 0.05)
dv_analysis$significant_taxa
#only ASV4 and ASV5

otu_to_taxonomy(OTU = da_analysis$significant_taxa, data = ps_corn)
da_analysis$p[1:5]
da_analysis$p_fdr[1:5]
plot(da_analysis)

which(is.na(da_analysis$p)) %>% names
#there are loads of taxa that could not be fit to a model

#subset phosphate treatment only
ps_corn <- ps %>%
  phyloseq::subset_samples(No_level_no_T0 %in% c("Phosphate")) %>%
  tax_glom("Genus")

set.seed(1)
da_analysis_phosphate <- differentialTest(formula = ~ Exposure_weeks,
                                phi.formula = ~ Exposure_weeks,
                                formula_null = ~ 1,
                                phi.formula_null = ~ Exposure_weeks,
                                test = "Wald", boot = FALSE,
                                data = ps_corn,
                                fdr_cutoff = 0.05)
da_analysis_phosphate$significant_taxa

plot(da_analysis_phosphate, level = c("Genus"))
ggsave(filename="~/Mote_nutrient_experiment/plots/da_phosphate.png", plot=plot, device="png",  width = 20, height = 3, dpi=500)


ps_corn_HL <- ps %>%
  phyloseq::subset_samples(Level %in% c("H", "L")) %>%
  tax_glom("Genus")
da_analysis <- differentialTest(formula = ~ Level,
                                phi.formula = ~ Level,
                                formula_null = ~ 1,
                                phi.formula_null = ~ Level,
                                test = "Wald", boot = FALSE,
                                data = ps_corn_HL,
                                fdr_cutoff = 0.05)
da_analysis$significant_taxa #no taxa significantly differentially abundant between treatment levels
