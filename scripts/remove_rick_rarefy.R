library("dada2")
library("seqinr")
library("biomformat")
library("data.table")
library("ggplot2")
library("ggthemes")
library("ampvis2")
library("cowplot")
library("phyloseq")

load(file = "ps_g50.RData")
ps = subset_samples(ps, SampleID != "Negative")

allTaxa = taxa_names(ps)
allTaxa <- allTaxa[!(allTaxa %in% "ASV_1")]
ps_no_rick = prune_taxa(allTaxa, ps)
# new phyloseq object with just the taxa you kept.

save(ps_no_rick, file = "ps_no_rick.RData")
rm(allTaxa)

#rarefy without rick...have to rarefy to 351 which sucks

# Need to convert from phyloseq to ampvis
av2_otutable <- data.frame(OTU = rownames(t(phyloseq::otu_table(ps)@.Data)),
                           t(phyloseq::otu_table(ps)@.Data),
                           phyloseq::tax_table(ps)@.Data,
                           check.names = F
)

#Extract metadata from the phyloseq object:
av2_metadata <- data.frame(phyloseq::sample_data(ps), 
                           check.names = F
)

av2_metadata <- cbind(rownames(av2_metadata), av2_metadata)

#Load the data with amp_load:
av2_obj <- amp_load(av2_otutable, av2_metadata)

rare_plot_amp <- amp_rarecurve(data = av2_obj, color_by = "Phenotype")
rare_curve_plot <- rare_plot_amp + ylab('Observed ASVs (count)') + 
  geom_vline(xintercept=min(sample_sums(ps)), linetype='dashed') +
  scale_colour_colorblind() +
  xlim(c(0, 35000))
plot(rare_plot_amp)
rare_curve_plot

min(sample_sums(ps_no_rick))
ntaxa(ps_no_rick) #1068

ps_rarefied <- rarefy_even_depth(ps_no_rick, sample.size = 351, rngseed = 999) #365 OTUs removed
ps_rarefied #703 taxa
sum(sample_sums(ps_rarefied))
sample_sums(ps_rarefied)

save(ps_rarefied, file = "ps_no_rick_rare.RData")

#CLR transform using microbiome package for beta div analyses
#CLR transform applies a pseudocount of min(relative abundance)/2 
#to exact zero relative abundance entries in OTU table before taking logs.
clr_no_rick <- microbiome::transform(ps_no_rick, 'clr')
save(clr_no_rick, file = "clr_no_rick.RData") #renamed, clr transformed data with no aquarickettsia