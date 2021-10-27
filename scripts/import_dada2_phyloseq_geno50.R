library("dada2")
library("seqinr")
library("data.table")
library("ggplot2")
library("ggthemes")
library("ampvis2")
library("phyloseq")
library("microbiome")

setwd("~/Mote_nutrient_experiment/data")

# functions
fast_melt = function(physeq){
  # supports "naked" otu_table as `physeq` input.
  otutab = as(otu_table(physeq), "matrix")
  if(!taxa_are_rows(physeq)){otutab <- t(otutab)}
  otudt = data.table(otutab, keep.rownames = TRUE)
  setnames(otudt, "rn", "taxaID")
  # Enforce character taxaID key
  otudt[, taxaIDchar := as.character(taxaID)]
  otudt[, taxaID := NULL]
  setnames(otudt, "taxaIDchar", "taxaID")
  # Melt count table
  mdt = melt.data.table(otudt, 
                        id.vars = "taxaID",
                        variable.name = "SampleID",
                        value.name = "count")
  # Remove zeroes, NAs
  mdt <- mdt[count > 0][!is.na(count)]
  # Calculate relative abundance
  mdt[, RelativeAbundance := count / sum(count), by = SampleID]
  if(!is.null(tax_table(physeq, errorIfNULL = FALSE))){
    # If there is a tax_table, join with it. Otherwise, skip this join.
    taxdt = data.table(as(tax_table(physeq, errorIfNULL = TRUE), "matrix"), keep.rownames = TRUE)
    setnames(taxdt, "rn", "taxaID")
    # Enforce character taxaID key
    taxdt[, taxaIDchar := as.character(taxaID)]
    taxdt[, taxaID := NULL]
    setnames(taxdt, "taxaIDchar", "taxaID")
    # Join with tax table
    setkey(taxdt, "taxaID")
    setkey(mdt, "taxaID")
    mdt <- taxdt[mdt]
  }
  return(mdt)
}

# phyloseq object output from Dada2
ps <- readRDS("ps_object_50.rds")
seqtab_50 <- readRDS("seqtab_50.rds")
# sequences object made from Dada2

# exporting sequences to a fasta file for import into qiime
uniqueSeqs <- as.list(colnames(seqtab_50))
write.fasta(uniqueSeqs, uniqueSeqs, "uniqueSeqs_50.fasta")

# import metadata and merge into phyloseq object
mapfile = "Mapping-file-full-renamed.txt"
map = import_qiime_sample_data(mapfile)
sample_data(ps) <- map

# export taxonomy to import into qiime2
tax<-as(tax_table(ps),"matrix")
tax_cols <- c("Kingdom", "Phylum", "Class", "Order","Family","Genus", "Species")
tax<-as.data.frame(tax)
tax$taxonomy<-do.call(paste, c(tax[tax_cols], sep=";"))
for(co in tax_cols) tax[co]<-NULL
write.table(tax, "taxonomy_50.txt", quote=FALSE, col.names=FALSE, sep="\t")

# summary of data
ps
ps@sam_data

ntaxa(ps)
nsamples(ps)
rank_names(ps)
sample_names(ps)[1:5]
sample_variables(ps)

# remove mitochondria and chloroplasts, is.na important becuase if not included
# this command will also remove all Family = NA or Order = NA
ps_with_mito = subset_taxa(ps, (Order!="Chloroplast") | is.na(Order)) #lost 74 taxa
ps_no_mito = subset_taxa(ps_with_mito, (Family!="Mitochondria") | is.na(Family)) #lost 14 taxa
ps_no_Eukaryota = subset_taxa(ps_no_mito, (Kingdom!="Eukaryota") | is.na(Kingdom)) #lost 11 taxa

sum(sample_sums(ps)) - sum(sample_sums(ps_no_Eukaryota))
ps = ps_no_Eukaryota #lost 99 taxa


summary(sample_sums(ps))
#mean was 21,901, min was 5098

summary(taxa_sums(ps)) # 5 was first quantile, mean was 2122
# Filter on prevalence or total counts
pst = fast_melt(ps)
prevdt = pst[, list(Prevalence = sum(count > 0), 
                    TotalCounts = sum(count)),
             by = taxaID]
uniques = prevdt[(Prevalence <2 & TotalCounts >5), taxaID] #580 taxa present in only one sample
keepTaxa = prevdt[(Prevalence >=0 & TotalCounts >5), taxaID]
ps_pruned = prune_taxa(keepTaxa,ps)
summary(taxa_sums(ps_pruned)) #new min is 6, new mean is 3072
ntaxa(ps_pruned) #1069

ps_pruned
sample_sums(ps_pruned)
summary(sample_sums(ps_pruned)) #mean was 21890, min stayed the same
sum(sample_sums(ps_pruned)) #3283526
sum(sample_sums(ps)) - sum(sample_sums(ps_pruned)) #only 1611 reads lost by pruning

save(ps, file = "ps_unfiltered_geno50.RData")

ps <- ps_pruned
save(ps, file = "ps_pruned_geno50.RData")

summary_tab <- data.frame(init=(as.matrix(sample_sums(ps_50)))[,1],
                           chloros_removed=(as.matrix(sample_sums(ps_with_mito)))[,1],
                           mitos_removed=(as.matrix(sample_sums(ps_no_mito)))[,1],
                           euks_removed=(as.matrix(sample_sums(ps_no_Eukaryota)))[,1],
                           pruned=(as.matrix(sample_sums(ps_pruned)))[,1])
write.table(summary_tab, "reads_lost_phyloseq_geno50.txt", quote=FALSE, col.names=TRUE, sep="\t")
rm(ps_with_mito, ps_no_Eukaryota, ps_no_mito, ps_50)

#plot otus by sequencing depth
observed <- estimate_richness(ps, measures = c('Observed'))
explore.df <- cbind(observed, sample_sums(ps), sample_data(ps)$Nutrient_no_level)
colnames(explore.df) <- c('Observed', 'Sample_Sums', 'Nutrient_no_level')
observed_mean <- mean(explore.df$Observed)
sample_sum_mean <- mean(explore.df$Sample_Sums)
ggplot(data = explore.df, aes(x = Sample_Sums, y = Observed, color = Nutrient_no_level)) + 
  geom_point() +
  geom_smooth(method="auto", se=TRUE, fullrange=FALSE, level=0.95, 
              inherit.aes = F, mapping = aes(Sample_Sums, Observed),
              data = explore.df) +
  ylab("Observed OTUs") +
  scale_colour_colorblind()

# Let's use ampvis2 again so we can easily make a rarefaction curve

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

# RARE CURVE
rare_plot_amp <- amp_rarecurve(data = av2_obj, color_by = "Nutrient_no_level")
rare_curve_plot <- rare_plot_amp + ylab('Observed ASVs (count)') + 
  geom_vline(xintercept=min(sample_sums(ps)), linetype='dashed') +
  scale_colour_colorblind() +
  xlim(c(0, 45000))
plot(rare_plot_amp)
rare_curve_plot

#If we rarefy to 5098
rare_plot_amp <- amp_rarecurve(data = av2_obj, color_by = "Nutrient_no_level")
rare_curve_plot <- rare_plot_amp + ylab('Observed ASVs (count)') + 
  geom_vline(xintercept=c(5098), linetype='dashed') +
  scale_colour_colorblind() +
  xlim(c(0, 45000))
plot(rare_plot_amp)
rare_curve_plot
ggsave(filename="rare_curve_plot.png", plot=rare_curve_plot, device="png", dpi=500)


summary(explore.df$Sample_Sums)
length(explore.df$Sample_Sums)

#Actual rarefying done after renaming, see taxa_rename script
