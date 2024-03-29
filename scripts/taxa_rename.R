library("phyloseq")

setwd("~/Mote_nutrient_experiment/data")

load(file = "ps_unfiltered_geno50.RData") #unpruned

mapfile = "Mapping-file-full-renamed.txt"
map = import_qiime_sample_data(mapfile)
sample_data(ps) <- map
ps
#remove negatives
ps = subset_samples(ps, SampleID != "Negative")
tax.clean <- data.frame(tax_table(ps))
write.table(tax.clean, "taxonomy_unrare_fullseq.txt", sep = "\t")

taxa_names(ps) <- paste0("ASV_", seq(ntaxa(ps)))
tax.clean <- data.frame(tax_table(ps))
write.table(tax.clean, "taxonomy_unrare.txt", sep = "\t")

new_names <- paste0("ASV_", seq(ntaxa(ps)))

for (i in 1:7){ tax.clean[,i] <- as.character(tax.clean[,i])}
tax.clean[is.na(tax.clean)] <- ""

for (i in 1:nrow(tax.clean)){

   if (tax.clean[i,6] == ""){
    family <- paste("Family_", tax.clean[i,5], sep = "")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == ""){
    tax.clean$Species[i] <- paste("Genus",tax.clean$Genus[i], sep = "_")
  }
}

for (i in 1:nrow(tax.clean)){
  
  if (tax.clean[i,6] == "Family_"){
    family <- paste("", tax.clean[i,5], sep = "")
    tax.clean[i, 6:7] <- family
  }
  if (tax.clean[i,6] == "MD3-55"){
    tax.clean[i, 6:7] <- "Aquarickettsia"
  }
}

for (i in 1:nrow(tax.clean)){
  
  if (tax.clean[i,6] == ""){
    family <- paste0("Unclassified_",new_names[i], sep = "") 
    tax.clean[i, 5:7] <- family
  }
}


tax_table(ps) <- as.matrix(tax.clean)
summary(sample_sums(ps)) #min is 5098
ps_rarefied <- rarefy_even_depth(ps, sample.size = 5098, rngseed = 999)
ps_rarefied #1259 taxa
sum(sample_sums(ps_rarefied)) 
sample_sums(ps_rarefied)
save(ps_rarefied, file = "ps_rare_g50.RData") #rarefied, unpruned

load(file = "ps_pruned_geno50.RData") #pruned

mapfile = "Mapping-file-full-renamed.txt"
map = import_qiime_sample_data(mapfile)
sample_data(ps) <- map
ps
#remove negatives
ps = subset_samples(ps, SampleID != "Negative")
tax.clean <- data.frame(tax_table(ps))
#write.table(tax.clean, "taxonomy_unrare_fullseq.txt", sep = "\t")

taxa_names(ps) <- paste0("ASV_", seq(ntaxa(ps)))
tax.clean <- data.frame(tax_table(ps))
#write.table(tax.clean, "taxonomy_unrare.txt", sep = "\t")

new_names <- paste0("ASV_", seq(ntaxa(ps)))

for (i in 1:7){ tax.clean[,i] <- as.character(tax.clean[,i])}
tax.clean[is.na(tax.clean)] <- ""

for (i in 1:nrow(tax.clean)){
  
  if (tax.clean[i,6] == ""){
    family <- paste("Family_", tax.clean[i,5], sep = "")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == ""){
    tax.clean$Species[i] <- paste("Genus",tax.clean$Genus[i], sep = "_")
  }
}

for (i in 1:nrow(tax.clean)){
  
  if (tax.clean[i,6] == "Family_"){
    family <- paste("", tax.clean[i,5], sep = "")
    tax.clean[i, 6:7] <- family
  }
  if (tax.clean[i,6] == "MD3-55"){
    tax.clean[i, 6:7] <- "Aquarickettsia"
  }
}

for (i in 1:nrow(tax.clean)){
  
  if (tax.clean[i,6] == ""){
    family <- paste0("Unclassified_",new_names[i], sep = "") 
    tax.clean[i, 5:7] <- family
  }
}


tax_table(ps) <- as.matrix(tax.clean)

#CLR transform using microbiome package for beta div analyses
#CLR transform applies a pseudocount of min(relative abundance)/2 
#to exact zero relative abundance entries in OTU table before taking logs.
clr <- microbiome::transform(ps, 'clr') #pruned, unrarefied, clr transformed

ps_rarefied <- rarefy_even_depth(ps, sample.size = 5098, rngseed = 999) 
ps_rarefied 
sum(sample_sums(ps_rarefied)) 
sample_sums(ps_rarefied)

#relative abundance transform
ps_rel <- microbiome::transform(ps_rarefied, "compositional")


# save ps objects
save(ps, file = "ps_rename_g50.RData") #unrarefied, pruned
save(ps_rel, file = "ps_rel_g50.RData") #rel abundance, rarefied, pruned
save(clr, file = "ps_clr_g50.RData") #clr transformed, pruned
