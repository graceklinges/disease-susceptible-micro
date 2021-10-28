library(dada2); packageVersion("dada2")
library(phyloseq)

setwd("~/Mote_nutrient_experiment/raw_read_old_data")

# Location of raw reads
pathF <- "F_trimmed"
pathR <- "R_trimmed"

# Double check you've set the right paths - should list your read 1 files
list.files(pathF)

# Assumes forward and reverse fastq filenames have format: SAMPLENAME_R1.fq and SAMPLENAME_R2.fq
fnFs = sort(intersect(list.files(pathF, pattern = "-50-"), list.files(pathF, pattern = "_R1_001.fastq.gz")))
fnRs = sort(intersect(list.files(pathR, pattern = "-50-"), list.files(pathR, pattern = "_R2_001.fastq.gz")))

length(fnFs) #153
length(fnRs) #153

# Plot read quality- you expect read quality to drop off towards the end
# plotQualityProfile(fnFs[1:4])
# plotQualityProfile(fnRs[1:4])

sample.names <- sapply(strsplit(basename(fnFs), "trimmed_lane1-s[0-9]{3}-index-[A-Z]{8}-[A-Z]{8}-"), `[`, 2)
sample.names <- sapply(strsplit(sample.names, "_"), `[`, 1)

sample.namesR <- sapply(strsplit(basename(fnRs), "trimmed_lane1-s[0-9]{3}-index-[A-Z]{8}-[A-Z]{8}-"), `[`, 2)
sample.namesR <- sapply(strsplit(sample.names, "_"), `[`, 1)
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")

write.table(sample.names, file = "names.txt", sep="\t")

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(pathF, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(pathR, "filtered", paste0(sample.namesR, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names


# Filter and trim
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(250,200),
                     maxN=0, maxEE=c(2,2), truncQ=10, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE, matchIDs = TRUE)
write.table(out, file = "filter_out.txt", sep="\t")

keep <- out[,"reads.out"] > 30 # Or other cutoff. This removed 3 samples
filtFs <- filtFs[keep]
filtRs <- filtRs[keep]
filtFs <- filtFs[!is.na(filtFs)]
filtRs <- filtRs[!is.na(filtRs)]

sample.names <- names(filtFs)

set.seed(100)
# Learn forward error rates
errF <- learnErrors(filtFs, nbases=1e8, multithread=TRUE)
# Learn reverse error rates
errR <- learnErrors(filtRs, nbases=1e8, multithread=TRUE)

# If you like, visualize these error rates. You want the estimated error rates (black lines) to be a good 
# fit to the observed error rates (points) and for the error rates drop with increased quality
# plotErrors(errF, nominalQ=TRUE)
# plotErrors(errR, nominalQ=TRUE)

# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names

# Note: if this loop doesn't work, you've done something wrong. You've probably made an error earlier in the code, 
# probably with the location of your forward and reverse reads and it's generated duplicate sample names
# Normal to get warning message that there are duplicate sequences
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE) # core sample interference algorithm
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE) # core sample interference algorithm
  merger <- mergePairs(ddF, derepF, ddR, derepR) # Merge paired end reads
  mergers[[sam]] <- merger
}

# Make sequence table from merged reads
st.all <- makeSequenceTable(mergers) # Normal to get warning message saying the sequences being tabled vary in length

# Remove any ASVs that are considerably off target length
seqtab_trimmed <- st.all[,nchar(colnames(st.all)) %in% seq(250,254)] #removed 73 ASVs

# Inspect distribution of read lengths after removal of off-target reads
table(nchar(getSequences(seqtab_trimmed)))
#250  252  253  254 
#61   65 1676   73 

# Remove chimeric sequences
seqtab_50 <- removeBimeraDenovo(seqtab_trimmed, method="consensus", multithread=TRUE, verbose = T) #228 bimeras
sum(seqtab_trimmed)-sum(seqtab_50) # How many chimeric reads were removed? 68,101
sum(seqtab_50)/sum(seqtab_trimmed) #0.9803926

# Save chimera-free ASV table as downstream tasks may cause R to crash
saveRDS(seqtab_50, "seqtab_50.rds")

# Assign taxonomy based on silva reference database at genus level, you must have the appropriate Silva database downloaded
tax_silva_50 <- assignTaxonomy(seqtab_50, "/Volumes/Grace\ External\ 2/Mote_Nutrient_Exp/silva_nr_v132_train_set.fa.gz", multithread=TRUE)

# Assign taxonomy based on silva reference database at species (100%) level
silva_sp_50 <- addSpecies(tax_silva_50, "/Volumes/Grace\ External\ 2/Mote_Nutrient_Exp/silva_species_assignment_v132.fa.gz")

# Export sequence table with genus and species assignments as phyloseq objects
ps_object_50 <- phyloseq(otu_table(seqtab_50, taxa_are_rows=FALSE), tax_table(silva_sp_50))

# Save as RDS object
saveRDS(ps_object_50, file = "ps_object_50.rds")
