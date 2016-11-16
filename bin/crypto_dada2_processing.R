library(dada2)
library(ShortRead)
library(ggplot2)
library(phyloseq)
library(PMCMR)

path <- "~/Documents/Projects/crypto/data/reads"
fns <- list.files(path)
fns

### Load forward and reverse reads
fastqs <- fns[grepl(".fastq$", fns)]
fastqs <- sort(fastqs) # Sort ensures forward/reverse reads are in same order
fnFs <- fastqs[grepl("_R1", fastqs)] # Just the forward read files
fnRs <- fastqs[grepl("_R2", fastqs)] # Just the reverse read files
# Get sample names from the first part of the forward read filenames
sample.names <- sapply(strsplit(fnFs, "_"), `[`, 1)
# Fully specify the path for the fnFs and fnRs
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)

# visualize quality profile. can also be done w/ reverse reads
plotQualityProfile(fnFs[[15]])

### Filtering and trimming
# Make directory and filenames for the filtered fastqs
filt_path <- file.path(path, "filtered")
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))
# Filter
for(i in seq_along(fnFs)) {
  fastqFilter(fnFs[i], filtFs[i],
                    trimLeft=20, truncLen=200, 
                    maxN=0, maxEE=2, truncQ=2, 
                    compress=TRUE, verbose=TRUE)
}

# In the above function, trimleft describes the number of bases to skip at the beagining of each read as (forward, reverse).
# truncLen is the same, except describes the base number to truncate at. Choice for truncLen is justified by the point at which
# error rates tend to start increasing.

### Dereplication
derepFs <- derepFastq(filtFs, verbose=TRUE)
#derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
#names(derepRs) <- sample.names

derepFs[[1]]

### Sample inference

dadaFs <- dada(derepFs, err=NULL, selfConsist = TRUE,multithread = TRUE)
plotErrors(dadaFs[[2]], nominalQ=TRUE)

### construct sequence table from forward reads
seqtab <- makeSequenceTable(dadaFs)
dim(seqtab)
table(nchar(colnames(seqtab)))

### Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)


### Assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, paste(path,"/rdp_train_set_14.fa.gz",sep = ""))
taxa.plus <- addSpecies(taxa, paste(path,"/rdp_species_assignment_14.fa.gz",sep = ""))
colnames(taxa.plus) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus","Species")
unname(taxa.plus)


### analyze w/ phyloseq
# load sample metadata
sample_data <- read.table('~/Documents/Projects/crypto/data/crypto_sample_metadata.csv',sep="\t",header=TRUE)
# Pass data from DADA2 objects into phyloseq
rownames(sample_data) <- sample_data$SAMPLE_ID
samples.out <- rownames(seqtab.nochim)
sample_data <- sample_data[samples.out,] #re-order dataframe to match seqtab.nochim sample order
diet_description = sample_data$Diet_and_des
days_on_diet = sample_data$Days.on.diet

# save seqtab.nochim and taxa.plus for analysis w/ metagenomeSeq
write.table(t(seqtab.nochim), file='~/Documents/Projects/crypto/data/seqtab.nochim.tsv', quote=FALSE, sep='\t', col.names = NA)
write.table(taxa.plus, file='~/Documents/Projects/crypto/data/taxa.plus.tsv', quote=FALSE, sep='\t', col.names = NA)





samdf <- data.frame(sample=samples.out, description=diet_description, day=days_on_diet)
rownames(samdf) <- samples.out

seqtab.nochim.taxa <- seqtab.nochim # copy to create a new object that has taxa information
colnames(seqtab.nochim.taxa) <- 1:nrow(taxa.plus) # add the taxa information
taxa.for.table <- taxa.plus
rownames(taxa.for.table) <- 1:nrow(taxa.plus)
ps <- phyloseq(otu_table(seqtab.nochim.taxa, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa.for.table))
ps

# Filter sequences that don't have more than 10 counts in any sample
ps.count_thresh <- filter_taxa(ps, function(x) max(x) > 10, TRUE)
ps.count_thresh

library(DESeq2)

# convert day values to factor for DESeq
sample_data(ps.count_thresh)$day <- factor(sample_data(ps.count_thresh)$day)
# Create separate object for comparisons at each time point
d13_subset <- prune_samples(sample_data(ps.count_thresh)$day == 13, ps.count_thresh)
d14_subset <- prune_samples(sample_data(ps.count_thresh)$day == 14, ps.count_thresh)
d20_subset <- prune_samples(sample_data(ps.count_thresh)$day == 20, ps.count_thresh)

# D13 comparison
deseq2_input_d13 <- phyloseq_to_deseq2(d13_subset,~description)
deseq2_output_d13 <- DESeq(deseq2_input_d13, test="Wald", fitType="parametric")
deseq2_results_d13 <- results(deseq2_output_d13, cooksCutoff = FALSE)

alpha = 0.05
sigtab_d13 = deseq2_results[which(deseq2_results_d13$padj < alpha), ]
sigtab_d13 = cbind(as(sigtab_d13, "data.frame"), as(tax_table(ps.count_thresh)[rownames(sigtab_d13), ], "matrix"))
head(sigtab_d13)

# D14 comparison
deseq2_input_d14 <- phyloseq_to_deseq2(d14_subset,~description)
deseq2_output_d14 <- DESeq(deseq2_input_d14, test="Wald", fitType="parametric")
deseq2_results_d14 <- results(deseq2_output_d14, cooksCutoff = FALSE)

alpha = 0.05
sigtab_d14 = deseq2_results_d14[which(deseq2_results_d14$padj < alpha), ]
sigtab_d14 = cbind(as(sigtab_d14, "data.frame"), as(tax_table(ps.count_thresh)[rownames(sigtab_d14), ], "matrix"))
head(sigtab_d14)

# D20 comparison
deseq2_input_d20 <- phyloseq_to_deseq2(d20_subset,~description)
deseq2_output_d20 <- DESeq(deseq2_input_d20, test="Wald", fitType="parametric")
deseq2_results_d20 <- results(deseq2_output_d20, cooksCutoff = FALSE)

alpha = 0.05
sigtab_d20 = deseq2_results_d20[which(deseq2_results_d20$padj < alpha), ]
sigtab_d20 = cbind(as(sigtab_d20, "data.frame"), as(tax_table(ps.count_thresh)[rownames(sigtab_d20), ], "matrix"))
head(sigtab_d20)


### Get sequences that were significantly different in each comparison to generate ecological heatmaps
sig_seqs <- union(rownames(sigtab_d13),rownames(sigtab_d14))
sig_seqs <- union(sig_seqs,rownames(sigtab_d20))
# Remove weaned and 6d on-diet mice
ps.remove_weaned <- prune_samples(sample_data(ps.count_thresh)$description != "dPD.weaned",ps.count_thresh)
ps.remove_pre_inf <- prune_samples(sample_data(ps.remove_weaned)$description != "dPD.6d.post.weaning",ps.remove_weaned)

# Transform sequence table from counts to relative abundance
ps.remove_pre_inf.rel <- transform_sample_counts(ps.remove_pre_inf, function(x) x/sum(x))

# Generate OTU tables for each comparison, jsut in case. Here, we don't use these.
d13_OTU_table <- otu_table(transform_sample_counts(d13_subset, function(x) x/sum(x)))[,rownames(sigtab_d13)]
colnames(d13_OTU_table) <- tax_table(d13_subset)[,6][rownames(sigtab_d13),]
d14_OTU_table <- otu_table(transform_sample_counts(d14_subset, function(x) x/sum(x)))[,rownames(sigtab_d14)]
colnames(d14_OTU_table) <- tax_table(d14_subset)[,6][rownames(sigtab_d14),]
d20_OTU_table <- otu_table(transform_sample_counts(d20_subset, function(x) x/sum(x)))[,rownames(sigtab_d20)]
colnames(d20_OTU_table) <- tax_table(d20_subset)[,6][rownames(sigtab_d20),]

### Generate OTU tables for all samples, keeping any OTU that was significantly differentially abundance in any of the 3 comparisons.
sig.thresh_OTU_table <- otu_table(transform_sample_counts(ps.count_thresh, function(x) x/sum(x)))[,sig_seqs]

### Generate heatmaps for paper
# Heatmap for reference for Jordi during cross-correlation analysis
p <- plot_heatmap(prune_taxa(sig_seqs,ps.remove_pre_inf.rel),"NMDS","bray",sample.label = "description","Genus", first.sample = "Plate2-B4")
ggsave(file = "~/Documents/Projects/crypto/results/seq_heatmap.svg", plot=p, width = 10, height = 10)

# Extract sample and taxa ordering as determined by plot_heatmap()
taxa_order <- names(p$plot_env$labvec)
sample_order <- p$plot_env$sample.order
# Save heatmap data for Jordi's correlation analysis
sig_seq_data <- prune_taxa(sig_seqs,ps.remove_pre_inf.rel)
final_seqs <- otu_table(sig_seq_data)
final_taxa <- tax_table(sig_seq_data)
final_metadata <- sample_data(sig_seq_data)

final_seqs <- final_seqs[sample_order,taxa_order]
final_taxa <- final_taxa[taxa_order,]
final_metadata <- final_metadata[sample_order,]

write.table(final_seqs, file = "~/Documents/Projects/crypto/results/crypto_final_seq_relative_abundance.tsv", sep = "\t", quote = FALSE)
write.table(final_taxa, file = "~/Documents/Projects/crypto/results/crypto_final_taxa.tsv", sep = "\t", quote = FALSE)
write.table(final_metadata, file = "~/Documents/Projects/crypto/results/crypto_final_metadata.tsv", sep = "\t", quote = FALSE)

### plotting summary of logfold changes for one of the comparisons: example w/ d20 comparison:

# 
# library("ggplot2")
# theme_set(theme_bw())
# scale_fill_discrete <- function(palname = "Set1", ...) {
#   scale_fill_brewer(palette = palname, ...)
# }
# # Phylum order
# x = tapply(sigtab_d20$log2FoldChange, sigtab_d20$Phylum, function(x) max(x))
# x = sort(x, TRUE)
# sigtab_d20$Phylum = factor(as.character(sigtab_d20$Phylum), levels=names(x))
# # Genus order
# x = tapply(sigtab_d20$log2FoldChange, sigtab_d20$Genus, function(x) max(x))
# x = sort(x, TRUE)
# sigtab_d20$Genus = factor(as.character(sigtab_d20$Genus), levels=names(x))
# ggplot(sigtab_d20, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
#   theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

