library(metagenomeSeq)

### An attempt at using metagenomeSeq. Abandoned in favor of DESeq2.

dataDirectory <- "~/Documents/Projects/crypto/data"
seqs <- loadMeta(file.path(dataDirectory, "seqtab.nochim.tsv"),sep="\t")
taxa <- read.table(file.path(dataDirectory, "taxa.plus.tsv"),
                  stringsAsFactors = FALSE, header=TRUE) #Important: read.delim, as in tutorial, will get rid of row names.

# load metadata
metadata = loadPhenoData(file.path(dataDirectory, "crypto_sample_metadata.csv"),
                     tran = TRUE, sep = "\t")
ord = match(colnames(seqs$counts), rownames(metadata))
metadata = metadata[ord, ]

# Create MRexperiment object
phenotypeData = AnnotatedDataFrame(metadata)
OTUdata = AnnotatedDataFrame(taxa)
obj = newMRexperiment(seqs$counts,phenoData=phenotypeData,featureData=OTUdata)

# Filter based on depth and presence
obj.filt <- filterData(obj, present = 4, depth = 1000) # present in at least 4 samples (# of samples in the smallest group),
                                                       # at least 1000 seqs per sample

# Normalization
p <- cumNormStatFast(obj.filt) # calculates percentile to be used for normalization
obj.filt <- cumNorm(obj.filt, p = p)

# export normalized counts and sample statistics:
mat = MRcounts(obj.filt, norm = TRUE, log = TRUE)[1:5, 1:5]
exportMat(mat, file = file.path(dataDirectory, "normalized_seq_counts.tsv"))
exportStats(obj.filt[, 1:5], file = file.path(dataDirectory,
                                              "sample_stats.tsv"))

pd <- pData(obj.filt)
settings = zigControl(maxit = 100, verbose = FALSE)
diet.des <- pd$Diet_and_des
mod <- model.matrix(~diet.des)
colnames(mod) <- levels(diet.des)

# fit ZIG model
res <- fitZig(obj = obj.filt, mod = mod, control = settings)
eff.samples <- calculateEffectiveSamples(res)
obj.filt.effective <- obj.filt[which(eff.samples>mean(eff.samples)),]
effect.p <- cumNormStat(obj.filt.effective, pFlag = TRUE, main = "Trimmed crypto data")
obj.filt.effective <- cumNorm(obj.filt.effective, p = effect.p)

diet.des <- pData(obj.filt.effective)
diet.des <- diet.des$Diet_and_des
mod <- model.matrix(~diet.des)
colnames(mod) <- levels(diet.des)

res <- fitZig(obj = obj.filt.effective, mod = mod, control = settings)
zigFit <- res$fit
finalMod <- res$fit$design

contrast.matrix = makeContrasts(d20 = dPD.20d.post.weaning.13d.post.crypto - dPD.20d.post.weaning,
                                d14 = dPD.14d.post.weaning.7d.post.crypto - dPD.14d.post.weaning,
                                d13 = dPD.13d.post.weaning.6d.post.crypto - dPD.13d.post.weaning,
                                levels = finalMod)
fit2 = contrasts.fit(zigFit, contrast.matrix)
fit2 = eBayes(fit2)

# get info from fit model
tbl.d13 <- limma::topTable(fit2, ceof=c(1,2),number=dim(fit2)[1], adjust.method = "fdr")#, resort.by = "p", p.value = 0.05)


### model export

taxa = sapply(fData(obj.filt.effective)[,4],
              function(i) {
                i[length(i)]
              })
head(MRcoefs(res, taxa = taxa, coef = 5))
# Make MRtable!!!

final.table <- MRtable(fit2,by=d20,number=72)
final.table
