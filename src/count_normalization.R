library(edgeR)
library(GenomicFeatures)

# read in counts
df <- read.csv("../outputs/raw_counts.csv", row.names=1)
# transform into TPM
# - load in GTF file
txdb <- makeTxDbFromGFF("../data/reference/mm10.refGene.gtf", format="gtf")
# - collect the exons per gene id
exons.list.per.gene <- exonsBy(txdb, by="gene")
# - for each gene, reduce all the exons to a set of non overlapping exons,
#   calculate their lengths (widths) and sum then
exonic.gene.sizes <- sum(width(reduce(exons.list.per.gene)))
exonic.gene.sizes <- exonic.gene.sizes[row.names(df)]  # convert to same order
# transform
dge <- DGEList(counts=df, genes=row.names(df))
dge$genes$Length <- exonic.gene.sizes[dge$genes$genes]
# transform into CPM - seq.depth
values.cpm <- cpm(dge, log=FALSE)  # count / library size
# transform into RPKM - seq.depth, gene.length
values.rpkm <- rpkm(dge, log=FALSE)  # count / library size / gene length
# transform into TPM - seq.depth, gene.length (better than RPKM)
values.tpm <- (as.matrix(dge) / (exonic.gene.sizes / 1e3))  # get counts / kilobase = transcripts
values.tpm <- t(t(values.tpm) / colSums(values.tpm)) * 1e6  # get RPK / million transcripts
# transform into TMM - RNA.composition
dge <- calcNormFactors(dge, method = "TMM")
values.tmm <- t(t(as.matrix(dge)) / (dge$samples$lib.size * dge$samples$norm.factors) * 1e6)
# transform into TPM TMM - seq.depth, gene.length, RNA.composition
values.tpm_tmm <- t(t(values.tpm) / dge$samples$norm.factors)
# write data
write.csv(values.cpm, file="../outputs/normalized_counts_CPM.csv")
write.csv(values.rpkm, file="../outputs/normalized_counts_RPKM.csv")
write.csv(values.tpm, file="../outputs/normalized_counts_TPM.csv")
write.csv(values.tmm, file="../outputs/normalized_counts_TMM.csv")
write.csv(values.tpm_tmm, file="../outputs/normalized_counts_TPM_TMM.csv")


