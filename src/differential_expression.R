library(edgeR)

## SETUP DATA
# - read in raw data
df <- read.csv("../outputs/raw_counts.csv", row.names=1)
colnames(df) <- substr(colnames(df), 2, 6)
# - transform into dgelist class
dge <- DGEList(counts=df, genes=row.names(df))

## FILTER DATA
# - transform into CPM (sequencing depth normalization)
values.cpm <- cpm(dge, log=FALSE)  # count / library size * 1e6
# - filter for genes in all samples with cpm >= 1
valid.genes <- which(rowSums(values.cpm >= 1) >= dim(dge$samples)[1])
dge <- dge[valid.genes,]  # actually filter
values.cpm <- cpm(dge, log=FALSE) 
# transform into TMM (RNA composition)
dge <- calcNormFactors(dge, method = "TMM")
values.tmm <- t(t(as.matrix(dge)) / (dge$samples$lib.size * dge$samples$norm.factors) * 1e6)

## PLOT PCA
library(ggplot2)
library(ggrepel)
# function
analyze_pca <- function(data, label) {
  pcs <- prcomp(t(data), center=TRUE, rank=50)
  pcs.df <- as.data.frame(pcs$x)  # copy over
  pcs.df$n_counts <- colSums(df)  # add on an n-counts column
  pcs.df$sample <- row.names(pcs.df)  # add on sample name column
  var.perc <- data.frame(sapply(pcs$sdev, function(stddev) {  # get perc var
    return(stddev**2/sum((pcs$sdev)**2)*100)
  }), row.names = colnames(pcs$x))
  names(var.perc) <- c("var.explain")
  var.perc$pc <- row.names(var.perc)
  var.perc$pc.num <- as.numeric(substr(var.perc$pc, 3, 4))
  p1 <-
    ggplot(data=var.perc, mapping=aes(x=pc.num, y=var.explain)) +
    geom_col() + geom_point() + geom_line() + 
    scale_x_discrete(name ="Principal Components", limits=row.names(var.perc)) +
    ylab("Percent Variance Explained") +
    ggtitle(paste0(label,"-Based PCA Variance Explained Analysis")) +
    theme(panel.background=element_rect(fill="white", colour="white"),
          axis.line=element_line(colour="black", size=0.5))
  print(p1)
  p2 <- 
    ggplot(data=pcs.df, mapping=aes(x=PC1, y=PC2, color=n_counts,
                                    label=sample)) +
    geom_point() + geom_text_repel() + 
    xlab(paste0("PC1 (", round(var.perc["PC1","var.explain"], 2), "%)")) +
    ylab(paste0("PC2 (", round(var.perc["PC2","var.explain"], 2), "%)")) +
    scale_color_gradient2(low="forestgreen", mid="darkgray", high="red",
                          midpoint=median(pcs.df$n_counts)*0.8) +
    labs(color="N-Counts") + ggtitle(paste0(label,"-Based PCA")) +
    theme(panel.background=element_rect(fill="white", colour="white"),
          axis.line=element_line(colour="black", size=0.5))
  print(p2)
}
# plots
analyze_pca(values.cpm, "CPM")
analyze_pca(values.tmm, "TMM")

## CALL DIFFERENTIAL
# set up samples
valid.samples <- c("4306L", "4301L", "4202L", "4397L")
dge <- dge[,valid.samples]  # subset data
treatments <- as.factor(c(rep("cancer:DKO+KI",2), rep("cancer:DKO+RV",2)))
dge$samples$group <- treatments
designMat <- model.matrix(~treatments)
# calculate dispersions
dge <- estimateDisp(dge, design=designMat)
dge <- estimateCommonDisp(dge, design=designMat)
dge <- estimateTagwiseDisp(dge, design=designMat)
plotBCV(dge)
# perform test
comparison <- c("cancer:DKO+RV","cancer:DKO+KI")
test <- exactTest(dge, pair=comparison)
test.table <- topTags(test, n=nrow(test$table))$table  # add on fdr
# create plotting dataframe
# - get data
test.plot <- test.table  # copy over
test.plot$nlog10pval <- -log10(test.plot$PValue)
test.plot$nlog10fdr <- -log10(test.plot$FDR)
# - find significant
test.plot.de <- test.plot[test.plot$FDR < 0.001,]
valid.genes <- c("Tbx21","Ctla4","Ifng","Gzma","Cxcr3")
test.plot.de <- test.plot.de[valid.genes,]
test.plot$names <- ""
test.plot[valid.genes,"names"] <- valid.genes
p1 <-
  ggplot() + ggtitle(paste0("Volcano Plot of ",comparison[2]," - ",comparison[1])) +
  geom_point(data=test.plot, mapping=aes(x=logFC, y=nlog10pval, color=nlog10fdr)) +
  geom_text_repel(data=test.plot, mapping=aes(x=logFC, y=nlog10pval, label=names),
                  max.overlaps=Inf, box.padding=0.5, force_pull=3, force=5) +
  geom_point(data=test.plot.de, mapping=aes(x=logFC, y=nlog10pval), color="red") +
  ylab("-log10(P-Value") + xlab(paste0("logFC(",comparison[2],"/",comparison[1],")")) +
  theme(panel.background=element_rect(fill="white", colour="white"),
        axis.line=element_line(colour="black", size=0.5)) + labs(color="-log10(FDR)")
p1
# save results
write.csv(test.table,"../outputs/diffHL.CDkoKi.CDkoRv.csv")






