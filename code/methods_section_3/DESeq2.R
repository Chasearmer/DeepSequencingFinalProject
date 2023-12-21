# initialize
rm(list=ls()) ; gc();
set.seed(1337)

# load library
library(DESeq2)
library(tibble)
library(readr)
library(dplyr)

# load the count data
countdata <- read.csv(
  "/path/to/gene_count_matrix_2.csv",
  header = TRUE,
  sep = ",",
  row.names = "gene_id"
)
# load the metadata
coldata <-  read.csv(
  "/path/to/metadata_GSE147849.csv",
  header = TRUE,
  sep = ",",
)

coldata <- coldata[, colnames(coldata) %in% c("Run","treatment")]
coldata <- subset(coldata, Run %in% colnames(countdata))

# create DESeqDataSet
dds <- DESeqDataSetFromMatrix(
    countData = countdata, 
    colData = coldata,
    design = ~ treatment
    )

dds$treatment <- relevel(dds$treatment, ref="WT")

# quick check of what we have
as.data.frame(colData(dds))

# pre-filtering
smallestGroupSize <- 2
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

dds_norm_prep <- estimateSizeFactors(dds)
dds_norm_counts <- counts(dds_norm_prep, normalized=TRUE)
dds_norm_counts <- as.data.frame(dds_norm_counts)
dds_norm_counts <- rownames_to_column(dds_norm_counts, var = "gene_id")
# write.table(dds_norm_counts, file=paste0("/path/to/normalized_count_matrix_2.csv"), sep=",", quote=F, row.names=FALSE)

dds <- DESeq(dds)

# extract the DE results
save_res <- function(treatment1, treatment0){

    # extract the DE results
    # we need to (and are able to) extract comparison of any two levels of a variable
    # specify three values:  the name of the variable, the name of the level in the numerator, and the name of the level in the denominator
    res <- as.data.frame(results(dds, contrast = c("treatment", treatment1, treatment0)))

    # convert the rownames into the first column
    res <- rownames_to_column(res, var = "gene_id")
    # simplify the gene name
    # original name from gene count matrix: ENSG00000223764.2|LINC02593
    # simplified name: LINC02593
    library(stringr)
    res$gene_id <- str_split_i(res$gene_id, "[|]", 2)
    # replace the NA value in DESeq results
    res$padj[is.na(res$padj)] <- 1
    # sort the data
    res <- res[order(res$padj, res$log2FoldChange, decreasing = c(FALSE, TRUE)), ]

    # mark the significant changes
    res[,"sig"] <- "non-sig"
    res[which(res$log2FoldChange >= 0 & res$padj < 0.05),'sig'] <- 'up'
    res[which(res$log2FoldChange <= 0 & res$padj < 0.05),'sig'] <- 'down'

    res_sig <- subset(res, sig %in% c('up', 'down'))

    # save the DE results to csv
    library(readr)
    write_csv(res, paste0("/path/to/de_", treatment1, "_vs_", treatment0, ".csv"))
    write_csv(res_sig, paste0("/path/to/de_", treatment1, "_vs_", treatment0, "_sig.csv"))
}

treatments <- c("Ythdf1_KO","Ythdf2_KO","Ythdf3_KO","Mettl3_KO","YTHDF_Tri_KO")
for (i in treatments){
    save_res(i, "WT")
}