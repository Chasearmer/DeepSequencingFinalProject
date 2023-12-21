# initialize
rm(list=ls()) ; gc();
set.seed(1337)

# load library
library("rtracklayer")

subs <- c("Ythdf1_KO","Ythdf2_KO","Ythdf3_KO","Mettl3_KO","YTHDF_Tri_KO")

for (sub in subs){
# load the data
deg <- read.csv(
    paste0("/path/to/de_", sub, "_vs_WT.csv"), 
    header = TRUE, 
    sep = ","
    )

anno_gtf <- rtracklayer::import('/path/to/gencode.vM23.annotation.gtf')
anno_gtf <- as.data.frame(anno_gtf)

# retrieve the gene_name:gene_type relationship
# keep the genes, not the transcripts
anno_gtf <- anno_gtf[which(anno_gtf$type == "gene"),c(11,12)]
name_type <- as.list(anno_gtf$gene_type)
names(name_type) <- anno_gtf$gene_name

# create a new column in deg to specify gene_type according to gene_id
deg$gene_type <- unlist(lapply(deg$gene_id, function(x)name_type[[x]]))
head(deg)

# save the dataset
library(readr)
write_csv(deg, paste0("/path/to/de_", sub, "_vs_WT.csv"))
}