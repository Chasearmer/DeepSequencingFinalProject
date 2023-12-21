# initialize
rm(list=ls()) ; gc();
set.seed(1337)

library(dplyr)
library(readr)

treatments <- c(
    "Mettl3_KO",
    "Ythdf1_KO",
    "Ythdf2_KO",
    "Ythdf3_KO",
    "YTHDF_Tri_KO"
)

grid <- data.frame(
    treatment = character(0),
    all = integer(0),
    all_padj = integer(0),
    all_padj_m6A = integer(0),
    coding = integer(0),
    coding_padj = integer(0),
    coding_padj_m6A = integer(0),
    lnc = integer(0),
    lnc_padj = integer(0),
    lnc_padj_m6A = integer(0)
)

for (treatment in treatments){

deg <- read.csv(
    paste0("/path/to/de_", treatment, "_vs_WT_m6A.csv"), 
    header = TRUE, 
    sep = ","
    )

new_row <- data.frame(
    treatment = treatment,
    all = nrow(deg),
    all_padj = nrow(subset(deg, padj < 0.05)),
    all_padj_m6A = nrow(subset(deg, (padj < 0.05) & (m6A_ESC_mm == "m6A"))),
    coding = nrow(subset(deg, gene_type == "protein_coding")),
    coding_padj = nrow(subset(deg, (gene_type == "protein_coding") & (padj < 0.05))),
    coding_padj_m6A = nrow(subset(deg, (gene_type == "protein_coding") & (padj < 0.05) & (m6A_ESC_mm == "m6A"))),
    lnc = nrow(subset(deg, gene_type == "lncRNA")),
    lnc_padj = nrow(subset(deg, (gene_type == "lncRNA") & (padj < 0.05))),
    lnc_padj_m6A = nrow(subset(deg, (gene_type == "lncRNA") & (padj < 0.05) & (m6A_ESC_mm == "m6A")))
)

grid <- rbind(grid, new_row)

}

write_csv(grid,"/path/to/impacted_gene_gird.csv")