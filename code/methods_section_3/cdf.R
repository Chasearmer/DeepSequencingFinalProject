# initialize
rm(list=ls()) ; gc();
set.seed(1337)

library(ggplot2)

paths <- c("Mettl3_KO_vs_WT","Mettl3_KO_vs_WT_prelim","YTHDF_Tri_KO_vs_WT")
path_titles <- as.list(c("Mettl3-KO vs WT","Mettl3-KO vs WT prelim","YTHDF1/2/3-Triple-KO vs WT"))
names(path_titles) <- paths

for (path in paths){

deg <- read.csv(
    paste0("/path/to/de_",path,"_m6A.csv"), 
    header = TRUE, 
    sep = ","
    )

path_title <- path_titles[[path]]

# cumulative distribution plot: m6A
cdf_plot <- function(df,c,title_add){
    p <- ggplot(
        df,
        aes(log2FoldChange, colour = c)
        ) + 
        stat_ecdf(pad = FALSE) + 
        theme_bw() + 
        xlab("log2FoldChange(KO/WT)") + 
        ylab("Cumulative fraction") + 
        xlim(c(-10,10)) + 
        ggtitle(paste0(path_title, " ", title_add)) + 
        theme(
            legend.title=element_blank(), 
            legend.position = c(.15, .70),
            plot.title = element_text(hjust = 0.5, size = 12)
            )
}
p1 <- cdf_plot(deg, deg$m6A_ESC_mm_bin, "(all genes)")
ggsave(paste0("/path/to/de_",path,"/de_cdf_m6A_all.png"), p1, width = 5, height = 3)

# cumulative distribution plot: coding/lncRNAs
deg_main <- deg[which(deg$gene_type %in% c("protein_coding","lncRNA")),]

## coding VS lncRNAs, all
p2 <- cdf_plot(deg_main, deg_main$gene_type, "(all genes)")
ggsave(paste0("/path/to/de_",path,"/de_cdf_codinglnc_all.png"), p2, width = 5, height = 3)

## coding VS lncRNAs, m6A
deg_main_m6A <- deg_main[which(deg_main$m6A_ESC_mm == "m6A"),]

p3 <- cdf_plot(deg_main_m6A, deg_main_m6A$gene_type, "(all genes with at least 1 m6A)")
ggsave(paste0("/path/to/de_",path,"/de_cdf_codinglnc_m6A_all.png"), p3, width = 5, height = 3)

## m6A/non-m6A lncRNA
deg_lnc <- deg[which(deg$gene_type == "lncRNA"),]

p4 <- cdf_plot(deg_lnc, deg_lnc$m6A_ESC_mm_bin, "(lncRNA)")
ggsave(paste0("/path/to/de_",path,"/de_cdf_lnc_all.png"), p4, width = 5, height = 3)

## m6A/non-m6A coding
deg_coding <- deg[which(deg$gene_type == "protein_coding"),]

p5 <- cdf_plot(deg_coding, deg_coding$m6A_ESC_mm_bin, "(protein coding)")
ggsave(paste0("/path/to/de_",path,"/de_cdf_coding_all.png"), p5, width = 5, height = 3)

}