# initialize
rm(list=ls()) ; gc();
set.seed(1337)

library(ggplot2)
library(plyr)

paths <- c("Mettl3_KO_vs_WT","YTHDF_Tri_KO_vs_WT")
path_titles <- as.list(c("Mettl3-KO vs WT","YTHDF1/2/3-Triple-KO vs WT"))
names(path_titles) <- paths

for (path in paths){

deg <- read.csv(
    paste0("/path/to/de_",path,"_m6A.csv"), 
    header = TRUE, 
    sep = ","
    )

path_title <- path_titles[[path]]

density_plot <- function(df,c,cname,title_add){
    mu <- ddply(df, cname, summarise, logFC.mean=mean(log2FoldChange))
    colnames(mu) <- c("group","logFC.mean")

    p <- ggplot(
        df, 
        aes(x=log2FoldChange, color=c)
        ) + 
        geom_density() + 
        #geom_vline(data=mu, aes(xintercept=logFC.mean, color=group),linetype="dashed") + 
        geom_vline(aes(xintercept=0), color="grey", linetype="dashed") + 
        theme_classic() + 
        xlim(c(-5,5)) + 
        ylim(c(0,1)) + 
        xlab("log2FoldChange(KO/WT)") + 
        ylab("Density") + 
        ggtitle(paste0(path_title, " ", title_add)) + 
        theme(
            legend.title=element_blank(), 
            legend.position = c(.15, .70),
            plot.title = element_text(hjust = 0.5, size = 12)
            )

    return(p)
}

deg_main <- deg[which(deg$gene_type %in% c("protein_coding","lncRNA")),]
## coding VS lncRNAs, m6A
deg_main_m6A <- deg_main[which(deg_main$m6A_ESC_mm == "m6A"),]
p3 <- density_plot(deg_main_m6A, deg_main_m6A$gene_type, "gene_type", "(all genes with at least 1 m6A)")
ggsave(paste0("/path/to/de_",path,"/de_density_codinglnc_m6A_all.png"), p3, width = 5, height = 3)

## m6A/non-m6A lncRNA
deg_lnc <- deg[which(deg$gene_type == "lncRNA"),]
p4 <- density_plot(deg_lnc, deg_lnc$m6A_ESC_mm_bin, "m6A_ESC_mm_bin", "(lncRNA)")
ggsave(paste0("/path/to/de_",path,"/de_density_lnc_all.png"), p4, width = 5, height = 3)

## m6A/non-m6A coding
deg_coding <- deg[which(deg$gene_type == "protein_coding"),]
p5 <- density_plot(deg_coding, deg_coding$m6A_ESC_mm_bin, "m6A_ESC_mm_bin", "(protein coding)")
ggsave(paste0("/path/to/de_",path,"/de_density_coding_all.png"), p5, width = 5, height = 3)

}