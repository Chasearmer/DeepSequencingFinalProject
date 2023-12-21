# initialize
rm(list=ls()) ; gc();
set.seed(1337)

# load library
library(stringr)
library(readr)

subs <- c("Ythdf1_KO","Ythdf2_KO","Ythdf3_KO","Mettl3_KO","YTHDF_Tri_KO")

for (sub in subs){
# load the data
deg <- read.csv(
    paste0("/path/to/de_", sub, "_vs_WT.csv"), 
    header = TRUE, 
    sep = ","
    )
m6A_all <- read.csv(
    "/path/to/mm10_CL_Tech.tsv",
    header = TRUE,
    sep = "\t"
    )

# filter the m6A_all data  
# ideally we only need m6A sites detected in mESC
# m6A_all$Cell_Line: brain;dorsal root ganglion?(DRG);Skin;ESC
ctoi <- "ESC"

# get the $Cell_Line
cell_lines <- m6A_all$Cell_Line
# for each in $Cell_Line: split the character with sep=";", return a vector
total_rows <- length(cell_lines)
m6A_i_vec <- c()
for (i in 1:total_rows){
    i_row_ct <- cell_lines[i]
    i_row_ct_split <- str_split_1(i_row_ct, ";")
    # judge whether "ESC" is in the vector
    if (ctoi %in% i_row_ct_split){
        # if it is, return the row index
        m6A_i_vec <- append(m6A_i_vec, i)
        }
    }
# only keep the rows in m6A_all whose index is returned
m6A_ESC <- m6A_all[m6A_i_vec,]

# number of m6A sites of each gene  
library(dplyr)
m6A_ESC_number <- m6A_ESC %>% group_by(m6A_ESC$Gene_Name) %>% summarize(count = n())

# create a new column to specify m6A sites number per gene
m6A_ESC_number.list <- as.list(m6A_ESC_number$count)
names(m6A_ESC_number.list) <- m6A_ESC_number$`m6A_ESC$Gene_Name`

# get the $Gene_Name from the filtered m6A data
m6A_ESC.genes <- m6A_ESC$Gene_Name
# filter the deg table with m6A genes in ESC / or just create a new column to specify if it's m6A_ESC_mm or not
deg[which(deg$gene_id %in% m6A_ESC.genes),"m6A_ESC_mm"] <- "m6A"
deg$m6A_ESC_mm[is.na(deg$m6A_ESC_mm)] <- "non-m6A"

deg[which(deg$m6A_ESC_mm == "m6A"), "m6A_ESC_mm_num"] <- unlist(lapply(deg[which(deg$m6A_ESC_mm == "m6A"),]$gene_id, function(x)m6A_ESC_number.list[[x]]))
deg$m6A_ESC_mm_num[is.na(deg$m6A_ESC_mm_num)] <- 0

# create a new column to specify bins
deg[which(deg$m6A_ESC_mm_num == 0), "m6A_ESC_mm_bin"] <- "0 m6A"
deg[which(deg$m6A_ESC_mm_num == 1), "m6A_ESC_mm_bin"] <- "1 m6A"
deg[which(deg$m6A_ESC_mm_num >= 5), "m6A_ESC_mm_bin"] <- "5+ m6A"
deg$m6A_ESC_mm_bin[is.na(deg$m6A_ESC_mm_bin)] <- "2-4 m6A"

# save the dataset
write_csv(deg, paste0("/path/to/de_", sub, "_vs_WT_m6A.csv"))
}