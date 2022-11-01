library("seqinr")
library("phangorn")
library("msa")
library("ggplot2")
library("reshape2")

# set working directory
WORK_DIR <- "E:/Kuliah/Semester 7/Bioinformatik & AI/Week 9/Project UTS/Phylogenetic"
setwd(WORK_DIR)

# reading amino acid data
AA_sequence <- readAAStringSet("COX1 Human vs Quadruped Animals Dataset (without species name).fasta")
# sequence <- system.file("DATASET/Project UTS", "COX1 Human vs Quadruped Animals Dataset (without species name).fasta", package="msa")
# AA_sequence <- readAAStringSet(sequence)

# compute multiple sequence alignment score and the matrix distance
msa_result <- msa(AA_sequence, "ClustalW")
msa_result <- msaConvert(msa_result, "seqinr::alignment")
matrix_dist <- dist.alignment(msa_result)

# UPGMA clustering
upgma_cluster <- upgma(matrix_dist)
plot(upgma_cluster, 
     main="UPGMA Cluster: Human vs Quadruped Animals", 
     sub="Phylogenetic Analysis")

# convert matrix distance to matrix data
matrix_data <- as.matrix(matrix_dist)
# convert matrix data to data frame
matrix_df <- data.frame(matrix_data)
# create matrix correlation rounded to 4 integers
matrix_cor <- round(cor(matrix_df), 4)
# melt the matrix to feed it to ggplot
melted_matrix <- melt(matrix_cor)
# heatmap visualization
ggplot(melted_matrix, aes(x=Var1, y=Var2, fill=value)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
             axis.title.y = element_blank(),
             axis.title.x = element_blank()
          #    legend.title = element_blank(),
         ) + 
    geom_tile()


