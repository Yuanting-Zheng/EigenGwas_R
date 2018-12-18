install.packages("BiocManager")
BiocManager::install("ComplexHeatmap", version = "3.8")

biocLite

options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/")

biocLite("ComplexHeatmap")
a

library(ComplexHeatmap)






install.packages("pacman")
pacman::p_load(ComplexHeatmap, circlize)

set.seed(7)

mat <- cbind(rbind(matrix(rnorm(16, -1),4), matrix(rnorm(32, 1), 8)), rbind(matrix(rnorm(24, 1), 4), matrix(rnorm(48, -1), 8)))

mat <- mat[sample(nrow(mat), nrow(mat)), sample(ncol(mat), ncol(mat))]

rownames(mat) <- paste0("R", 1:12)

colnames(mat) <- paste0("C", 1:10)

Heatmap(mat)












