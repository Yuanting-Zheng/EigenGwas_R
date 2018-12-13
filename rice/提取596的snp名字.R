setwd("rice")
plink='/Users/mac/tool/plink'
data <- read.table("596.bim",header = F)
head(data)
snp_name <- data$V2
head(snp_name)
write.table(snp_name,"snp_name.txt",col.names = F,row.names = F,quote = FALSE) # 中文的要去除引号
