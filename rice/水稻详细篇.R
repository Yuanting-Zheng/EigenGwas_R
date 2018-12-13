# 此脚本针对水稻自交系
setwd("rice")
plink='/Users/mac/tool/plink'
#  输入你的项目名称 针对bed bim fam 
name = "404" 

# 检查一下bim 文件
check_bim = paste0("head ", name,".bim") 
system(check_bim)

# 使用plink 的 freq 检查一下
check_freq = paste0(plink," --bfile ", name," --freq --out freq_analysis ") 
system(check_freq)

frq <- read.table("freq_analysis.frq",header = T) #读取数据
maf_range = range(frq$MAF)    # 求出最大最小值
maf_name = paste0("MAF ",maf_range[1],"-",maf_range[2])
tiff(filename = "frq.tiff",
     width = 2500, height = 2500, units = "px",res = 300,compression="lzw")
hist(frq$MAF, main = maf_name, breaks = 50, xlim=c(0,0.5))  # 画个频率分布直方图 看看是否别人前期处理过了
dev.off() # 保存图片


#制作grm 矩阵
make_grm = paste0(plink," --bfile ", name," --make-grm-gz --allow-extra-chr") 
system(make_grm)
unzip_grm = "gunzip plink.grm.gz"
system(unzip_grm)

#Total genotyping rate is 0.955783.
#404388 variants and 3024 people pass filters and QC.

#> system(check_eigenval)
#471.149
#242.477
#146.372
#113.96
#104.472



#读取 grm 矩阵 并计算 Ne Me 画一幅图  按照inbred的做 修改了outbred的参数除以二变为除以1 
grm <- read.table("plink.grm",header = F) #读取数据
grm_trim = grm[which(grm[,1]!=grm[,2]),]  # 已经提取了左右不相等的个体
grm_trim$data_by2 = grm_trim[,4]/2 # 增加一个第五列，值是第四列的一半  这里inbred我除了2
grm_inbred <- grm_trim[,c(1,2,3,5)] # 做一个inbred的举证

Ne = -1/mean((grm_inbred[,4]))# inbred line 要除以2的矩阵 有效个体数
Me = 1/var((grm_inbred[,4]))# inbred line 要除以2的矩阵 有效标记数
Ne_str <- as.character(Ne)
Me_str <- as.character(Me)

Ne_text = paste("Ne is ",Ne_str," ")
Me_text = paste("Me is ",Me_str," ")

tiff(filename = "grm.tiff",
     width = 2500, height = 2500, units = "px",res = 300,compression="lzw")
hist(grm_inbred$data_by2,breaks = 25,xlab = "GRM score",main = paste("GRM distribution",Ne_text,Me_text)) # 对第四列的值进行统计 画图
dev.off() # 保存图片

# 计算pca 与 eigenval 我这里计算5个为标准

pca_number = 5 # 这里设置计算的pca 数量

make_pca = paste0(plink," --bfile ", name, " ", "--pca ", pca_number, " --out pca --allow-extra-chr") 
system(make_pca)

check_eigenval = "head pca.eigenval"
system(check_eigenval)

# 以5个主成分为y值，做单标记回归
do_Regression = paste0(plink," --bfile ", name, " --linear --pheno pca.eigenvec --all-pheno --allow-no-sex --allow-extra-chr") 
system(do_Regression)


# 检查一下 取出P1的结果
P1 <- read.table("plink.P1.assoc.linear",header = TRUE, stringsAsFactors=F) #读取数据

pdata=median(P1$P) # 取出P值的中位数
pdata
which(is.na(P1$P)) # 没有NA

#计算GC值（每个主成分的）画出比较的图
library(qqman)
library(vcd)
library(grid)
P1 <- read.table("plink.P1.assoc.linear",header = T) #读取数据
pdata=median(P1$P) # 取出P值的中位数
chi=qchisq(pdata, 1, lower.tail = F)  # 按照这个P值的中位数去 chi 的值
gc1=chi/qchisq(0.5, 1, lower.tail = F)  #  按照这个分布的chi值 除以 理论为P = 0.5的qchisq 这个与理论的比值就是 GC
print(gc1) 

P2 <- read.table("plink.P2.assoc.linear",header = T) #读取数据
pdata=median(P2$P) # 取出P值的中位数
chi=qchisq(pdata, 1, lower.tail = F)  # 按照这个P值的中位数去 chi 的值
gc2=chi/qchisq(0.5, 1, lower.tail = F)  #  按照这个分布的chi值 除以 理论为P = 0.5的qchisq 这个与理论的比值就是 GC
print(gc2) 

P3 <- read.table("plink.P3.assoc.linear",header = T) #读取数据
pdata=median(P3$P) # 取出P值的中位数
chi=qchisq(pdata, 1, lower.tail = F)  # 按照这个P值的中位数去 chi 的值
gc3=chi/qchisq(0.5, 1, lower.tail = F)  #  按照这个分布的chi值 除以 理论为P = 0.5的qchisq 这个与理论的比值就是 GC
print(gc3) 

P4 <- read.table("plink.P4.assoc.linear",header = T) #读取数据
pdata=median(P4$P) # 取出P值的中位数
chi=qchisq(pdata, 1, lower.tail = F)  # 按照这个P值的中位数去 chi 的值
gc4=chi/qchisq(0.5, 1, lower.tail = F)  #  按照这个分布的chi值 除以 理论为P = 0.5的qchisq 这个与理论的比值就是 GC
print(gc4) 

P5 <- read.table("plink.P5.assoc.linear",header = T) #读取数据
pdata=median(P5$P) # 取出P值的中位数
chi=qchisq(pdata, 1, lower.tail = F)  # 按照这个P值的中位数去 chi 的值
gc5=chi/qchisq(0.5, 1, lower.tail = F)  #  按照这个分布的chi值 除以 理论为P = 0.5的qchisq 这个与理论的比值就是 GC
print(gc5) 


gc_all <- c(gc1,gc2,gc3,gc4,gc5)
eigenval_gc <- read.table("pca.eigenval",header = F)/2 #inbred
eigenval_gc$V2 = gc_all
names(eigenval_gc) <- c("eigenval","gc")
par(mfrow=c(1,2))

tiff(filename = "VS.tiff",
     width = 2500, height = 2500, units = "px",res = 300,compression="lzw")
barplot(t(as.matrix(eigenval_gc)),main = "Mean VS Middle",names.arg = rownames(eigenval_gc), xlab = "PCA",ylab = "value",col = c("red","yellow"),legend=colnames(eigenval_gc),beside = TRUE)
dev.off() # 保存图片
eigenval_gc


#> eigenval_gc
#eigenval         gc
#1 235.5745 164.047491
#2 121.2385  39.477044
#3  73.1860   4.694614
#4  56.9800  17.504229
#5  52.2360  19.436042


# 开始画曼哈顿图，画曼哈顿图之前，需要检查P值的范围，这里只以第一个主成分为例，介绍矫正前和用plink矫正后的曼哈顿图情况
plink_adjust = paste0(plink," --bfile ", name, " --linear --pheno pca.eigenvec --all-pheno --adjust --allow-no-sex --allow-extra-chr") 
system(plink_adjust)

library(qqman) # 载入包
P1_for_man <- read.table("plink.P1.assoc.linear",header = TRUE) #读取数据
P1_adjust_for_man <- read.table("plink.P1.assoc.linear.adjusted",header = TRUE) #读取数据

head(P1_for_man)
head(P1_adjust_for_man)
# 检查是否有P值等于O的点
which(P1_for_man$P==0)
which(P1_adjust_for_man$GC=='Inf')
P1_for_man[14677,]
P1_for_man[3275343,]
P1_adjust_for_man[P1_adjust_for_man$SNP=="14677",]
P1_adjust_for_man[2,]

# 先把0点去掉画图  
rm_zero = P1_for_man[-which(P1_for_man$P==0),] 

head(rm_zero)
which(rm_zero$P==0)

man_P1_final <- rm_zero[,c(1,2,3,9)]
colnames(man_P1_final) <- c("CHR","SNP","BP","P")


rm_zero = P1_for_man[-which(P1_for_man$P==0),] 
rm_inf = P1_adjust_for_man[-which(P1_adjust_for_man$GC=="Inf"),] 
rm_combine = merge(man_P1_final,P1_adjust_for_man,by.x = 'SNP',by.y = 'SNP')   # 这里做了修改

head(rm_combine) 
man_P1 <- data.frame(rm_combine$CHR.x,rm_combine$SNP,rm_combine$BP,rm_combine$P)
man_P1_adjust <- data.frame(rm_combine$CHR.x,rm_combine$SNP,rm_combine$BP,rm_combine$GC)
man_P1 <- na.omit(man_P1) # 删除含有NA的整行
man_P1_adjust <- na.omit(man_P1_adjust) # 删除含有NA的整行
colnames(man_P1) <- c("CHR","SNP","BP","P")
colnames(man_P1_adjust) <- c("CHR","SNP","BP","P")


head(man_P1_adjust)
# 矫正前画图
par(cex=0.8) #设置点的大小
tiff(filename = "Manhadun_PC1.tiff",
     width = 2500, height = 2500, units = "px",res = 300,compression="lzw")

manhattan(man_P1_final,main="Manhattan Plot PC1",col = c("#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494","#66C2A5","#FC8D62"),suggestiveline = FALSE,annotatePval = 0.01) #suggestiveline = FALSE 更加显著
dev.off() # 保存图片

tiff(filename = "QQ_PC1.tiff",
     width = 2500, height = 2500, units = "px",res = 300,compression="lzw")
qq(man_P1_final$P, main = "Q-Q plot of GWAS p-values PC1", xlim = c(0, 7), ylim = c(0, 12), pch = 18, col = "blue4", cex = 1.5, las = 1)
dev.off() # 保存图片

#矫正后画图
tiff(filename = "Manhadun_PC1_adjusted_plink_GC.tiff",
     width = 2500, height = 2500, units = "px",res = 300,compression="lzw")

manhattan(man_P1_adjust,main="Manhattan Plot PC1_adjusted",col = c("#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494","#66C2A5","#FC8D62"),suggestiveline = FALSE,annotatePval = 0.01) #suggestiveline = FALSE 更加显著
dev.off() # 保存图片

tiff(filename = "QQ_PC1_adjusted_plink_GC.tiff",
     width = 2500, height = 2500, units = "px",res = 300,compression="lzw")
qq(man_P1_adjust$P, main = "Q-Q plot of GWAS p-values PC1_adjusted", xlim = c(0, 7), ylim = c(0, 12), pch = 18, col = "blue4", cex = 1.5, las = 1)
dev.off() # 保存图片





