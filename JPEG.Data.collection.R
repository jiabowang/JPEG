
# Based on 565 hapmap file to create numeric genotype data
rm(list=ls())
setwd("/home/jiabowang/data/Lanjuan_Rice/hapmap")

source("/home/jiabowang/Code/GAPIT.library.R")
source("/home/jiabowang/Code/gapit_functions.txt")
library(data.table)
myG1=data.table::fread("hmp565.txt",
                          header = F,
                          na.strings = c("NA", "NaN"),
                          data.table = F)
GAPIT0=GAPIT(G=myG1,PCA.total=3,file.output=F)
myGD=GAPIT0$GD
myGM=GAPIT0$GM
myCV=GAPIT0$PC

fwrite(myGD,"565.GD.txt",row.names=F,quote=FALSE,col.names=T)
fwrite(myGM,"565.GM.txt",row.names=F,quote=FALSE,col.names=T)
fwrite(myCV,"565.CV.txt",row.names=F,quote=FALSE,col.names=T)

setwd("/home/jiabowang/data/Lanjuan_Rice")
myGD=data.table::fread("565.GD.txt",
	header=T)
taxa=as.character(as.matrix(as.data.frame(myGD[,1])))
X=as.matrix(myGD[,-1])
pp=GAPIT.PCA(X=X,taxa=taxa,file.output=TRUE,PCA.total=length(taxa))
pc=pp$PCs
pcs=cbind(as.data.frame(taxa),pc[,-1])
write.csv(pcs,"GAPIT.PCA.csv",quote=F,row.names=F)
system("mv GAPIT.PCA.csv 565.PCA.csv")
system("mv GAPIT.PCA.eigenValue.pdf 565.PCA.eigenValue.pdfs")

setwd("/home/jiabowang/data/Lanjuan_Rice")
myGD=data.table::fread("115.GD.txt",
	header=T)
taxa=as.character(as.matrix(as.data.frame(myGD[,1])))
X=as.matrix(myGD[,-1])
pp=GAPIT.PCA(X=X,taxa=taxa,file.output=TRUE,PCA.total=length(taxa))
system("mv GAPIT.PCA.csv 115.PCA.csv")
system("mv GAPIT.PCA.eigenValue.pdf 115.PCA.eigenValue.pdfs")

# Based on 118 hapmap file to create numeric genotype data

rm(list=ls())
source("/home/jiabowang/Code/GAPIT.library.R")
source("/home/jiabowang/Code/gapit_functions.txt")
setwd("/home/jiabowang/data/Lanjuan_Rice/hapmap")
library(data.table)
myG2=data.table::fread("hmp118.txt",
                          header = F,
                          na.strings = c("NA", "NaN"),
                          data.table = F)
GAPIT0=GAPIT(G=myG2,PCA.total=0,file.output=F)
myGD=GAPIT0$GD
myGM=GAPIT0$GM
myCV=GAPIT0$PC
rm(myG2,GAPIT0)
fwrite(myGD,"118.GD.txt",row.names=F,quote=FALSE,col.names=T)
fwrite(myGM,"118.GM.txt",row.names=F,quote=FALSE,col.names=T)
fwrite(myCV,"118.CV.txt",row.names=F,quote=FALSE,col.names=T)
