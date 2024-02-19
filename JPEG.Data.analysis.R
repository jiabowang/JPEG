
# Based on 565 hapmap file to create numeric genotype data
rm(list=ls())
setwd("/home/jiabowang/data/Lanjuan_Rice/hapmap")

source("/home/jiabowang/Code/GAPIT.library.R")
source("/home/jiabowang/Code/gapit_functions.txt")
library(data.table)
# read hapmap genotype file
myG1=data.table::fread("hmp565.txt",
                          header = F,
                          na.strings = c("NA", "NaN"),
                          data.table = F)
GAPIT0=GAPIT(G=myG1,PCA.total=3,file.output=F)
myGD=GAPIT0$GD
myGM=GAPIT0$GM
myCV=GAPIT0$PC
# output numeric genotype files
fwrite(myGD,"565.GD.txt",row.names=F,quote=FALSE,col.names=T)
fwrite(myGM,"565.GM.txt",row.names=F,quote=FALSE,col.names=T)
fwrite(myCV,"565.CV.txt",row.names=F,quote=FALSE,col.names=T)

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


################## Joint Analysis

rm(list=ls())

source("/home/jiabowang/Code/GAPIT.library.R")
source("/home/jiabowang/Code/gapit_functions.txt")
library(data.table)



myGD1=data.table::fread("565.GD.txt",
	header=T)
myGM1=read.table("565.GM.txt",sep=",",head=T)

myGD2=data.table::fread("118.GD.txt",
	header=T)
myGM2=read.table("118.GM.txt",sep=",",head=T)


index=as.character(myGM2[,1])%in%as.character(myGM1[,1])
GD=as.matrix(myGD2)
newGD2=GD[,c(TRUE,index)]


rm(GD,myGD1,myGD2)
fwrite(newGD2,"118.GD.txt",row.names=F,quote=FALSE,col.names=T)
fwrite(myGM1,"118.GM.txt",row.names=F,quote=FALSE,col.names=T)



#Inbred

#inbred phenotype
setwd("/home/jiabowang/data/Lanjuan_Rice/Phenotype")
myY1=read.table("V.traits.txt",head=T)
myY2=read.table("GC.traits.txt",head=T)

myGD2=data.table::fread("118.GD.txt",header=T)
myGM2=read.table("118.GM.txt",sep=",",head=T)
library(data.table)
taxav=paste("V_",as.character(as.matrix(as.data.frame(myGD2[,1]))),sep="")
taxagc=paste("GC_",as.character(as.matrix(as.data.frame(myGD2[,1]))),sep="")

X=as.matrix(myGD2[,-1])

taxall=append(taxav,taxagc)
y0=myY1[,-1]
y1=myY2[,-1]
y00=apply(y0,2,function(one) (one-mean(one,na.rm=T))/sqrt(var(one,na.rm=T)))
y11=apply(y1,2,function(one) (one-mean(one,na.rm=T))/sqrt(var(one,na.rm=T)))

y3=rbind(y00,y11)
inbredY=cbind(as.data.frame(taxall),y3)
write.table(inbredY,"inbred.Y.txt",quote=F,row.names=F)


#hybrid phenotype
setwd("/home/jiabowang/data/Lanjuan_Rice/Phenotype")
myY1=read.table("TC.traits.txt",head=T)
myY2=read.table("HMP.traits.txt",head=T)

taxav=paste("TC_",as.character(myY1[,1]),sep="")
taxagc=paste("HMP_",as.character(myY2[,1]),sep="")
taxall=append(taxav,taxagc)
y0=myY1[,-1]
y1=myY2[,-1]
y00=apply(y0,2,function(one) (one-mean(one,na.rm=T))/sqrt(var(one,na.rm=T)))
y11=apply(y1,2,function(one) (one-mean(one,na.rm=T))/sqrt(var(one,na.rm=T)))

y3=rbind(y00,y11)
hybridY=cbind(as.data.frame(taxall),y3)
write.table(hybridY,"hybrid.Y.txt",quote=F,row.names=F)

#combind phenotype
setwd("/home/jiabowang/data/Lanjuan_Rice/Phenotype")
myY1=read.table("inbred.Y.txt",head=T)
myY2=read.table("hybrid.Y.txt",head=T)

taxav=as.character(myY1[,1])
taxagc=as.character(myY2[,1])
taxall=append(taxav,taxagc)
y0=myY1[,-1]
y1=myY2[,-1]
y00=apply(y0,2,function(one) (one-mean(one,na.rm=T))/sqrt(var(one,na.rm=T)))
y11=apply(y1,2,function(one) (one-mean(one,na.rm=T))/sqrt(var(one,na.rm=T)))

y3=rbind(y00,y11)
combinY=cbind(as.data.frame(taxall),y3)
write.table(combinY,"combin.Y.txt",quote=F,row.names=F)

# inbred genotype
rm(list=ls())
setwd("/home/jiabowang/data/Lanjuan_Rice/hapmap")

myGD2=data.table::fread("118.GD.txt",header=T)
myGM2=read.table("118.GM.txt",sep=",",head=T)
library(data.table)
taxav=paste("V_",as.character(as.matrix(as.data.frame(myGD2[,1]))),sep="")
taxagc=paste("GC_",as.character(as.matrix(as.data.frame(myGD2[,1]))),sep="")

X=as.matrix(myGD2[,-1])

taxa=c(taxav,taxagc)
Xin=rbind(X,X)
GDin=cbind(as.data.frame(taxa),Xin)
str=c(rep(1,length(taxav)),rep(0,length(taxav)))
CVin=cbind(as.data.frame(taxa),str)

fwrite(myGM2,"inbred.GM.txt",row.names=F,quote=FALSE,col.names=T)
fwrite(CVin,"inbred.CV.txt",row.names=F,quote=FALSE,col.names=T)
rm(X,Xin)
rm(myGD2)
fwrite(GDin,"inbred.GD.txt",row.names=F,quote=FALSE,col.names=T)

# hybrid genotype
rm(list=ls())
setwd("/home/jiabowang/data/Lanjuan_Rice/hapmap")

myGD1=data.table::fread("565.GD.txt",header=T)
myGM1=read.table("565.GM.txt",sep=",",head=T)
library(data.table)
taxav=paste("TC_",as.character(as.matrix(as.data.frame(myGD1[,1]))),sep="")
taxagc=paste("HMP_",as.character(as.matrix(as.data.frame(myGD1[,1]))),sep="")
taxa=c(taxav,taxagc)

SNP=as.character(myGM1[,1])
chr=as.numeric(myGM1[,2])
SNPin=c(SNP,paste(SNP,"_2",sep=""))
chrin=c(chr,chr+max(chr))
Posin=c(as.numeric(myGM1[,3]),as.numeric(myGM1[,3]))


str=c(rep(1,length(taxav)),rep(0,length(taxav)))
CVhy=cbind(as.data.frame(taxa),str)
GMhy=cbind(as.data.frame(SNPin),chrin,Posin)


rm(X,myGD1)
rm(Xhy1,Xhy2,Xhy)
fwrite(GMhy,"hybrid.GM.txt",row.names=F,quote=FALSE,col.names=T)
fwrite(CVhy,"hybrid.CV.txt",row.names=F,quote=FALSE,col.names=T)
rm(GMhy,CVhy)

rm(list=ls())
setwd("/home/jiabowang/data/Lanjuan_Rice/hapmap")

myGD1=data.table::fread("565.GD.txt",header=T)
library(data.table)
taxav=paste("TC_",as.character(as.matrix(as.data.frame(myGD1[,1]))),sep="")
taxagc=paste("HMP_",as.character(as.matrix(as.data.frame(myGD1[,1]))),sep="")
taxa=c(taxav,taxagc)
X=as.matrix(myGD1[,-1])
Xhy1=cbind(X,abs(abs(X-1)-1))
Xhy2=cbind(X*0,abs(abs(X-1)-1))
Xhy=rbind(Xhy1,Xhy2)
GDhy=cbind(as.data.frame(taxa),Xhy)
blink_GD=t(GDhy[,-1])
setwd("/home/jiabowang/data/Lanjuan_Rice/Phenotype")
fwrite(blink_GD,"hybrid.dat",quote=F,sep=" ",col.names=F,row.names=F)

# Combin genotype
rm(list=ls())
setwd("/home/jiabowang/data/Lanjuan_Rice/hapmap")

myGM1=read.table("inbred.GM.txt",sep=",",head=T)
myGM2=read.table("hybrid.GM.txt",sep=",",head=T)
myCV1=read.table("inbred.CV.txt",sep=",",head=T)
setwd("/home/jiabowang/data/Lanjuan_Rice/Phenotype")
myCV2=read.table("hybrid.cov",head=T)



library(data.table)
taxain=as.character(as.matrix(as.data.frame(myCV1[,1])))
taxahy=as.character(as.matrix(as.data.frame(myCV2[,1])))
taxaall=c(taxain,taxahy)
cv0=matrix(0,nrow(myCV1),4)
cv0[grep("V116",taxain),1]=1
cv0[grep("V117",taxain),2]=1
cv0[grep("V118",taxain),3]=1
cv0[grep("V119",taxain),4]=1

cv1=cbind(myCV1,cv0)
colnames(cv1)=colnames(myCV2)
blink_CV=rbind(cv1,myCV2)
blink_GM=myGM2
myY1=read.table("inbred.Y.txt",head=T)
myY2=read.table("hybrid.Y.txt",head=T)
myY1[is.na(myY1)]="NaN"
myY2[is.na(myY2)]="NaN"
colnames(myY1)[1]="taxa"
colnames(myY1)[-1]=paste("inbred_",colnames(myY1)[-1],sep="")
colnames(myY2)[-1]=paste("hybrid_",colnames(myY2)[-1],sep="")
fwrite(myY1,"inbred.txt",sep=" ",quote=F,col.names=T,row.names=F)
fwrite(myY2,"hybrid.txt",sep=" ",quote=F,col.names=T,row.names=F)




blink_Y=rbind(myY1,myY2)
blink_Y[is.na(blink_Y)]="NaN"
colnames(blink_Y)[1]="taxa"
colnames(blink_Y)[-1]=paste("combin_",colnames(blink_Y)[-1],sep="")
library(data.table)


fwrite(blink_GM,"combin.map",sep=" ",quote=F,col.names=T,row.names=F)
fwrite(blink_Y,"combin.txt",sep=" ",quote=F,col.names=T,row.names=F)
fwrite(blink_CV,"combin.cov",sep=" ",quote=F,col.names=T,row.names=F)


rm(list=ls())
setwd("/home/jiabowang/data/Lanjuan_Rice/hapmap")

myGD1=data.table::fread("inbred.GD.txt",header=T)
myGD2=data.table::fread("/home/jiabowang/data/Lanjuan_Rice/Phenotype/hybrid.dat",header=F)
X=as.matrix(myGD1[,-1])
Xhy=as.matrix(myGD2)
Xcom1=cbind(X,X*0)
Xcom=cbind(t(Xcom1),Xhy)
rm(myGD1,myGD2,X,Xhy,Xcom1)
setwd("/home/jiabowang/data/Lanjuan_Rice/Phenotype")
library(data.table)
fwrite(Xcom,"combin.dat",row.names=F,sep=" ",quote=FALSE,col.names=F)

### Run Blink-C with each population

list=`seq -s ' ' -w 1 12` ;for i in $list;do ./blink_linux --gwas --file inbred --numeric --trait $i;done;
list=`seq -s ' ' -w 1 12` ;for i in $list;do ./blink_linux --gwas --file hybrid --numeric --trait $i;done;
list=`seq -s ' ' -w 1 12` ;for i in $list;do ./blink_linux --gwas --file combin --numeric --trait $i;done;


########### Estimated H2 for TC population with BGLR


setwd("/home/jiabowang/data/Lanjuan_Rice/hapmap")
set.seed(198521)
GD0=data.table::fread("565.GD.txt",header=T)
GM0=read.table("565.GM.txt")
Y0=read.table("/home/jiabowang/data/Lanjuan_Rice/Phenotype/TC.traits.txt",head=T)

X=as.matrix(GD0[,-1])
XA=X
XD=abs(abs(X-1)-1)

library(BGLR)

source("/home/jiabowang/Code/GAPIT.library.R")
source("/home/jiabowang/Code/gapit_functions.txt")

KA=GAPIT.kinship.VanRaden(snps=XA)
KD=GAPIT.kinship.VanRaden(snps=XD)

taxa.g=as.character(as.matrix(GD0[,1]))
taxa.y=as.character(Y0[,1])
index=match(taxa.g,taxa.y)


ETA<-list(
list(K=KA,model='RKHS'),
list(K=KD,model='RKHS'))



ni=450000
bi=400000
h2=NULL
Var.all=NULL
for(i in 2:ncol(Y0))
{
y=Y0[,i]

fm<-BGLR(y=y,ETA=ETA, nIter=ni, burnIn=bi,saveAt='LZ.rice_')
v1=as.numeric(fm$ETA[[1]]$varU)
v2=as.numeric(fm$ETA[[2]]$varU)
ve=as.numeric(fm$varE)
h2.0=c(v1/(v1+v2+ve),v2/(v1+v2+ve))
Var.0=c(v1,v2,ve)
h2=rbind(h2,h2.0)
Var.all=rbind(Var.all,Var.0)
}
h2=cbind(colnames(Y0)[-1],h2)
colnames(h2)=c("Traits","H2.ADD","H2.DOM")
rownames(h2)=colnames(Y0)[-1]
colnames(Var.all)=c("Var.Add","Var.Dom","Var.res")
rownames(Var.all)=colnames(Y0)[-1]
write.table(h2,"H2.ADD.DOM.TC.txt",quote=F,row.names=F)
write.table(Var.all,"Var.ADD.DOM.TC.txt",quote=F,row.names=F)




############ Estimaed H2 with BGLR using A2D convert additve X to dominant X
A2D<-function(one)
{
    fre=sum(one)/(2*length(one))
    one2=one
    one2[one2==1]=2*round(fre,2)
    one2[one2==2]=4*round(fre,2)-2
    return(one2)
}

setwd("/home/jiabowang/data/Lanjuan_Rice/hapmap")
set.seed(198521)
GD0=data.table::fread("565.GD.txt",header=T)
GM0=read.table("565.GM.txt")
Y0=read.table("/home/jiabowang/data/Lanjuan_Rice/Phenotype/TC.traits.txt",head=T)

X=as.matrix(GD0[,-1])
XA=X
XD=apply(X,2,A2D)
library(BGLR)

source("/home/jiabowang/Code/GAPIT.library.R")
source("/home/jiabowang/Code/gapit_functions.txt")

KA=GAPIT.kinship.VanRaden(snps=XA)
KD=GAPIT.kinship.VanRaden(snps=XD)
n=nrow(KD)

ni=450000
bi=400000

h2=NULL
Var.all=NULL
for(i in 2:ncol(Y0))
{
y=Y0[,i]

fm<-BGLR(y=y,ETA=ETA, nIter=ni, burnIn=bi,saveAt='LZ.A2D_')

v1=as.numeric(fm$ETA[[1]]$varU)
v2=as.numeric(fm$ETA[[2]]$varU)
ve=as.numeric(fm$varE)
h2.0=c(v1/(v1+v2+ve),v2/(v1+v2+ve))
Var.0=c(v1,v2,ve)
h2=rbind(h2,h2.0)
Var.all=rbind(Var.all,Var.0)
}

h2=cbind(colnames(Y0)[-1],h2)
colnames(h2)=c("Traits","H2.Add.a2d","H2.Dom.a2d")
rownames(h2)=colnames(Y0)[-1]
colnames(Var.all)=c("Var.Add.a2d","Var.Dom.a2d","Var.res")
rownames(Var.all)=colnames(Y0)[-1]
write.table(h2,"H2.ADD.DOM.TC.a2d.txt",quote=F,row.names=F)
write.table(Var.all,"Var.ADD.DOM.TC.a2d.txt",quote=F,row.names=F)



########### Estimated H2 for V population with BGLR


setwd("/home/jiabowang/data/Lanjuan_Rice/hapmap")

GD0=data.table::fread("118.GD.txt",header=T)
GM0=read.table("118.GM.txt")
Y0=read.table("/home/jiabowang/data/Lanjuan_Rice/Phenotype/V.traits.txt",head=T)

X=as.matrix(GD0[,-1])

library(BGLR)

source("/home/jiabowang/Code/GAPIT.library.R")
source("/home/jiabowang/Code/gapit_functions.txt")

K=GAPIT.kinship.VanRaden(snps=X)

taxa.g=as.character(as.matrix(GD0[,1]))
taxa.y=as.character(Y0[,1])
index=match(taxa.g,taxa.y)


ETA<-list(
list(K=K,model='RKHS'))



ni=22000
bi=20000
h2=NULL
for(i in 2:ncol(Y0))
{
y=Y0[,i]

fm<-BGLR(y=y,ETA=ETA, nIter=ni, burnIn=bi,saveAt='LZ.rice_')
v1=as.numeric(fm$ETA[[1]]$varU)
# v2=as.numeric(fm$ETA[[2]]$varU)
ve=as.numeric(fm$varE)
h2.0=c(v1/(v1+ve))
h2=rbind(h2,h2.0)
}
h2=cbind(colnames(Y0)[-1],h2)
colnames(h2)=c("Traits","H2")
rownames(h2)=colnames(Y0)[-1]
write.table(h2,"H2.V.txt",quote=F,row.names=F)

################ GS using BLINK with parents and half offsprints to predict remaining

rm(list=ls())
set.seed(99163)
setwd("/home/jiabowang/data/Lanjuan_Rice/Phenotype")
myY=read.table("combin.Y.txt",head=T)
com.cv=read.table("combin.cov",head=T)
traits=as.character(colnames(myY)[-1])
cvindex=com.cv[,3:6]
cvindex[,5]=apply(cvindex,1,function(one) abs(1-sum(one)))
cvindex[1:236,1:5]=0
cv2=cbind(com.cv[,1:2],cvindex)
nrep=30
nfold=5
library(data.table)
setwd("/home/jiabowang/data/Lanjuan_Rice/GS/GS.5P.half")

y5=NULL
for(i in 1:12)
{
        y0=myY[,c(1,i+1)]
        y5=NULL
        for(j in 1:nrep)
        {
        # AA=sample(237:801,282)
        sets=sample(cut(1:565,nfold,labels=FALSE),565)
        y4=NULL
        for(f in 1:nfold)
        {
        AA=sets==f
    	y_cv=c(rep(FALSE,236),AA,AA)# NA index for each rep
    	y=y0
    	y[y_cv,2]=NA
    	y3=y
    	y3[is.na(y)]="NaN"
    	y4=cbind(y4,y3[,2])
        }
        y5=cbind(y5,y4)

        }
        y6=cbind(as.character(y0[,1]),y5)
        colnames(y6)=c("Taxa",paste(traits[i],1:(nrep*nfold),sep=""))
        fwrite(y6,paste("test.T.",i,".txt",sep=""),sep=" ",quote=F,col.names=T,row.names=F)
}


## Run Blink-C in linux system

cp test.T.1.txt combin2.txt #done 1:12
list=`seq -s ' ' -w 1 150` ;for i in $list;do ./blink_linux --gwas --file combin2 --numeric --trait $i;done;





