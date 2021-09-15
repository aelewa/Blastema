### Script for blastema ribosome activity paper ###

source('/Users/salnewt/Documents/02_stockholm/12_scripts/WABI/Rscripts/functions.R')
source('/Users/salnewt/Documents/02_stockholm/12_scripts/WABI/Rscripts/pca_plotting_v2.R')
library(plotrix)
setwd('/Users/salnewt/Documents/02_stockholm/07_projects/01_collaborations/Brito_miR-10/')
df_uA <- read.table('matrix_uA',header=T)
df_mA <- read.table('matrix_mA',header=T)
df_aA <- df_uA + df_mA

# Identify data for ERCC spike-ins and non-coding RNAs
rna_aA  <-grep("Infernal",rownames(df_aA))
ercc_aA <-grep("^ERCC-[[:digit:]]",rownames(df_aA))

# Number of detected genes and ERCC spike-ins
nD.genes      <- apply(df_aA[-c(ercc_aA,rna_aA),],2,detect.genes)
nD.ercc       <- apply(df_aA[ercc_aA,],2,detect.genes)

# Filter out cells with less than 10 ERCC
passed <- (nD.ercc >=10)

# Colors
popAcolor <- rgb((252/255),(169/255),(133/255),0.5)
popBcolor <- rgb((133/255),(202/255),(93/255),0.5)
ribocolor <- rgb((1/255),(100/255),(255/255),0.5)
TEcolor   <- rgb((255/255),(1/255),(100/255),0.5)
# Keep protein coding genes with more than one mapped read
df_aA2 <- df_aA[-c(ercc_aA,rna_aA),]
df_aA2 <- df_aA2[rownames(df_aA2)!='UNKNOWN',]
df_aA2 <- df_aA2[rownames(df_aA2)!='unknown',]
keep<-which(rowSums(df_aA2)>1)

# Normalize using kept genes and passed cells
df_aAn <- sweep(df_aA2[keep,passed], 2, colSums(df_aA2[keep,passed]), FUN="/") *1000000

# Define two populations
popA <- colnames(df_aA2[keep,nD.genes >= 2000 & nD.ercc >=10])
popB <- colnames(df_aA2[keep,nD.genes < 2000  & nD.ercc >=10])

# Genes of interest
# Prod1 expressing cells (i.e. with two or more reads mapping to Cluster@@ - based on prior analysis)
prod1 <- c('P1_B18','P1_B23','P1_B5','P1_D11','P1_D13','P1_D9','P1_H10','P1_J15','P1_K6','P1_L16','P1_L9',
           'P2_C6','P2_G21','P2_H9','P2_I18','P2_I19','P2_I21','P2_I23','P2_L21','P2_N17','P2_O21')
ribogenes <- c('RL10','RL10A','RL10L','RL11','RL12','RL13','RL13A','RL14','RL15','RL17','RL18','RL18A',
               'RL19','RL21','RL22','RL22L','RL23','RL23A','RL24','RL26','RL26L','RL27','RL27A','RL28',
               'RL29','RL3','RL30','RL31','RL32','RL34','RL35','RL35A','RL36','RL36A','RL37','RL37A','RL38',
               'RL39','RL4','RL40','RL4B','RL5','RL5B','RL6','RL7','RL7A','RL7L','RL8','RL9','RLA0','RLA1',
               'RLA2','RS10','RS11','RS12','RS13','RS14','RS15','RS15A','RS16','RS17','RS18','RS19','RS2',
               'RS20','RS21','RS23','RS24','RS25','RS26','RS27','RS27A','RS27L','RS28','RS29','RS3','RS30',
               'RS31','RS3A','RS3AB','RS4','RS4X','RS4Y1','RS5','RS6','RS7','RS8','RS9')
TEs <- c('DPOL','ENR1','ENV','ERVV1','ERVV2','GAG','HARB1','LIN1','LITD1','LORF1','LORF2','POL',
         'POL1','POL4','RTL1','Transposon_cacta','Transposon_Crypton','Transposon_DDE_1','Transposon_gypsy',
         'Transposon_hAT','Transposon_helitronORF','Transposon_ISb','Transposon_ISC1316','Transposon_LINE',
         'Transposon_ltr_Roo','Transposon_mariner','Transposon_MuDR_A_B','Transposon_P_element',
         'Transposon_piggybac','Transposon_TY1_Copia')
# NDUs
NDUs <- rownames(df_aAn[substr(rownames(df_aAn),1,3)=='NDU',])
NDUtop8 <- c('NDUAB','NDUA3','NDUBB','NDUV2','NDUAD','NDUB9','NDUA6','NDUA1')

# PCA (all cells)
PC_aAn <-run.pca(df_aAn)

# PC contributions and top contributers to each PC (all cells)
PC_aAn.vars<- PC_aAn$sdev^2
PC_aAn.vars<- PC_aAn.vars/sum(PC_aAn.vars)
PC_aAn.aload <- abs(PC_aAn$rotation[,1:5])
PC_aAn.contr<- sweep(PC_aAn.aload, 2, colSums(PC_aAn.aload), "/")

# PCA (popB)
PC_aAnB <-run.pca(df_aAn[,popB])

# PC contributions and top contributers to each PC (popB)
PC_aAnB.vars<- PC_aAnB$sdev^2
PC_aAnB.vars<- PC_aAnB.vars/sum(PC_aAnB.vars)
PC_aAnB.aload <- abs(PC_aAnB$rotation[,1:5])
PC_aAnB.contr<- sweep(PC_aAnB.aload, 2, colSums(PC_aAnB.aload), "/")

# Negative correlation between ribogenes and TEs
cor.test(colMeans(df_aAn[TEs,]),colMeans(df_aAn[ribogenes,]))
cor.test(colMeans(df_aAn[TEs,popA]),colMeans(df_aAn[ribogenes,popA]))
cor.test(colMeans(df_aAn[TEs,popB]),colMeans(df_aAn[ribogenes,popB]))

cor.test(colMeans(df_aAn[TEs,popB][substr(popB,1,2)=='P1']),colMeans(df_aAn[ribogenes,popB][substr(popB,1,2)=='P1']))
substr(popB,1,2)

############################
##        Figures         ##
############################

# Fig 1b.Number of detected genes vs ERCC
pdf('20210324_Fig1b.pdf')
par(mar=c(4,4,3,20))
plot(nD.ercc[passed],nD.genes[passed],main='',
     pch=20,
     col=ifelse(names(nD.genes[passed])%in%popA,rgb((252/255),(169/255),(133/255),0.5),rgb((133/255),(202/255),(93/255),0.5)),
     ylim=c(0,8500),ylab = 'detected genes',xlab='detected ERCC', las=1,cex.axis=1)
legend('topleft',legend=c('kiloscript','hectoscript'),
       col=(c(rgb((252/255),(169/255),(133/255),0.5),rgb((133/255),(202/255),(93/255),0.5))),pch=20,bty='n')
dev.off()


# Fig 1c. PC1 (all cells) colored by popA vs popB
pdf('20210324_Fig1c.pdf')
par(mar = c(10,10,14,10))
plot(PC_aAn$x[,'PC1'],nD.genes[passed],pch=20,
     col=ifelse(names(nD.genes[passed])%in%popA,popAcolor,popBcolor),
     main='',las=2,cex.axis=1,ylab='detected genes',xlab='PC1')
abline(v=0,lty=3)
legend('topleft',legend='kiloscript',bty='n',cex=1.4)
legend('topright',legend='hectoscript',bty='n',cex=1.4)
dev.off()

# Fig 1d. MKNK2 is a marker for popB
pdf('/Users/salnewt/Documents/02_stockholm/07_projects/01_collaborations/Brito_miR-10/20210501_Manuscript/Fig1f.pdf')
par(mar = c(10,10,14,10),mgp=c(2,0.6,0))
plot(PC_aAn$x[,'PC1'],as.numeric(df_aAn['MKNK2',]),
     col=ifelse(names(nD.genes[passed])%in%popA,popAcolor,popBcolor),yaxt='n',
     pch=20,main='',cex.axis=0.7,ylab='MKNK2 expression CPM/1000',xlab='PC1')
axis(2,at=c(50000,100000,150000,200000),labels=c('50','100','150','200'),las=2,cex.axis=0.7)
abline(v=0,lty=3)
text(-200,200000,'Kiloscript')
text(105,200000,'Hectoscript')
dev.off()

# Fig 1e. Top contributors to PC1 (all cells)
pdf('20210324_Fig1e.pdf')
par(mar = c(10,10,14,10))
nPlot <- 20 #number of genes top plot per pc.
for (i in 1:1) {
  top<-order(PC_aAn$rotation[,i],decreasing=T)[1:nPlot]
  bottom<-order(PC_aAn$rotation[,i],decreasing=F)[1:nPlot]
  barplot(PC_aAn.contr[top,i],main='',ylab='',
          las=2,horiz=T,cex.axis=1,cex.names=0.7)
  #barplot(PC_aAn.contr[bottom,i],main='',ylab='',
  #        las=2,horiz=T,cex.axis=0.5,cex.names=0.5)
}
dev.off()

# Fig 1f. MKNK2 top expressor
par(mar = c(10,10,14,10),mgp=c(3,0.6,0))
barplot(rowMeans(df_aAn[c('SMUF2','MKNK2','TPD52','WDR26','QRIC1','CHSS1','MARCS','KS6B1'),popB]),las=2,
        ylab='Hectoscript expression (Mean CPM)')
dev.off()


pdf('/Users/salnewt/Documents/02_stockholm/07_projects/01_collaborations/Brito_miR-10/20210501_Manuscript/Fig1e.pdf')
par(mar = c(10,10,14,10),mgp=c(3.5,0.6,0))
boxplot(as.numeric(df_aAn['SMUF2',]),as.numeric(df_aAn['MKNK2',]),as.numeric(df_aAn['TPD52',]),
        as.numeric(df_aAn['WDR26',]),as.numeric(df_aAn['QRIC1',]),as.numeric(df_aAn['CHSS1',]),
        as.numeric(df_aAn['MARCS',]),as.numeric(df_aAn['KS6B1',]),outline = F,las=2,
        ylab='Hectoscript expression CPM',xaxt='n',ylim=c(0,200000),frame.plot=T,cex.axis=0.7)
axis(1,at=1:8,labels = c('SMURF2','MKNK2','TPD52','WDR26','QRIC1','CHSS1','MARCS','KS6B1'),las=2,cex.axis=0.7)
        #'MKNK2','TPD52','WDR26','QRIC1','CHSS1','MARCS','KS6B1'),popB]),las=2,
        #ylab='Hectoscript expression (Mean CPM)')
dev.off()


# Fig 1g. TE branch with ribogene marker
pdf('20210324_Fig1g.pdf')
par(mar = c(10,10,14,10))
pca.plot(PC_aAnB,pch=16,main='',
         col=color.scale(log2(colMeans(df_aAn[TEs,popB])),c(0,1,1),c(1,1,0),0),xrange(0,14),las=2,cex.axis=1)
dev.off()

# Fig 1h. Ribosome branch with ribogene marker
pdf('20210324_Fig1h.pdf')
par(mar = c(10,10,14,10))
pca.plot(PC_aAnB,pch=16,main='',
         col=color.scale(log(colMeans(df_aAn[ribogenes,popB])),c(0,1,1),c(1,1,0),0),xrange(0,14),las=2,cex.axis=1)
dev.off()

# legend
pdf('20210324_Fig1g-h_legend.pdf')
par(mar = c(10,10,14,10))
plot(0:14,rep(1,15),pch=15,ylim=c(1,100),cex=10,col=color.scale(0:14,c(0,1,1),c(1,1,0),0),yaxt='n',bty='n',ylab='',xlab='')
dev.off()

# Fig 1i.rise of ribo, fall of TE
pdf('/Users/salnewt/Documents/02_stockholm/07_projects/01_collaborations/Brito_miR-10/20210501_Manuscript/Fig1h.pdf')
par(mar = c(10,10,14,10),mgp=c(1.5,0.6,0))
plot(PC_aAn$x[,'PC1'],log(colMeans(df_aAn[ribogenes,])),
     ylim=c(4.5,10.5),
     col=ribocolor,
     pch=20,main='',las=1,cex.axis=0.7,ylab='ln(meanCPM)',xlab='PC1')

points(PC_aAn$x[,'PC1'],log(colMeans(df_aAn[TEs,])),
       col=TEcolor,pch=16)
splineTE <- smooth.spline(PC_aAn$x[,'PC1'],log(colMeans(df_aAn[TEs,])), df=4)
predictTE <- predict(splineTE, PC_aAn$x[,'PC1'])
splineRibogene <- smooth.spline(PC_aAn$x[,'PC1'],log(colMeans(df_aAn[ribogenes,])), df=4)
predictRibogene <- predict(splineRibogene, PC_aAn$x[,'PC1'])
lines(splineTE,col=TEcolor,lwd=6)
lines(splineRibogene,col=ribocolor,lwd=6)
legend('topleft',pch=c(20,16),col=c(ribocolor,TEcolor),c('Ribogene','TE'),bty='n')
abline(v=0,lty=3)
dev.off()

#par(mar = c(4,4,4,1),mgp=c(2.1,0.5,0))
#pca.plot(PC_aAnB,pch=20,main='',
#         col=color.scale(log(colMeans(df_aAn[ribogenes,popB])),c(0,1,1),c(1,1,0),0),xrange(0,14),las=2,cex.axis=0.6)


#par(mar = c(4,4,4,1),mgp=c(2.1,0.5,0))
#pca.plot(PC_aAnB,pch=20,main='',
#         col=color.scale(log(colMeans(df_aAn[TEs,popB])),c(0,1,1),c(1,1,0),0),xrange(0,14),las=2,cex.axis=0.6)


## Supplemetnary Figure 1
dev.off()
pdf('Manuscript 202103/FigSupp1a.pdf')
par(mfrow=c(1,1))
par(mar=c(10,5,10,5))
colbatch = c('chartreuse3','darkgoldenrod')
pca.plot(PC_aAn,pch=16,main='PCA and Batch Effect',col=colbatch[as.numeric(substr((colnames(df_aAn)),2,2))],
         ylim=c(-220,220),xlim=c(-420,300),las=1)
dev.off()
pdf('Manuscript 202103/FigSupp1b.pdf')
par(mar=c(4,4,4,4))
nPlot <- 20 #number of genes top plot per pc.
par(mfrow=c(2,2),mar=c(2,2,1,1),oma=c(1,5,1,1))
for (i in 1:1) {
  top<-order(PC_aAn$rotation[,i],decreasing=T)[1:nPlot]
  bottom<-order(PC_aAn$rotation[,i],decreasing=F)[1:nPlot]
  #barplot(PC_aAn.contr[top,i],main='',ylab='',
  #        las=2,horiz=T,cex.axis=0.5,cex.names=0.5)
  barplot(PC_aAn.contr[bottom,i],main='Top Contributors PC1-ve',ylab='',
          las=2,horiz=T,cex.axis=1,cex.names=1)
}
dev.off()

# PCA (all cells) colored by number of detected genes
#cc<-colorRampPalette(c("yellow","red","black"))
#col.nD<-convert.to.color(nD.genes[passed],cc)
#pca.plot(PC_aAn,pch=16,main='',col=col.nD$cols)

# Fig 1c. PC1 (all cells) colored by popA vs popB
plot(PC_aAn$x[,'PC1'],nD.genes[passed],pch=20,
     col=ifelse(names(nD.genes[passed])%in%popA,popAcolor,popBcolor),
     main='',las=2,cex.axis=0.6,ylab='',xlab='')

# Fig 1d. MKNK2 is a marker for popB
#plot(PC_aAn$x[,'PC1'],as.numeric(df_aAn['MKNK2',]),
#     col=ifelse(names(nD.genes[passed])%in%popA,popAcolor,popBcolor),
#     pch=20,main='',las=2,cex.axis=0.6,ylab='',xlab='')




#  Contributions of PCs (all cells)
#barplot(PC_aAn.vars[1:9],names.arg=1:9,main="",cex.axis = 0.6,cex.names = 0.6,las=1,
#        ylab='',xlab='')

#  Contributions of PCs (popB)
#barplot(PC_aAnB.vars[1:9],names.arg=1:9,main="",cex.axis = 0.6,cex.names = 0.6, las=1,
#        ylab='',xlab='')

# Fig 1e. Top contributors to PC1 (all cells)
nPlot <- 10 #number of genes top plot per pc.
#par(mfrow=c(2,2),mar=c(2,2,1,1),oma=c(1,5,1,1))
for (i in 1:1) {
  top<-order(PC_aAn$rotation[,i],decreasing=T)[1:nPlot]
  bottom<-order(PC_aAn$rotation[,i],decreasing=F)[1:nPlot]
  barplot(PC_aAn.contr[top,i],main='',ylab='',
          las=2,horiz=T,cex.axis=0.5,cex.names=0.5)
  #barplot(PC_aAn.contr[bottom,i],main='',ylab='',
  #        las=2,horiz=T,cex.axis=0.5,cex.names=0.5)
}

# Fig 1f. Ribosome branch with ribogene marker
pca.plot(PC_aAnB,pch=20,main='',
         col=color.scale(log(colMeans(df_aAn[ribogenes,popB])),c(0,1,1),c(1,1,0),0),xrange(0,14),las=2,cex.axis=0.6)



# Fig 1g.TE branch with TE marker
pca.plot(PC_aAnB,pch=16,main='',
         col=color.scale(log2(colMeans(df_aAn[TEs,popB])),c(0,1,1),c(1,1,0),0),xrange(0,14),las=2,cex.axis=0.6)






# Fig 2d. qRT-PCR of miR-10 in Nv and Pw
pn <- read.table('fromAnoop/temp.txt',header=TRUE)
pdf('plotpPCRv2.pdf')
par(mfrow=c(3,3))
boxplot(pn[,c(5:7,1:2,4)],las=1,cex.axis=0.9,names=c(0,3,18,0,3,18),ylim=c(0.3,1.2))
abline(v=3.5,lty=2)
dev.off()
t.test(pn[,5],pn[,6]) # Nv 0vs3 t = 9.107, df = 3, p-value = 0.002798
t.test(pn[,6],pn[,7]) # Nv 3vs18 t = -4.7566, df = 6.9933, p-value = 0.002073

t.test(pn[,1],pn[,2]) # Pw 0vs3 t = 32.002, df = 5.7606, p-value = 1.039e-07
t.test(pn[,2],pn[,4]) # Pw 3vs18 t = -36.516, df = 6.499, p-value = 9.111e-09


# Fig 2g miR-10b mimic and ribogene qRT-PCR
df <- read.csv('temp5RP.csv')
unique(df$rp)
#"RL4"   "RL15"  "RL27A" "RS29"  "RL30L"
t.test(df[df$condition=='control','foldchange'],df[df$condition=='mimic','foldchange'])
x <- c(rep(1,length(df[df$condition=='control','foldchange'])),rep(2,length(df[df$condition=='mimic','foldchange'])))
y <- c(df[df$condition=='control','foldchange'],df[df$condition=='mimic','foldchange'])
par(mfrow=c(2,4))
cols <- c(rep(rainbow(1, s = 0.5),25),rep(rainbow(2, s = 0.5)[2],25))
plot(x,y,pch=16,ylim=c(0,2),xlim=c(0,3),las=1,col=cols)

############################
## Supplementary Figures  ##

# PC1 axis (all cells) MKNK2 and KS6B1 decline
plot(PC_aAn$x[,'PC1'],log(as.numeric(df_aAn['MKNK2',])),pch=16,main='KS6B1 vs MKNK2',las=1,cex.axis=0.7,
     ylim=c(2,13),ylab='CPM/cell',xlab='PC1',
     col=rgb((255/255),(1/255),(100/255),0.5))
points(PC_aAn$x[,'PC1'],log(as.numeric(df_aAn['KS6B1',])),pch=16,col=rgb((1/255),(100/255),(255/255),0.5))
points(PC_aAn$x[,'PC1'],log(colMeans(df_aAn[ribogenes,])),pch=16,col=ribocolor)

# PC2 axis (popB) with ribogene and TE log mean CPM
plot(PC_aAnB$x[,'PC2'],log(colMeans(df_aAn[TEs,popB])),pch=16,main='ribogenes vs TEs',las=1,cex.axis=0.7,
     ylim=c(4,11),ylab='mean CPM/cell',xlab='PC2',
     col=TEcolor)
points(PC_aAnB$x[,'PC2'],log(colMeans(df_aAn[ribogenes,popB])),pch=16,col=ribocolor)


# Basic QC
par(mfrow=c(2,2),mar=c(3,3,1,1),cex=0.6)
hist(nD.genes,n=100,main="Detected Genes")
plot(nD.ercc,nD.genes,main="ERCC vs Genes")
plot(nD.ercc,log10(nD.genes),main="ERCC vs Genes")

# PC1vsPC2 with batch identification
pca.plot(PC_aAn,pch=16,main='PC_aAn',col=as.numeric(substr((colnames(df_aAn)),2,2)))
par(mfrow=c(2,2))
pca.plot(PC_aAn,pch=16,main='RAC2',
         col=color.scale(log2(as.numeric(df_aAn['RAC2',]+1)),c(0,1,1),c(1,1,0),0),xrange(0,14))
pca.plot(PC_aAn,pch=16,main='GMFG',
         col=color.scale(log2(as.numeric(df_aAn['GMFG',]+1)),c(0,1,1),c(1,1,0),0),xrange(0,14))
pca.plot(PC_aAn,pch=16,main='ARPC3',
         col=color.scale(log2(as.numeric(df_aAn['ARPC3',]+1)),c(0,1,1),c(1,1,0),0),xrange(0,14))
pca.plot(PC_aAn,pch=16,main='CD53',
         col=color.scale(log2(as.numeric(df_aAn['CD53',]+1)),c(0,1,1),c(1,1,0),0),xrange(0,14))




# RAC2
plot(PC_aAn$x[,'PC1'],as.numeric(df_aAn['RAC2',]),
     col=ifelse(names(nD.genes[passed])%in%popA,popAcolor,popBcolor),
     pch=16,main='RAC2',las=1,cex.axis=0.7,ylab='',xlab='')

plot(nD.ercc[passed],nD.genes[passed],main="Blastema scRNAseq",
     pch=ifelse(names(nD.genes[passed])%in%prod1,15,20),
     col=color.scale(log2(as.numeric(df_aAn['RAC2',]+1)),c(0,1,1),c(1,1,0),0),
     xlab = "Detected ERCC",ylab="Detected genes", las=1,cex.axis=0.8)

# PC1 of MKNK2+ cells
par(mfrow=c(3,3))
top<-order(PC_aAnB$rotation[,1],decreasing=T)[1:9]
bottom<-order(PC_aAnB$rotation[,1],decreasing=F)[1:9]
goi <- rownames(PC_aAnB.contr[top,])
for (i in goi){
  plot(PC_aAnB$x[,'PC1'],as.numeric(df_aAn[i,popB]),pch=16,main=i,las=1,cex.axis=0.7,ylab='',xlab='')
}
goi <- rownames(PC_aAnB.contr[bottom,])
for (i in goi){
  plot(PC_aAnB$x[,'PC1'],as.numeric(df_aAn[i,popB]),pch=16,main=i,las=1,cex.axis=0.7,ylab='',xlab='')
}

# PC2 of MKNK2+ cells
par(mfrow=c(3,3))
top<-order(PC_aAnB$rotation[,2],decreasing=T)[1:9]
bottom<-order(PC_aAnB$rotation[,2],decreasing=F)[1:9]
goi <- rownames(PC_aAnB.contr[top,])
for (i in goi){
  plot(PC_aAnB$x[,'PC2'],as.numeric(df_aAn[i,popB]),pch=16,main=i,las=1,cex.axis=0.7,ylab='',xlab='')
}
goi <- rownames(PC_aAnB.contr[bottom,])
for (i in goi){
  plot(PC_aAnB$x[,'PC2'],as.numeric(df_aAn[i,popB]),pch=16,main=i,las=1,cex.axis=0.7,ylab='',xlab='')
}

# PCA plot PC1
par(mfrow=c(3,3))
top<-order(PC_aAnB$rotation[,1],decreasing=T)[1:3]
bottom<-order(PC_aAnB$rotation[,1],decreasing=F)[1:6]
goi <- rownames(PC_aAnB.contr[top,])
for (i in goi){
  pca.plot(PC_aAnB,pch=16,main=paste('PC_aAnB',i),
           col=color.scale(log2(as.numeric(df_aAn[i,popB]+1)),c(0,1,1),c(1,1,0),0),xrange(0,14))
}
goi <- rownames(PC_aAnB.contr[bottom,])
for (i in goi){
  pca.plot(PC_aAnB,pch=16,main=paste('PC_aAnB',i),
           col=color.scale(log2(as.numeric(df_aAn[i,popB]+1)),c(0,1,1),c(1,1,0),0),xrange(0,14))
}

# PCA plot PC2
par(mfrow=c(3,3))
top<-order(PC_aAnB$rotation[,2],decreasing=T)[1:3]
bottom<-order(PC_aAnB$rotation[,2],decreasing=F)[1:6]
goi <- rownames(PC_aAnB.contr[top,])
for (i in goi){
  pca.plot(PC_aAnB,pch=16,main=paste('PC_aAnB',i),
           col=color.scale(log2(as.numeric(df_aAn[i,popB]+1)),c(0,1,1),c(1,1,0),0),xrange(0,14))
}
goi <- rownames(PC_aAnB.contr[bottom,])
for (i in goi){
  pca.plot(PC_aAnB,pch=16,main=paste('PC_aAnB',i),
           col=color.scale(log2(as.numeric(df_aAn[i,popB]+1)),c(0,1,1),c(1,1,0),0),xrange(0,14))
}



par(mfrow=c(1,1))
plot(PC_aAn$x[,'PC1'],log(colMeans(df_aAn[NDUs,])),
     ylim=c(1,9),
     col=rgb((1/255),(100/255),(255/255),0.2),
     pch=16,main='NDUs',las=1,cex.axis=0.7,ylab='',xlab='')

points(PC_aAn$x[,'PC1'],log(colMeans(df_aAn[ribogenes,])),
       col=rgb((255/255),(1/255),(100/255),0.2),
       pch=16)
splineNDUs <- smooth.spline(PC_aAn$x[,'PC1'],log(colMeans(df_aAn[NDUs,]+1)), df=4)
splineRibogene <- smooth.spline(PC_aAn$x[,'PC1'],log(colMeans(df_aAn[ribogenes,]+1)), df=4)
lines(splineNDUs,col=rgb((1/255),(100/255),(255/255),0.5),lwd=6)
lines(splineRibogene,col=rgb((255/255),(1/255),(100/255),0.5),lwd=6)

cor.test(colMeans(df_aAn[ribogenes,]),colMeans(df_aAn[NDUs,]))
cor.test(splineNDUs$y,splineRibogene$y)

# PCA (popB) colored by number of detected genes
#col.nDB<-convert.to.color(nD.genes[names(nD.genes)%in%popB],cc)
#pca.plot(PC_aAnB,pch=16,main='PC_aAnB',col=col.nDB$cols)

# Ribosome branch with RS23 marker
#pca.plot(PC_aAnB,pch=16,main='RS23',
#         col=color.scale(log2(as.numeric(df_aAn['RS23',popB]+1)),c(0,1,1),c(1,1,0),0),xrange(0,14))
# TE branch with PiggyBac marker
#pca.plot(PC_aAnB,pch=16,main='PiggyBac',
#         col=color.scale(log2(as.numeric(df_aAn['Transposon_piggybac',popB]+1)),c(0,1,1),c(1,1,0),0),xrange(0,14))


# PC1 axis (popB) MKNK2 and KS6B1 decline
plot(PC_aAnB$x[,'PC1'],log(as.numeric(df_aAn['MKNK2',popB])),pch=16,main='KS6B1 vs MKNK2',las=1,cex.axis=0.7,
     ylim=c(1.8,12.3),ylab='CPM/cell',xlab='PC1',
     col='grey')
points(PC_aAnB$x[,'PC1'],log(as.numeric(df_aAn['KS6B1',popB])),pch=16,col='turquoise')
points(PC_aAnB$x[,'PC1'],log(colMeans(df_aAn[ribogenes,popB])),pch=16,col=ribocolor)
points(PC_aAnB$x[,'PC1'],log(colMeans(df_aAn[TEs,popB])),pch=16,col=TEcolor)
splineMKNK2 <- smooth.spline(PC_aAnB$x[,'PC1'],log(as.numeric(df_aAn['MKNK2',popB]+1)), df=4)
splineKS6B1 <- smooth.spline(PC_aAnB$x[,'PC1'],log(as.numeric(df_aAn['KS6B1',popB]+1)), df=4)
splineRibo  <- smooth.spline(PC_aAnB$x[,'PC1'],log(colMeans(df_aAn[ribogenes,popB])), df=4)
splineTEs   <- smooth.spline(PC_aAnB$x[,'PC1'],log(colMeans(df_aAn[TEs,popB])), df=4)

lines(splineMKNK2,col='grey',lwd=6)
lines(splineKS6B1,col='turquoise',lwd=6)
lines(splineRibo,col=ribocolor,lwd=6)
lines(splineTEs,col=TEcolor,lwd=6)

pdf('OPPhistogram.pdf')
par(mfrow=c(2,2))
hist(as.numeric(df$`Control-mimic`),col=rgb((252/255),(169/255),(133/255),0.5),10,ylim=c(0,20),
     las=1,xlab='OPP intensity',main='miR-10 mimic reduces OPP intensity')
hist(as.numeric(df$`10b-mimic`),add=T,col=rgb((133/255),(202/255),(93/255),0.5))
legend('right',legend=c('control','mimic'),
       col=c(rgb((252/255),(169/255),(133/255),0.5),rgb((133/255),(202/255),(93/255),0.5)),
       pch=15,bty='n',cex=1.5)
dev.off()

df_brdu_mimic <- read.table('brducount.txt')
df_brdu_mimic <- df_brdu_mimic /38.1
pdf('20210324_Manuscript/20210324_Fig4i.pdf')
par(mfrow=c(2,1))
boxplot(rev(df_brdu_mimic),las=1,xlab='Relative BrdU intensity',cex.axis=1.5,horizontal = T,ylim=c(0.3,1.5))
dev.off()
t.test(df_brdu_mimic$Control,df_brdu_mimic$Mimic)
#Welch Two Sample t-test

#data:  df_brdu_mimic$Control and df_brdu_mimic$Mimic
#t = 3.2519, df = 11.46, p-value = 0.00733
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  0.1189014 0.6095107
#sample estimates:
#  mean of x mean of y 
#0.9769357 0.6127297 

df_opp_mimic <- read.csv('/Users/salnewt/Documents/02_stockholm/07_projects/01_collaborations/Brito_miR-10/mir10b_OPP data.csv')
colnames(df_opp_mimic) <- c('Control', 'Mimic')
pdf('/Users/salnewt/Documents/02_stockholm/07_projects/01_collaborations/Brito_miR-10/20210324_Manuscript/20210324_Fig4c.pdf')
par(mfrow=c(1,2))
boxplot(df_opp_mimic, outline = F,las=1,cex.axis=1.5,ylab='fluorescence intensity (arbitrary units)')
dev.off()
t.test(df_opp_mimic$Control,df_opp_mimic$Mimic)
#Welch Two Sample t-test

#data:  df_opp_mimic$Control and df_opp_mimic$Mimic
#t = 3.1214, df = 93.444, p-value = 0.002395
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  4.062324 18.266883
#sample estimates:
#  mean of x mean of y 
#27.56787  16.40327 

df_opp_cerco <- read.table('/Users/salnewt/Documents/02_stockholm/07_projects/01_collaborations/Brito_miR-10/20210324_Manuscript/cercoOPP.tab',header=TRUE)
pdf('/Users/salnewt/Documents/02_stockholm/07_projects/01_collaborations/Brito_miR-10/20210324_Manuscript/20210324_Fig4cerco.pdf')
par(mfrow=c(1,2))
boxplot(df_opp_cerco, outline = F,las=1,cex.axis=1.5,ylab='fluorescence intensity (arbitrary units)')
dev.off()
t.test(df_opp_cerco$Control,df_opp_cerco$Cerco)
#Welch Two Sample t-test

#data:  df_opp_cerco$Control and df_opp_cerco$Cerco
#t = -6.1426, df = 1094.4, p-value = 1.134e-09
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -11.652208  -6.010283
#sample estimates:
#  mean of x mean of y 
#44.47299  53.30424 


df_brdu_cerco <- read.table('/Users/salnewt/Documents/02_stockholm/07_projects/01_collaborations/Brito_miR-10/20210324_Manuscript/cercoBrdU.tab',header=TRUE)
pdf('/Users/salnewt/Documents/02_stockholm/07_projects/01_collaborations/Brito_miR-10/20210324_Manuscript/20210324_Fig2Brdu.pdf')
par(mfrow=c(1,2))
boxplot(df_brdu_cerco, outline = F,las=1,cex.axis=1.5,ylab='fluorescence intensity (arbitrary units)')
dev.off()
t.test(df_brdu_cerco$Control,df_brdu_cerco$Cerco)

df <- matrix(nrow=6,ncol=2,c(0.7,0.3,0.7,0.74,23.7,0.4,155.3,147.1,347.2,123.5,247.5,59.3))
pdf('/Users/salnewt/Documents/02_stockholm/07_projects/01_collaborations/Brito_miR-10/20210324_Manuscript/Supp3.pdf')
par(mar=c(4,4,4,15))
boxplot(df,las=1,xaxt='n',ylab='Relative expression (fold change)',ylim=c(0,400))
axis(1,at=c(1,2),labels=c('control','mimic'),cex.axis=1.2)
segments(1,380,2,380)
segments(1,380,1,375)
segments(2,380,2,375)
text(1.5,385,'**',cex=1.5)
dev.off()


# Does number of detected genes define PC1 in hectoscripts?
pdf('20210501_Manuscript/PC1ofPC_aAnB.pdf')
par(mar = c(10,10,14,10))
plot(PC_aAnB$x[,'PC1'],nD.genes[popB],pch=20,
     col=ifelse(names(nD.genes[popB])%in%popA,popAcolor,popBcolor),
     main='PCA of Hectoscripts',las=2,cex.axis=1,ylab='detected genes',xlab='PC1')
dev.off()

# Top contributers to PCs (popB)
pdf('20210501_Manuscript/TopContributorsofPC2.pdf')
nPlot <- 40 #number of genes top plot per pc.
par(mfrow=c(2,2),mar=c(2,2,1,1),oma=c(1,5,1,1))
for (i in 2:2) {
  top<-order(PC_aAnB$rotation[,i],decreasing=T)[1:nPlot]
  bottom<-order(PC_aAnB$rotation[,i],decreasing=F)[1:nPlot]
  barplot(PC_aAnB.contr[top,i],main='',ylab='',
          las=2,horiz=T,cex.axis=0.5,cex.names=0.5)
  barplot(PC_aAnB.contr[bottom,i],main='',ylab='',
          las=2,horiz=T,cex.axis=0.5,cex.names=0.5)
}
dev.off()
