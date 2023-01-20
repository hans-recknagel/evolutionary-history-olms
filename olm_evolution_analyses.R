<<<<<<< HEAD
# load libraries
library(adegenet)
library(Amelia)
library(ape)
library(bppr)
library(caret)
library(castor)
library(chopper)
library(corHMM)
library(corrplot)
library(dartR)
library(diversitree)
library(extRemes)
library(factoextra)
library(FactoMineR)
library(ggfortify)
library(gridExtra)
library(lmtest)
library(MASS)
library(MCMCtreeR)
library(mice)
library(missForest)
library(missMDA)
library(naniar)
library(OptM)
library(pals)
library(pegas)
library(phangorn)
library(pheatmap)
library(phylotools)
library(plot.matrix)
library(poppr)
library(radiator)
library(RColorBrewer)
library(Rcpp)
library(reshape2)
library(seqinr)
library(smatr)
library(SNPRelate)
library(tidyverse)
library(VIM)
library(visdat)
#install_local("C:/Users/hans_/Downloads/radiator-master.zip")

########### GENOTYPING ########### 
### error rate estimation
setwd("01_replicates")
proteus.error.wl <- read.PLINK("replicates_whitelist_1-3snps_r66.raw", parallel=FALSE)
source("C:/Users/hans_/Desktop/Projects/Proteus_phylogeny/04_analysis/08_replicates/PairsDiff.R")
source("C:/Users/hans_/Desktop/Projects/Proteus_phylogeny/04_analysis/08_replicates/SNPs_error.R") # read in the R functions, which also calls the needed packages
SNP_error(proteus.error.wl)
# transpose data to analyse in excel and find those loci that do not match
rep_256 <- read.table("256_replicates.txt", sep="\t", header=TRUE)
rep_256_t <- t(rep_256) 
rep_592 <- read.table("592_replicates.txt", sep="\t", header=TRUE)
rep_592_t <- t(rep_592) 
#write.table(rep_256_t,"rep_256_t.txt")
#write.table(rep_592_t,"rep_592_t.txt")


########### POPULATION GENETICS ########### 
#### ADMIXTURE ####
setwd("02_admixture")
n <- 10
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))

colours.ten<-sample(col_vector, n)
par(mfrow=c(1,1))
data<-read.table("proteus_p1_r66_minmaf05_locus_chr.8.Q")
barplot(t(as.matrix(data)), col=colours.ten,
        xlab="Individual #", ylab="Ancestry", border=NA,las=2)
barplot(t(as.matrix(data)), col=rainbow(5),
        xlab="Individual #", ylab="Ancestry", border=NA,las=2)

## create admixture plots 
# for Q = 2
tbl_2<-read.table("proteus_p1_r66_minmaf05_locus_chr.2.Q")
fam <- read.table("proteus_p1_r66_minmaf05_locus_chr.fam")
rownames(tbl_2) <- fam$V2
colnames(tbl_2) <- c("V01", "V02")
# for Q = 3
tbl_3<-read.table("proteus_p1_r66_minmaf05_locus_chr.3.Q")
fam <- read.table("proteus_p1_r66_minmaf05_locus_chr.fam")
rownames(tbl_3) <- fam$V2
colnames(tbl_3) <- c("V01", "V02", "V03")
# for Q = 4
tbl_4<-read.table("proteus_p1_r66_minmaf05_locus_chr.4.Q")
fam <- read.table("proteus_p1_r66_minmaf05_locus_chr.fam")
rownames(tbl_4) <- fam$V2
colnames(tbl_4) <- c("V01", "V02", "V03", "V04")
# for Q = 5
tbl_5<-read.table("proteus_p1_r66_minmaf05_locus_chr.5.Q")
fam <- read.table("proteus_p1_r66_minmaf05_locus_chr.fam")
rownames(tbl_5) <- fam$V2
colnames(tbl_5) <- c("V01", "V02", "V03", "V04", "V05")
# for Q = 6
tbl_6<-read.table("proteus_p1_r66_minmaf05_locus_chr.6.Q")
fam <- read.table("proteus_p1_r66_minmaf05_locus_chr.fam")
rownames(tbl_6) <- fam$V2
colnames(tbl_6) <- c("V01", "V02", "V03", "V04", "V05", "V06")
# for Q = 7
tbl_7<-read.table("proteus_p1_r66_minmaf05_locus_chr.7.Q")
fam <- read.table("proteus_p1_r66_minmaf05_locus_chr.fam")
rownames(tbl_7) <- fam$V2
colnames(tbl_7) <- c("V01", "V02", "V03", "V04", "V05", "V06", "V07")
# for Q = 8
tbl_8<-read.table("proteus_p1_r66_minmaf05_locus_chr.8.Q")
fam <- read.table("proteus_p1_r66_minmaf05_locus_chr.fam")
rownames(tbl_8) <- fam$V2
colnames(tbl_8) <- c("V01", "V02", "V03", "V04", "V05", "V06", "V07", "V08")
# for Q = 9
tbl_9<-read.table("proteus_p1_r66_minmaf05_locus_chr.9.Q")
fam <- read.table("proteus_p1_r66_minmaf05_locus_chr.fam")
rownames(tbl_9) <- fam$V2
colnames(tbl_9) <- c("V01", "V02", "V03", "V04", "V05", "V06", "V07", "V08", "V09")
# for Q = 10
tbl_10<-read.table("proteus_p1_r66_minmaf05_locus_chr.10.Q")
fam <- read.table("proteus_p1_r66_minmaf05_locus_chr.fam")
rownames(tbl_10) <- fam$V2
colnames(tbl_10) <- c("V01", "V02", "V03", "V04", "V05", "V06", "V07", "V08", "V09","V10")
# using tidyverse
plot_data_2 <- tbl_2 %>% 
  mutate(id = rownames(tbl_2)) %>% 
  gather('pop', 'prob', V01:V02) %>% 
  group_by(id) %>% 
  mutate(likely_assignment = pop[which.max(prob)],
         assingment_prob = max(prob)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))
plot_data_2$id <- factor(plot_data_2$id,levels = c("PA331", "PC622", "PC621", "PC623", "PA329", "PA332", "PA334", "PC580", "PC581", "PC616", "PC617", "PC618", "PA040", "PA047", "PA049", "PA341", "PC587", "PC592", "PC614", "PC631-r", "PA354", "PC630", "PA255", "PA256", "PA343", "PA346", "PC582", "PA363", "PC578", "PC579", "PC586", "PC593"))

plot_data_3 <- tbl_3 %>% 
  mutate(id = rownames(tbl_3)) %>% 
  gather('pop', 'prob', V01:V03) %>% 
  group_by(id) %>% 
  mutate(likely_assignment = pop[which.max(prob)],
         assingment_prob = max(prob)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))
plot_data_3$id <- factor(plot_data_3$id,levels = c("PA331", "PC622", "PC621", "PC623", "PA329", "PA332", "PA334", "PC580", "PC581", "PC616", "PC617", "PC618", "PA040", "PA047", "PA049", "PA341", "PC587", "PC592", "PC614", "PC631-r", "PA354", "PC630", "PA255", "PA256", "PA343", "PA346", "PC582", "PA363", "PC578", "PC579", "PC586", "PC593"))

plot_data_4 <- tbl_4 %>% 
  mutate(id = rownames(tbl_4)) %>% 
  gather('pop', 'prob', V01:V04) %>% 
  group_by(id) %>% 
  mutate(likely_assignment = pop[which.max(prob)],
         assingment_prob = max(prob)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))
plot_data_4$id <- factor(plot_data_4$id,levels = c("PA331", "PC622", "PC621", "PC623", "PA329", "PA332", "PA334", "PC580", "PC581", "PC616", "PC617", "PC618", "PA040", "PA047", "PA049", "PA341", "PC587", "PC592", "PC614", "PC631-r", "PA354", "PC630", "PA255", "PA256", "PA343", "PA346", "PC582", "PA363", "PC578", "PC579", "PC586", "PC593"))

plot_data_5 <- tbl_5 %>% 
  mutate(id = rownames(tbl_5)) %>% 
  gather('pop', 'prob', V01:V05) %>% 
  group_by(id) %>% 
  mutate(likely_assignment = pop[which.max(prob)],
         assingment_prob = max(prob)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))
plot_data_5$id <- factor(plot_data_5$id,levels = c("PA331", "PC622", "PC621", "PC623", "PA329", "PA332", "PA334", "PC580", "PC581", "PC616", "PC617", "PC618", "PA040", "PA047", "PA049", "PA341", "PC587", "PC592", "PC614", "PC631-r", "PA354", "PC630", "PA255", "PA256", "PA343", "PA346", "PC582", "PA363", "PC578", "PC579", "PC586", "PC593"))

plot_data_6 <- tbl_6 %>% 
  mutate(id = rownames(tbl_6)) %>% 
  gather('pop', 'prob', V01:V06) %>% 
  group_by(id) %>% 
  mutate(likely_assignment = pop[which.max(prob)],
         assingment_prob = max(prob)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))
plot_data_6$id <- factor(plot_data_6$id,levels = c("PA331", "PC622", "PC621", "PC623", "PA329", "PA332", "PA334", "PC580", "PC581", "PC616", "PC617", "PC618", "PA040", "PA047", "PA049", "PA341", "PC587", "PC592", "PC614", "PC631-r", "PA354", "PC630", "PA255", "PA256", "PA343", "PA346", "PC582", "PA363", "PC578", "PC579", "PC586", "PC593"))

plot_data_7 <- tbl_7 %>% 
  mutate(id = rownames(tbl_7)) %>% 
  gather('pop', 'prob', V01:V07) %>% 
  group_by(id) %>% 
  mutate(likely_assignment = pop[which.max(prob)],
         assingment_prob = max(prob)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))
plot_data_7$id <- factor(plot_data_7$id,levels = c("PA331", "PC622", "PC621", "PC623", "PA329", "PA332", "PA334", "PC580", "PC581", "PC616", "PC617", "PC618", "PA040", "PA047", "PA049", "PA341", "PC587", "PC592", "PC614", "PC631-r", "PA354", "PC630", "PA255", "PA256", "PA343", "PA346", "PC582", "PA363", "PC578", "PC579", "PC586", "PC593"))

plot_data_8 <- tbl_8 %>% 
  mutate(id = rownames(tbl_8)) %>% 
  gather('pop', 'prob', V01:V08) %>% 
  group_by(id) %>% 
  mutate(likely_assignment = pop[which.max(prob)],
         assingment_prob = max(prob)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))
plot_data_8$id <- factor(plot_data_8$id,levels = c("PA331", "PC622", "PC621", "PC623", "PA329", "PA332", "PA334", "PC580", "PC581", "PC616", "PC617", "PC618", "PA040", "PA047", "PA049", "PA341", "PC587", "PC592", "PC614", "PC631-r", "PA354", "PC630", "PA255", "PA256", "PA343", "PA346", "PC582", "PA363", "PC578", "PC579", "PC586", "PC593"))

plot_data_9 <- tbl_9 %>% 
  mutate(id = rownames(tbl_9)) %>% 
  gather('pop', 'prob', V01:V09) %>% 
  group_by(id) %>% 
  mutate(likely_assignment = pop[which.max(prob)],
         assingment_prob = max(prob)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))
plot_data_9$id <- factor(plot_data_9$id,levels = c("PA331", "PC622", "PC621", "PC623", "PA329", "PA332", "PA334", "PC580", "PC581", "PC616", "PC617", "PC618", "PA040", "PA047", "PA049", "PA341", "PC587", "PC592", "PC614", "PC631-r", "PA354", "PC630", "PA255", "PA256", "PA343", "PA346", "PC582", "PA363", "PC578", "PC579", "PC586", "PC593"))

plot_data_10 <- tbl_10 %>% 
  mutate(id = rownames(tbl_10)) %>% 
  gather('pop', 'prob', V01:V10) %>% 
  group_by(id) %>% 
  mutate(likely_assignment = pop[which.max(prob)],
         assingment_prob = max(prob)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))
plot_data_10$id <- factor(plot_data_10$id,levels = c("PA331", "PC622", "PC621", "PC623", "PA329", "PA332", "PA334", "PC580", "PC581", "PC616", "PC617", "PC618", "PA040", "PA047", "PA049", "PA341", "PC587", "PC592", "PC614", "PC631-r", "PA354", "PC630", "PA255", "PA256", "PA343", "PA346", "PC582", "PA363", "PC578", "PC579", "PC586", "PC593"))

q2<- ggplot(plot_data_2, aes(id, prob, fill = pop)) +
  geom_col() +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_fill_brewer(palette='Paired') +
  theme_classic()
q3<- ggplot(plot_data_3, aes(id, prob, fill = pop)) +
  geom_col() +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_fill_brewer(palette='Paired') +
  theme_classic()
q4<- ggplot(plot_data_4, aes(id, prob, fill = pop)) +
  geom_col() +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_fill_brewer(palette='Paired') +
  theme_classic()
q5<- ggplot(plot_data_5, aes(id, prob, fill = pop)) +
  geom_col() +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_fill_brewer(palette='Paired') +
  theme_classic()
q6<- ggplot(plot_data_6, aes(id, prob, fill = pop)) +
  geom_col() +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_fill_brewer(palette='Paired') +
  theme_classic()
q7<- ggplot(plot_data_7, aes(id, prob, fill = pop)) +
  geom_col() +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_fill_brewer(palette='Paired') +
  theme_classic()
q8<- ggplot(plot_data_8, aes(id, prob, fill = pop)) +
  geom_col() +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_fill_brewer(palette='Paired') +
  theme_classic()
q9<- ggplot(plot_data_9, aes(id, prob, fill = pop)) +
  geom_col() +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_fill_brewer(palette='Paired') +
  theme_classic()
q10<- ggplot(plot_data_10, aes(id, prob, fill = pop)) +
  geom_col() +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_fill_brewer(palette='Paired') +
  theme_classic()
grid.arrange(q2,q3,q4,q5,q6,q7,q8,q9,q10,nrow=9)


#### FINERADSTRUCTURE ####
### 1) EDIT THE FOLLOWING THREE LINES TO PROVIDE PATHS TO THE fineRADstructure OUTPUT 
setwd("03_fineradstructure/") ## The directory where the files are located
chunkfile<-"proteus_haplotypes_filtered_reordered_chunks.out" ## RADpainter output file
mcmcfile<-"proteus_haplotypes_filtered_reordered_chunks.mcmc.xml" ## finestructure mcmc file
treefile<-"proteus_haplotypes_filtered_reordered_chunks.mcmcTree.xml" ## finestructure tree file
### 2) EDIT THIS PATH TO WHERE YOU WANT THE PLOTS:
plotsFolder <- "03_fineradstructure/"
### 3) SET VALUES FOR THESE VARIABLES: "analysisName" will be included in output plots
analysisName <- "Proteus population structure";  maxIndv <- 10000; maxPop<-10000
### 4) EDIT THE PATH TO YOUR COPY of FinestructureLibrary.R
source("C:/Users/hans_/Desktop/Projects/Asellus_cave_evolution/04_analysis/02_fineradstructure/FinestructureLibrary.R", chdir = TRUE) # read in the R functions, which also calls the needed packages
### 5) EXECUTE THE CODE ABOVE AND THE REST OF THE CODE BELOW
## make some colours
some.colors<-MakeColorYRP() # these are yellow-red-purple
some.colorsEnd<-MakeColorYRP(final=c(0.2,0.2,0.2)) # as above, but with a dark grey final for capped values
### READ IN THE CHUNKCOUNT FILE
dataraw<-as.matrix(read.table(chunkfile,row.names=1,header=T,skip=1)) # read in the pairwise coincidence 
### READ IN THE MCMC FILES
mcmcxml<-xmlTreeParse(mcmcfile) ## read into xml format
mcmcdata<-as.data.frame.myres(mcmcxml) ## convert this into a data frame
### READ IN THE TREE FILES
treexml<-xmlTreeParse(treefile) ## read the tree as xml format
ttree<-extractTree(treexml) ## extract the tree into ape's phylo format
## Reduce the amount of significant digits printed in the posteror assignment probabilities (numbers shown in the tree):
ttree$node.label[ttree$node.label!=""] <-format(as.numeric(ttree$node.label[ttree$node.label!=""]),digits=2)
# convert to dendrogram format
tdend<-myapetodend(ttree,factor=1)
## Now we work on the MAP state
mapstate<-extractValue(treexml,"Pop") # map state as a finestructure clustering
mapstatelist<-popAsList(mapstate) # .. and as a list of individuals in populations
popnames<-lapply(mapstatelist,NameSummary) # population names IN A REVERSIBLE FORMAT (I.E LOSSLESS)
## NOTE: if your population labels don't correspond to the format we used (NAME<number>) YOU MAY HAVE TROUBLE HERE. YOU MAY NEED TO RENAME THEM INTO THIS FORM AND DEFINE YOUR POPULATION NAMES IN popnamesplot BELOW
popnamesplot<-lapply(mapstatelist,NameMoreSummary) # a nicer summary of the populations
names(popnames)<-popnamesplot # for nicety only
names(popnamesplot)<-popnamesplot # for nicety only
popdend<-makemydend(tdend,mapstatelist) # use NameSummary to make popdend
popdend<-fixMidpointMembers(popdend) # needed for obscure dendrogram reasons
popdendclear<-makemydend(tdend,mapstatelist,"NameMoreSummary")# use NameMoreSummary to make popdend
popdendclear<-fixMidpointMembers(popdendclear) # needed for obscure dendrogram reasons
###
## Plot 1: COANCESTRY MATRIX
fullorder<-labels(tdend) # the order according to the tree
datamatrix<-dataraw[fullorder,fullorder] # reorder the data matrix
tmpmat<-datamatrix 
tmpmat[tmpmat>maxIndv]<-maxIndv #  # cap the heatmap
pdf(file=paste(plotsFolder,analysisName,"-SimpleCoancestry.pdf",sep=""),height=25,width=25)
plotFinestructure(tmpmat,dimnames(tmpmat)[[1]],dend=tdend,cols=some.colorsEnd,cex.axis=1.1,edgePar=list(p.lwd=0,t.srt=90,t.off=-0.1,t.cex=1.2))
dev.off()
###
## Plot 2: POPULATIONS AND COANCESTRY AVERAGES
popmeanmatrix<-getPopMeanMatrix(datamatrix,mapstatelist)
tmpmat<-popmeanmatrix
tmpmat[tmpmat>maxPop]<-maxPop # cap the heatmap
pdf(file=paste(plotsFolder,analysisName,"-PopAveragedCoancestry.pdf",sep=""),height=20,width=20)
plotFinestructure(tmpmat,dimnames(tmpmat)[[1]],dend=tdend,cols=some.colorsEnd,cex.axis=1.1,edgePar=list(p.lwd=0,t.srt=90,t.off=-0.1,t.cex=1.2))
dev.off()
###
## Plot 3: POPULATIONS AND COANCESTRY AVERAGES WITH PERHAPS MORE INFORMATIVE LABELS
mappopcorrectorder<-NameExpand(labels(popdend))
mappopsizes<-sapply(mappopcorrectorder,length)
labellocs<-PopCenters(mappopsizes)
xcrt=0
ycrt=45
pdf(file=paste(plotsFolder,analysisName,"-PopAveragedCoancestry2.pdf",sep=""),height=25,width=25)
plotFinestructure(tmpmat,dimnames(tmpmat)[[1]],labelsx=labels(popdendclear),labelsatx=labellocs,xcrt=xcrt,cols=some.colorsEnd,ycrt=ycrt,dend=tdend,cex.axis=1.1,edgePar=list(p.lwd=0,t.srt=90,t.off=-0.1,t.cex=1.2),hmmar=c(3,0,0,1))
dev.off()


#### PCA ####
## adegenet package
setwd("04_pca/")
all.proteus <- read.PLINK("p6_r66_single_snp_minmaf05_all.raw", parallel=FALSE)
istria.rem.proteus <- read.PLINK("p6_r66_single_snp_minmaf05_rem_istria.raw", parallel=FALSE)
slo.proteus <- read.PLINK("p4_r66_single_snp_minmaf05_slovenia.raw", parallel=FALSE)
slo.light.proteus <- read.PLINK("p2_r66_single_snp_minmaf05_slo_light.raw", parallel=FALSE)
cro.bos.proteus <- read.PLINK("p2_r66_single_snp_minmaf05_cro_bos.raw", parallel=FALSE)

pca1 <- glPca(all.proteus, nf = 11)
pca2 <- glPca(istria.rem.proteus, nf = 11)
pca3 <- glPca(slo.proteus, nf = 6)
pca4 <- glPca(cro.bos.proteus, nf = 5)
# simple plot to check
scatter(pca1, posi="bottomright")
add.scatter.eig(pca1$eig[1:10],2,3,4, posi="topright", inset=.05, ratio=.3)
s.label(pca1$scores,xax=3,yax=5)
add.scatter.eig(ca1$eig,nf=3,xax=2,yax=3,posi="topleft")
scatter(pca2, posi="bottomright")
add.scatter.eig(pca2$eig[1:10],2,1,2, posi="topright", inset=.05, ratio=.3)
scatter(pca3, posi="topleft")
add.scatter.eig(pca3$eig[1:10],2,1,2, posi="topright", inset=.05, ratio=.3)
scatter(pca4, posi="topright")
add.scatter.eig(pca4$eig[1:10],2,1,2, posi="topright", inset=.05, ratio=.3)
# colour plot
myCol <- colorplot(pca1$scores,pca1$scores, transp=TRUE, cex=4)
abline(h=0,v=0, col="grey")
add.scatter.eig(pca1$eig[1:10],2,1,2, posi="topright", inset=.08, ratio=.2)
myCol <- colorplot(pca1$scores,pca1$scores, axes=1,3,4, transp=TRUE, cex=4)
abline(h=0,v=0, col="grey")
add.scatter.eig(pca1$eig[1:10],2,3,4, posi="topright", inset=.08, ratio=.2)
myCol <- colorplot(pca2$scores,pca2$scores, transp=TRUE, cex=4)
abline(h=0,v=0, col="grey")
add.scatter.eig(pca$eig[1:10],2,1,2, posi="topleft", inset=.08, ratio=.2)
myCol <- colorplot(pca3$scores,pca3$scores, transp=TRUE, cex=4)
abline(h=0,v=0, col="grey")
add.scatter.eig(pca3$eig[1:10],2,1,2, posi="topleft", inset=.08, ratio=.2)
colorplot(pca4$scores,pca4$scores, transp=TRUE, cex=4)
abline(h=0,v=0, col="grey")
add.scatter.eig(pca4$eig[1:10],2,1,2, posi="topleft", inset=.08, ratio=.2)

### Fst correlation plots ###
setwd("05_fst/")
dat <- read.csv("Proteus_fst_vs_dxy_matrix.tsv", sep = "\t")
data <- read.csv("Proteus_fst_vs_dxy_matrix_2.tsv", header=FALSE, sep = "\t")
rows<-(dat$X)
cols<-(dat$X)
colnames(data)<-cols
rownames(data)<-rows
names(data) <- cols
fst<-as.matrix(data)
dxy<-as.matrix(data)
dxy[upper.tri(dxy)]=NA
fst[lower.tri(fst)]=NA
pheatmap(dxy, display_numbers = T, number_format='%.3f', color = colorRampPalette(c('goldenrod2','firebrick3'))(100), cluster_rows = F, cluster_cols = F, fontsize_number = 15)
pheatmap(fst, display_numbers = T, number_format="%.2f", color = colorRampPalette(c('goldenrod2','firebrick3'))(100), cluster_rows = F, cluster_cols = F, fontsize_number = 15)

#### TREEMIX ####
setwd("06_treemix/")
source("plotting_funcs.R")
plot_tree("proteus_treemix")
plot_tree("proteus_treemix_sz_corr")
plot_tree("proteus_treemix_m1")
plot_tree("proteus_treemix_m2")
plot_tree("proteus_treemix_m3")
plot_tree("proteus_treemix_m4")
plot_tree("proteus_treemix_m5")
plot_resid("proteus_treemix", "poporder.txt")
plot_resid("proteus_treemix_sz_corr", "poporder.txt")
plot_resid("proteus_treemix_m1", "poporder.txt")
plot_resid("proteus_treemix_m2", "poporder.txt")
plot_resid("proteus_treemix_m3", "poporder.txt")
treemix_res <- optM("C:/Users/hans_/Desktop/Projects/Proteus_phylogeny/04_analysis/05_treemix/02_results/", method = "Evanno", thresh = 0.05)
plot_optM(treemix_res, method = "Evanno", plot = TRUE, pdf = "evanno_treemix.pdf")

#### BPP ####
#create random sample of 100 loci from the 29145 loci
ran_samples_100 <- sample(x = 1:29145, size = 100)
# load mcmc file and convert to data matrix
mcmc <- read.table("proteus_100_loci_mcmc.txt",header=TRUE)
mcmc_proteus <- as.data.frame(mcmc)
# load tree
#tree <- readMCMCtree("FigTree.tre")
tree <- read.tree("species_tree.tre")

mcmc.summary(mcmc, prob = 0.95)
test1 <- curve(bppr::dslnorm(x, shift=14, meanlog=0, sdlog=0.01), from=14, to=20, n=1000)
str(test1)
max(test1$y)
root_calmsc <- msc2time.t(mcmc=mcmc_proteus, 
                     node="10Para_litKrajinaLikaKras_CarsoSW_SloSticnaSE_SloparkeljIstria", 
                     calf=rslnorm, shift=15.8, meanlog=0.01, sdlog=0.1)
mcmc.summary(root_calmsc)
root_RAD_BSC2 <- msc2time.t(mcmc=mcmc_proteus, 
                          node="10Para_litKrajinaLikaKras_CarsoSW_SloSticnaSE_SloparkeljIstria", 
                          calf=rslnorm, shift=15.4, meanlog=0.1, sdlog=0.1)
mcmc.summary(root_RAD_BSC2)
root_mtDNA_starbeast2 <- msc2time.t(mcmc=mcmc_proteus, 
                          node="10Para_litKrajinaLikaKras_CarsoSW_SloSticnaSE_SloparkeljIstria", 
                          calf=rslnorm, shift=14.7, meanlog=0.1, sdlog=0.1)
mcmc.summary(root_mtDNA_starbeast2)

mcmc2densitree(tree, root_calmsc,"t_", thin=0.05, alpha=0.01)

root_calmsc_uni <- msc2time.t(mcmc=mcmc_proteus, 
                          node="10Para_litKrajinaLikaKras_CarsoSW_SloSticnaSE_SloparkeljIstria", 
                          calf=runif, shift=8.4, max=20.2)

plot(density(root_calmsc$t_10Para_litKrajinaLikaKras_CarsoSW_SloSticnaSE_SloparkeljIstria, adj=.1), xlab="Time (Ma)",
     main="root age")
rug(root_calmsc$t_10Para_litKrajinaLikaKras_CarsoSW_SloSticnaSE_SloparkeljIstria)

mcmc2densitree(tree, root_calmsc,"t_", thin=0.05, alpha=0.01)
title(xlab="Divergence time (Ma)")


mcmc.summary(root_calmsc)
apply(root_calmsc, 2, mean)
coda::HPDinterval(coda::as.mcmc(root_calmsc))
curve(bppr::dslnorm(x, shift=3, meanlog=0.5, sdlog=0.5), from=0, to=15, n=1e3)
root_uplift_calmsc_ln <- msc2time.t(mcmc=mcmc_proteus, 
                     node="14LikaKras_Carso", 
                     calf=runif, shift=3,  meanlog=0.5, sdlog=0.5)
root_uplift_calmsc_nd <- msc2time.t(mcmc=mcmc_proteus, 
                                 node="14LikaKras_Carso", 
                                 calf=runif, min=3, max=5)
str(root_calmsc_uni)
mcmc.summary(root_uplift_calmsc_ln)
mcmc.summary(root_uplift_calmsc_nd)
apply(root_uplift_calmsc_nd, 2, mean)
mcmc2densitree(tree, root_uplift_calmsc,"t_", thin=0.05, alpha=0.01)


### MCMCtree ###

## shortest HPD sum tree
setwd("08_MCMCtree/01_root_LK_BSC2")
phy <- readMCMCtree("FigTree.tre")
mcmc <-read.table("mcmc.txt",header=TRUE)
MCMCtree.posterior <- as.data.frame(mcmc)
MCMC.tree.plot(phy, build.tree =TRUE, MCMC.chain = MCMCtree.posterior, cex.tips = 1.5, 
               time.correction = 100, plot.type = "distributions", cex.age = 1.5, 
               cex.labels = 1.5, col.tree = "grey40", 
               scale.res = c("Epoch"), add.time.scale	= TRUE, relative.height = 0.08, add.abs.time =TRUE,
               density.col = "#00000050", density.border.col = "#00000080")

## only Lika split and root
setwd("08_MCMCtree/01_LK")
phy <- readMCMCtree("FigTree.tre")
#phy <- read.tree("timetree_abs_ages.tre")
mcmc <-read.table("mcmc.txt",header=TRUE)
MCMCtree.posterior <- as.data.frame(mcmc)
node.post <-read.table("node_posteriors.txt",header=TRUE,check.names = FALSE)
node.posterior <- as.data.frame(node.post)
str(node.posterior)
MCMC.tree.plot(phy, build.tree =TRUE, MCMC.chain = MCMCtree.posterior, cex.tips = 1.5, 
               time.correction = 100, plot.type = "distributions", cex.age = 1.5, 
               cex.labels = 1.5, col.tree = "grey40", 
               scale.res = c("Epoch"), add.time.scale	= TRUE, relative.height = 0.08, add.abs.time =TRUE,
               density.col = "#00000050", density.border.col = "#00000080")

## mitochondrial DNA tree
setwd("08_MCMCtree/03_mtDNA")
phy <- readMCMCtree("FigTree.tre")
mcmc <-read.table("mcmc.txt",header=TRUE)
MCMCtree.posterior <- as.data.frame(mcmc)
MCMC.tree.plot(phy, build.tree =TRUE, MCMC.chain = MCMCtree.posterior, cex.tips = 1.5, 
               time.correction = 100, plot.type = "distributions", cex.age = 1.5, 
               cex.labels = 1.5, col.tree = "grey40", 
               scale.res = c("Epoch"), add.time.scale	= TRUE, relative.height = 0.08, add.abs.time =TRUE,
               density.col = "#00000050", density.border.col = "#00000080")


########### MtDNA analysis ########### 
### Plot different indices of diversity against samples and locations
setwd("08_arlequin")
moldiv <- read.csv("proteus_mol_diversity2.csv", header=TRUE)
str(moldiv)
plot(moldiv)
moldiv <- read.csv("area_vs_diversity.csv", header=TRUE)
str(moldiv)
plot(Pi~log.area., data=moldiv, pch=20, cex=3)
abline(lm(Pi~log.area., data=moldiv))
summary(lm(Pi~log.area., data=moldiv))
plot(Pi_genomics~log.area., data=moldiv, pch=20, cex=3)
abline(lm(Pi_genomics~log.area., data=moldiv))
summary(lm(Pi_genomics~log.area., data=moldiv))


### Phylogenetic analyses of mtDNA
### MAXIMUM PARSIMONY ###
## read in the trees
setwd("C:/Users/hans_/Desktop/Projects/Proteus_phylogeny/04_analysis/09_mtDNA/04_MP")
jack <- read.nexus("MP_tree_500_jackknife2")
boot<- read.nexus("MP_tree_500_bootstraps")
plot(jack)
plot(boot)
collapsed = collapse_tree_at_resolution(boot, resolution=380,shorten=FALSE,rename_collapsed_nodes = TRUE)$tree
plot(collapsed)
coll.tree<- collapse_tree_at_resolution(jack, 
                            resolution             = 0, 
                            by_edge_count          = FALSE,
                            shorten                = TRUE,
                            rename_collapsed_nodes = FALSE,
                            criterion              = 'max_tip_depth')
plot(coll.tree)


### ANCESTRAL TRAIT RECONSTRUCTION ###
## read in the trees
#corHMM
setwd("11_ancestral_state_reconstruction/")
morphs <- read.table("Proteus_states_ances.txt", header=TRUE)
str(morphs)
tree=read.tree("RAxML_bestTree.proteus_all_1-3snps_r66_MD_root_ancestor.tre")
# ancestral trait reconstruction
# null model
ARD_morph_fixed=rayDISC(tree,morphs,ntraits=1,charnum=1,model="ARD",node.states="marginal",root.p=c(1,0))
# equal rates
ARD_morph_equal=rayDISC(tree,morphs,ntraits=1,charnum=1,model="ER",node.states="marginal",root.p=c(1,0))
# no reversal
ARD_Dollo=rayDISC(tree,morphs,ntraits=1,charnum=1, p=c(0,1), model="ARD",node.states="marginal",root.p=c(1,0))
# no reversal equal rates
ARD_Dollo_equal=rayDISC(tree,morphs,ntraits=1,charnum=1, p=c(0,1), model="ER",node.states="marginal",root.p=c(1,0))
# no reversal equal rates 2
ARD_Dollo_equal=rayDISC(tree,morphs,ntraits=1,charnum=1, p=c(0,1), root.p=c(1,0))

plotRECON(tree,ARD_Dollo_equal$states,piecolors=col.list,
          cex=1.1, pie.cex=0.7)
ARD_Dollo_equal$tip.states

# null model rates
ASR_morph_rates=rayDISC(tree,morphs,ntraits=1,charnum=1, p=c(0.709, 1.892),model="ARD",node.states="marginal",root.p=c(1,0))

round(ARD_morph_fixed$states,2)
round(ARD_morph_fixed$states,2)

## plot the trees
col.list<-c("#52453A","#DDAF92")
plotRECON(tree,ARD_morph_equal$states,piecolors=col.list,
          cex=1.1, pie.cex=0.7,height=9, width=9, file="ER2.pdf")
plotRECON(tree,ARD_Dollo$states,piecolors=col.list,
          cex=1.1, pie.cex=0.7,height=9, width=9, file="Dollo-test2.pdf")
plotRECON(tree,ARD_morph_fixed$states,piecolors=col.list,
          cex=1.1, pie.cex=0.7)
plotRECON(tree,ARD_morph_equal$states,piecolors=col.list,
          cex=1.1, pie.cex=0.7)
plotRECON(tree,ASR_morph_rates$states,piecolors=col.list,
          cex=1.1, pie.cex=0.7)
plotRECON(tree,ARD_Dollo_equal$states,piecolors=col.list,
          cex=1.1, pie.cex=0.7)
plotRECON(tree,ARD_Dollo$states,piecolors=col.list,
          cex=1.1, pie.cex=0.7)


# calculate AICs
str(ARD_morph_fixed)
str(ARD_Dollo)
ARD_morph_fixed$AICc
ARD_morph_equal$AICc
ARD_Dollo$AICc
ARD_Dollo_equal$AICc
ASR_morph_rates$AICc
lr.test(ARD_morph_equal$loglik,ARD_Dollo$loglik)
ARD_morph_fixed$solution
ARD_Dollo$solution

### use only species as tips
setwd("11_ancestral_state_reconstruction/")
morphs <- read.table("Proteus_states_ances_species.txt", header=TRUE)
str(morphs)
tree=read.tree("species_tree.tre")
## start ancestral trait reconstruction
# null model
ARD_morph_fixed=rayDISC(tree,morphs,ntraits=1,charnum=1,model="ARD",node.states="marginal",root.p=c(1,0))
# equal rates
ARD_morph_equal=rayDISC(tree,morphs,ntraits=1,charnum=1,model="ER",node.states="marginal",root.p=c(1,0))
# no reversal
ARD_Dollo=rayDISC(tree,morphs,ntraits=1,charnum=1, p=c(0,1), model="ARD",node.states="marginal",root.p=c(1,0))
# no reversal equal rates
ARD_Dollo_equal=rayDISC(tree,morphs,ntraits=1,charnum=1, p=c(0,1), model="ER",node.states="marginal",root.p=c(1,0))
# no reversal equal rates 2
ARD_Dollo_equal=rayDISC(tree,morphs,ntraits=1,charnum=1, p=c(0,1), root.p=c(1,0))
plotRECON(tree,ARD_Dollo_equal$states,piecolors=col.list,
          cex=1.1, pie.cex=0.7)
ARD_Dollo_equal$tip.states
# null model rates
ASR_morph_rates=rayDISC(tree,morphs,ntraits=1,charnum=1, p=c(0.00001, 0.5),model="ARD",node.states="marginal",root.p=c(1,0))
round(ARD_morph_fixed$states,2)
round(ARD_morph_fixed$states,2)
plotRECON(tree,ASR_morph_rates$states,piecolors=col.list,
          cex=1.1, pie.cex=0.7)
ASR_morph_rates$AICc

## plot the trees
col.list<-c("#52453A","#DDAF92")
plotRECON(tree,ARD_morph_equal$states,piecolors=col.list,
          cex=1.1, pie.cex=0.7,height=9, width=9, file="species_tree_ER2.pdf")
plotRECON(tree,ARD_Dollo$states,piecolors=col.list,
          cex=1.1, pie.cex=0.7,height=9, width=9, file="species_tree_Dollo-test2.pdf")
plotRECON(tree,ARD_morph_fixed$states,piecolors=col.list,
          cex=1.1, pie.cex=0.7)
plotRECON(tree,ARD_morph_equal$states,piecolors=col.list,
          cex=1.1, pie.cex=0.7)
plotRECON(tree,ASR_morph_rates$states,piecolors=col.list,
          cex=1.1, pie.cex=0.7)
plotRECON(tree,ARD_Dollo_equal$states,piecolors=col.list,
          cex=1.1, pie.cex=0.7)
plotRECON(tree,ARD_Dollo$states,piecolors=col.list,
          cex=1.1, pie.cex=0.7)
# calculate AICs
str(ARD_morph_fixed)
str(ARD_Dollo)
ARD_morph_fixed$AICc
ARD_morph_equal$AICc
ARD_Dollo$AICc
ARD_Dollo_equal$AICc
ASR_morph_rates$AICc
lr.test(ARD_morph_equal$loglik,ARD_Dollo$loglik)
ARD_morph_fixed$solution
ARD_Dollo$solution



### SPECIES DELIMITATION PLOTTING ###

## read in the trees
library(rmutil)
setwd("12_species_delimitation")
twelveS.tree=read.tree("UPGMA_tree_12S.nwk")
dloop.tree=read.tree("UPGMA_tree_dloop.nwk")
CO1.tree=read.tree("UPGMA_tree_CO1.nwk")
#CO1 tree has a taxon which has been excluded later (PA134), -> exclude this tip
CO1.tree.new <- drop.tip(CO1.tree, c("PA134"), trim.internal = TRUE, subtree = FALSE,
         root.edge = 0, rooted = is.rooted(CO1.tree))

## read in species delimitations
species_12S <- read.csv("species_12S.csv", header = T,stringsAsFactors = TRUE)
species_12S_m <- species_12S[,-1]
rownames(species_12S_m) <- species_12S[,1]
species_12S_m
species_dloop <- read.csv("species_dloop.csv", header = T,stringsAsFactors = TRUE)
species_dloop_m <- species_dloop[,-1]
rownames(species_dloop_m) <- species_dloop[,1]
species_dloop_m
species_CO1 <- read.csv("species_CO1.csv", header = T,stringsAsFactors = TRUE)
species_CO1_m <- species_CO1[,-1]
rownames(species_CO1_m) <- species_CO1[,1]
species_CO1_m

# plot species delimitations on trees
trait.plot(twelveS.tree, species_12S_m, cols = list(ABGD = c("#E3191C", "#693D99","#2F9E2A","#F89896","darkolivegreen4","#FE9418","#E3D2FF","#1F77B3","#B2DE89"), 
                                                    ASAP = c("#E3191C", "#693D99","#2F9E2A","#1F77B3","#F89896","#B2DE89","#FE9418")))
trait.plot(dloop.tree, species_dloop_m, cols = list(ABGD = c("#E3191C", "#693D99","#2F9E2A","#F89896","#B2DE89","#FE9418","#E3D2FF","#1F77B3"), 
                                                    ASAP = c("#E3191C", "#693D99","#2F9E2A","#F89896","tan4","#B2DE89","#FDBE6F","#FE9418","#E3D2FF","#1F77B3")))
trait.plot(CO1.tree.new, species_CO1_m, cols = list(ABGD = c("#E3191C", "#693D99","#2F9E2A","#F89896","#B2DE89","#FE9418","#1F77B3"), 
                                                    ASAP = c("#E3191C", "#693D99","mediumorchid1","#2F9E2A","#F89896","darkolivegreen4","#B2DE89","seagreen1","#FDBE6F","#FE9418","brown2","yellow2","#E3D2FF","#1F77B3")))

### Add tree based analyses
# read in tree
beast.tree=read.tree("proteus_constant_pop.tre")
newtree <- drop.tip(beast.tree, c("PA766"), trim.internal = TRUE, subtree = FALSE,
                         root.edge = 0, rooted = is.rooted(beast.tree))
str(newtree)
newtree2 <- keep.tip(newtree, c("PA331","PC623","PA300","PA263","PC540", "PA816","PA260",
                                "PA315", "PC617","PA276","PA297","PA157","PA321","PA348","PA033","PA255",
                                "PB600","PA011","PA004"))
plot(newtree2)
## create input files separately for each gene
# GMYC & PTP results
results <- read.csv("results_GMYC_PTP_no_lin.csv", header=TRUE)
str(results)
results_m <- results[,-1]
rownames(results_m) <- results[,1]
results_m

#all.results <- merge(all.data, results, by = 'ID',all.y=TRUE)
#write.table(all.results,"all.results.txt",sep="\t")
# plot species delimitations on trees
trait.plot(newtree, results_m, cols = list(BPTP = c("#E3191C", "#693D99","mediumorchid1","#2F9E2A","#F89896","darkolivegreen4","darkolivegreen3","#B2DE89","seagreen1","#FDBE6F","#FE9418","brown2","yellow2","#E3D2FF","#1F77B3"), 
                                           PTP = c("#E3191C", "#693D99","mediumorchid1","#2F9E2A","#F89896","darkolivegreen4","darkolivegreen3","#B2DE89","seagreen1","#FDBE6F","#FE9418","brown2","yellow2","#E3D2FF","#1F77B3"),
                                           GMYC = c("#E3191C", "#693D99","mediumorchid1","#2F9E2A","#F89896","darkolivegreen4","darkolivegreen3","#B2DE89","seagreen1","#FDBE6F","#FE9418","brown2","yellow2","#E3D2FF","#1F77B3")))

### SPATIAL ANALYSIS ###
setwd("13_spatial_analysis")
### Perform analysis separately based on lineages
### CO1 gene
## Read sequence data in fasta
Dolenjska.CO1 <- read.dna(file="Dolenjska_CO1.fas", format="fasta")
Ljubljanica.CO1 <- read.dna(file="Ljubljanica_CO1.fas", format="fasta")
Kras.CO1 <- read.dna(file="Kras_CO1.fas", format="fasta")
Paralittoral.CO1 <- read.dna(file="Paralittoral_CO1.fas", format="fasta")
Lika.CO1 <- read.dna(file="Lika_CO1.fas", format="fasta")
parkelj.CO1 <- read.dna(file="parkelj_CO1.fas", format="fasta")
Sticna.CO1 <- read.dna(file="Sticna_CO1.fas", format="fasta")
Krajina.CO1 <- read.dna(file="Krajina_CO1.fas", format="fasta")
Istria.CO1 <- read.dna(file="Istria_CO1.fas", format="fasta")
class(Dolenjska.CO1)
# combine data
Dolenjska.CO1.genind <- DNAbin2genind(x=Dolenjska.CO1)
Ljubljanica.CO1.genind <- DNAbin2genind(x=Ljubljanica.CO1)
Kras.CO1.genind <- DNAbin2genind(x=Kras.CO1)
Paralittoral.CO1.genind <- DNAbin2genind(x=Paralittoral.CO1)
Lika.CO1.genind <- DNAbin2genind(x=Lika.CO1)
parkelj.CO1.genind <- DNAbin2genind(x=parkelj.CO1)
Krajina.CO1.genind <- DNAbin2genind(x=Krajina.CO1)
Sticna.CO1.genind <- DNAbin2genind(x=Sticna.CO1)
Istria.CO1.genind <- DNAbin2genind(x=Istria.CO1)
# add coordinates
Dolenjska.CO1.coord <- read.csv("Dolenjska_CO1_coordinates.csv", header=TRUE,quote="", dec=".", row.names=1)
Ljubljanica.CO1.coord <- read.csv("Ljubljanica_CO1_coordinates.csv", header=TRUE,quote="", dec=".", row.names=1)
Kras.CO1.coord <- read.csv("Kras_CO1_coordinates.csv", header=TRUE,quote="", dec=".", row.names=1)
Paralittoral.CO1.coord <- read.csv("Paralittoral_CO1_coordinates.csv", header=TRUE,quote="", dec=".", row.names=1)
Lika.CO1.coord <- read.csv("Lika_CO1_coordinates.csv", header=TRUE,quote="", dec=".", row.names=1)
parkelj.CO1.coord <- read.csv("parkelj_CO1_coordinates.csv", header=TRUE,quote="", dec=".", row.names=1)
Sticna.CO1.coord <- read.csv("Sticna_CO1_coordinates.csv", header=TRUE,quote="", dec=".", row.names=1)
Krajina.CO1.coord <- read.csv("Krajina_CO1_coordinates.csv", header=TRUE,quote="", dec=".", row.names=1)
Istria.CO1.coord <- read.csv("Istria_CO1_coordinates.csv", header=TRUE,quote="", dec=".", row.names=1)
# Add coordinates - note identification of slots within object
Dolenjska.CO1.genind$other$xy <- Dolenjska.CO1.coord
Ljubljanica.CO1.genind$other$xy <- Ljubljanica.CO1.coord
Kras.CO1.genind$other$xy <- Kras.CO1.coord
Paralittoral.CO1.genind$other$xy <- Paralittoral.CO1.coord
Lika.CO1.genind$other$xy <- Lika.CO1.coord
parkelj.CO1.genind$other$xy <- parkelj.CO1.coord
Krajina.CO1.genind$other$xy <- Krajina.CO1.coord
Sticna.CO1.genind$other$xy <- Sticna.CO1.coord
Istria.CO1.genind$other$xy <- Istria.CO1.coord
# Geographical distance
Dolenjska.CO1.gdist <- dist(x=Dolenjska.CO1.genind$other$xy, method="euclidean",diag=TRUE, upper=TRUE)
Ljubljanica.CO1.gdist <- dist(x=Ljubljanica.CO1.genind$other$xy, method="euclidean",diag=TRUE, upper=TRUE)
Kras.CO1.gdist <- dist(x=Kras.CO1.genind$other$xy, method="euclidean",diag=TRUE, upper=TRUE)
Paralittoral.CO1.gdist <- dist(x=Paralittoral.CO1.genind$other$xy, method="euclidean",diag=TRUE, upper=TRUE)
Lika.CO1.gdist <- dist(x=Lika.CO1.genind$other$xy, method="euclidean",diag=TRUE, upper=TRUE)
parkelj.CO1.gdist <- dist(x=parkelj.CO1.genind$other$xy, method="euclidean",diag=TRUE, upper=TRUE)
Krajina.CO1.gdist <- dist(x=Krajina.CO1.genind$other$xy, method="euclidean",diag=TRUE, upper=TRUE)
Sticna.CO1.gdist <- dist(x=Sticna.CO1.genind$other$xy, method="euclidean",diag=TRUE, upper=TRUE)
Istria.CO1.gdist <- dist(x=Istria.CO1.genind$other$xy, method="euclidean",diag=TRUE, upper=TRUE)
# Euclidean distance for individuals (plain ordinary distance matrix)
Dolenjska.CO1.dist <- dist(x=Dolenjska.CO1.genind, method="euclidean", diag=T, upper=T)
Ljubljanica.CO1.dist <- dist(x=Ljubljanica.CO1.genind, method="euclidean", diag=T, upper=T)
Kras.CO1.dist <- dist(x=Kras.CO1.genind, method="euclidean", diag=T, upper=T)
Paralittoral.CO1.dist <- dist(x=Paralittoral.CO1.genind, method="euclidean", diag=T, upper=T)
Lika.CO1.dist <- dist(x=Lika.CO1.genind, method="euclidean", diag=T, upper=T)
parkelj.CO1.dist <- dist(x=parkelj.CO1.genind, method="euclidean", diag=T, upper=T)
Krajina.CO1.dist <- dist(x=Krajina.CO1.genind, method="euclidean", diag=T, upper=T)
Sticna.CO1.dist <- dist(x=Sticna.CO1.genind, method="euclidean", diag=T, upper=T)
Istria.CO1.dist <- dist(x=Istria.CO1.genind, method="euclidean", diag=T, upper=T)
# Mantel test
Dolenjska.CO1.mantel <- mantel.randtest(m1=Dolenjska.CO1.dist, m2=Dolenjska.CO1.gdist,nrepet=1000)
Ljubljanica.CO1.mantel <- mantel.randtest(m1=Ljubljanica.CO1.dist, m2=Ljubljanica.CO1.gdist,nrepet=1000)
Kras.CO1.mantel <- mantel.randtest(m1=Kras.CO1.dist, m2=Kras.CO1.gdist,nrepet=1000)
Paralittoral.CO1.mantel <- mantel.randtest(m1=Paralittoral.CO1.dist, m2=Paralittoral.CO1.gdist,nrepet=1000)
Lika.CO1.mantel <- mantel.randtest(m1=Lika.CO1.dist, m2=Lika.CO1.gdist,nrepet=1000)
parkelj.CO1.mantel <- mantel.randtest(m1=parkelj.CO1.dist, m2=parkelj.CO1.gdist,nrepet=1000)
Krajina.CO1.mantel <- mantel.randtest(m1=Krajina.CO1.dist, m2=Krajina.CO1.gdist,nrepet=1000)
Sticna.CO1.mantel <- mantel.randtest(m1=Sticna.CO1.dist, m2=Sticna.CO1.gdist,nrepet=1000)
Istria.CO1.mantel <- mantel.randtest(m1=Istria.CO1.dist, m2=Istria.CO1.gdist,nrepet=1000)
plot(Dolenjska.CO1.mantel, nclass=30)
plot(Ljubljanica.CO1.mantel, nclass=30)
plot(Kras.CO1.mantel, nclass=30)
plot(Paralittoral.CO1.mantel, nclass=30)
plot(Lika.CO1.mantel, nclass=30)
plot(parkelj.CO1.mantel, nclass=30)
plot(Krajina.CO1.mantel, nclass=30)
plot(Sticna.CO1.mantel, nclass=30)
plot(Istria.CO1.mantel, nclass=30)

### Dloop gene
## Read sequence data in fasta
Dolenjska.Dloop <- read.dna(file="Dolenjska_Dloop.fas", format="fasta")
Ljubljanica.Dloop <- read.dna(file="Ljubljanica_Dloop.fas", format="fasta")
Kras.Dloop <- read.dna(file="Kras_Dloop.fas", format="fasta")
Paralittoral.Dloop <- read.dna(file="Paralittoral_Dloop.fas", format="fasta")
Lika.Dloop <- read.dna(file="Lika_Dloop.fas", format="fasta")
parkelj.Dloop <- read.dna(file="parkelj_Dloop.fas", format="fasta")
Sticna.Dloop <- read.dna(file="Sticna_Dloop.fas", format="fasta")
Krajina.Dloop <- read.dna(file="Krajina_Dloop.fas", format="fasta")
Istria.Dloop <- read.dna(file="Istria_Dloop.fas", format="fasta")
class(Dolenjska.Dloop)
# combine data
Dolenjska.Dloop.genind <- DNAbin2genind(x=Dolenjska.Dloop)
Ljubljanica.Dloop.genind <- DNAbin2genind(x=Ljubljanica.Dloop)
Kras.Dloop.genind <- DNAbin2genind(x=Kras.Dloop)
Paralittoral.Dloop.genind <- DNAbin2genind(x=Paralittoral.Dloop)
Lika.Dloop.genind <- DNAbin2genind(x=Lika.Dloop)
parkelj.Dloop.genind <- DNAbin2genind(x=parkelj.Dloop)
Krajina.Dloop.genind <- DNAbin2genind(x=Krajina.Dloop)
Sticna.Dloop.genind <- DNAbin2genind(x=Sticna.Dloop)
Istria.Dloop.genind <- DNAbin2genind(x=Istria.Dloop)
# add coordinates
Dolenjska.Dloop.coord <- read.csv("Dolenjska_Dloop_coordinates.csv", header=TRUE,quote="", dec=".", row.names=1)
Ljubljanica.Dloop.coord <- read.csv("Ljubljanica_Dloop_coordinates.csv", header=TRUE,quote="", dec=".", row.names=1)
Kras.Dloop.coord <- read.csv("Kras_Dloop_coordinates.csv", header=TRUE,quote="", dec=".", row.names=1)
Paralittoral.Dloop.coord <- read.csv("Paralittoral_Dloop_coordinates.csv", header=TRUE,quote="", dec=".", row.names=1)
Lika.Dloop.coord <- read.csv("Lika_Dloop_coordinates.csv", header=TRUE,quote="", dec=".", row.names=1)
parkelj.Dloop.coord <- read.csv("parkelj_Dloop_coordinates.csv", header=TRUE,quote="", dec=".", row.names=1)
Sticna.Dloop.coord <- read.csv("Sticna_Dloop_coordinates.csv", header=TRUE,quote="", dec=".", row.names=1)
Krajina.Dloop.coord <- read.csv("Krajina_Dloop_coordinates.csv", header=TRUE,quote="", dec=".", row.names=1)
Istria.Dloop.coord <- read.csv("Istria_Dloop_coordinates.csv", header=TRUE,quote="", dec=".", row.names=1)
# Add coordinates - note identification of slots within object
Dolenjska.Dloop.genind$other$xy <- Dolenjska.Dloop.coord
Ljubljanica.Dloop.genind$other$xy <- Ljubljanica.Dloop.coord
Kras.Dloop.genind$other$xy <- Kras.Dloop.coord
Paralittoral.Dloop.genind$other$xy <- Paralittoral.Dloop.coord
Lika.Dloop.genind$other$xy <- Lika.Dloop.coord
parkelj.Dloop.genind$other$xy <- parkelj.Dloop.coord
Krajina.Dloop.genind$other$xy <- Krajina.Dloop.coord
Sticna.Dloop.genind$other$xy <- Sticna.Dloop.coord
Istria.Dloop.genind$other$xy <- Istria.Dloop.coord
# Geographical distance
Dolenjska.Dloop.gdist <- dist(x=Dolenjska.Dloop.genind$other$xy, method="euclidean",diag=TRUE, upper=TRUE)
Ljubljanica.Dloop.gdist <- dist(x=Ljubljanica.Dloop.genind$other$xy, method="euclidean",diag=TRUE, upper=TRUE)
Kras.Dloop.gdist <- dist(x=Kras.Dloop.genind$other$xy, method="euclidean",diag=TRUE, upper=TRUE)
Paralittoral.Dloop.gdist <- dist(x=Paralittoral.Dloop.genind$other$xy, method="euclidean",diag=TRUE, upper=TRUE)
Lika.Dloop.gdist <- dist(x=Lika.Dloop.genind$other$xy, method="euclidean",diag=TRUE, upper=TRUE)
parkelj.Dloop.gdist <- dist(x=parkelj.Dloop.genind$other$xy, method="euclidean",diag=TRUE, upper=TRUE)
Krajina.Dloop.gdist <- dist(x=Krajina.Dloop.genind$other$xy, method="euclidean",diag=TRUE, upper=TRUE)
Sticna.Dloop.gdist <- dist(x=Sticna.Dloop.genind$other$xy, method="euclidean",diag=TRUE, upper=TRUE)
Istria.Dloop.gdist <- dist(x=Istria.Dloop.genind$other$xy, method="euclidean",diag=TRUE, upper=TRUE)
# Euclidean distance for individuals (plain ordinary distance matrix)
Dolenjska.Dloop.dist <- dist(x=Dolenjska.Dloop.genind, method="euclidean", diag=T, upper=T)
Ljubljanica.Dloop.dist <- dist(x=Ljubljanica.Dloop.genind, method="euclidean", diag=T, upper=T)
Kras.Dloop.dist <- dist(x=Kras.Dloop.genind, method="euclidean", diag=T, upper=T)
Paralittoral.Dloop.dist <- dist(x=Paralittoral.Dloop.genind, method="euclidean", diag=T, upper=T)
Lika.Dloop.dist <- dist(x=Lika.Dloop.genind, method="euclidean", diag=T, upper=T)
parkelj.Dloop.dist <- dist(x=parkelj.Dloop.genind, method="euclidean", diag=T, upper=T)
Krajina.Dloop.dist <- dist(x=Krajina.Dloop.genind, method="euclidean", diag=T, upper=T)
Sticna.Dloop.dist <- dist(x=Sticna.Dloop.genind, method="euclidean", diag=T, upper=T)
Istria.Dloop.dist <- dist(x=Istria.Dloop.genind, method="euclidean", diag=T, upper=T)
# Mantel test
Dolenjska.Dloop.mantel <- mantel.randtest(m1=Dolenjska.Dloop.dist, m2=Dolenjska.Dloop.gdist,nrepet=1000)
Ljubljanica.Dloop.mantel <- mantel.randtest(m1=Ljubljanica.Dloop.dist, m2=Ljubljanica.Dloop.gdist,nrepet=1000)
Kras.Dloop.mantel <- mantel.randtest(m1=Kras.Dloop.dist, m2=Kras.Dloop.gdist,nrepet=1000)
Paralittoral.Dloop.mantel <- mantel.randtest(m1=Paralittoral.Dloop.dist, m2=Paralittoral.Dloop.gdist,nrepet=1000)
Lika.Dloop.mantel <- mantel.randtest(m1=Lika.Dloop.dist, m2=Lika.Dloop.gdist,nrepet=1000)
parkelj.Dloop.mantel <- mantel.randtest(m1=parkelj.Dloop.dist, m2=parkelj.Dloop.gdist,nrepet=1000)
Krajina.Dloop.mantel <- mantel.randtest(m1=Krajina.Dloop.dist, m2=Krajina.Dloop.gdist,nrepet=1000)
Sticna.Dloop.mantel <- mantel.randtest(m1=Sticna.Dloop.dist, m2=Sticna.Dloop.gdist,nrepet=1000)
Istria.Dloop.mantel <- mantel.randtest(m1=Istria.Dloop.dist, m2=Istria.Dloop.gdist,nrepet=1000)
plot(Dolenjska.Dloop.mantel, nclass=30)
plot(Ljubljanica.Dloop.mantel, nclass=30)
plot(Kras.Dloop.mantel, nclass=30)
plot(Paralittoral.Dloop.mantel, nclass=30)
plot(Lika.Dloop.mantel, nclass=30)
plot(parkelj.Dloop.mantel, nclass=30)
plot(Krajina.Dloop.mantel, nclass=30)
plot(Sticna.Dloop.mantel, nclass=30)
plot(Istria.Dloop.mantel, nclass=30)
# test if dimensions are the same in case of errors
#dim(as.matrix(Ljubljanica.Dloop.dist))
#dim(as.matrix(Ljubljanica.Dloop.gdist))
## All Mantel tests with sufficient sample sizes:
Dolenjska.CO1.mantel
Ljubljanica.CO1.mantel
Kras.CO1.mantel
Paralittoral.CO1.mantel
Lika.CO1.mantel
Krajina.CO1.mantel
Sticna.CO1.mantel 

Dolenjska.Dloop.mantel
Ljubljanica.Dloop.mantel
Kras.Dloop.mantel
Paralittoral.Dloop.mantel
Lika.Dloop.mantel
Istria.Dloop.mantel 
=======
# load libraries
library(adegenet)
library(Amelia)
library(ape)
library(bppr)
library(caret)
library(castor)
library(chopper)
library(corHMM)
library(corrplot)
library(dartR)
library(diversitree)
library(extRemes)
library(factoextra)
library(FactoMineR)
library(ggfortify)
library(gridExtra)
library(lmtest)
library(MASS)
library(MCMCtreeR)
library(mice)
library(missForest)
library(missMDA)
library(naniar)
library(OptM)
library(pals)
library(pegas)
library(phangorn)
library(pheatmap)
library(phylotools)
library(plot.matrix)
library(poppr)
library(radiator)
library(RColorBrewer)
library(Rcpp)
library(reshape2)
library(seqinr)
library(smatr)
library(SNPRelate)
library(tidyverse)
library(VIM)
library(visdat)
#install_local("C:/Users/hans_/Downloads/radiator-master.zip")

########### GENOTYPING ########### 
### error rate estimation using three different files, unfiltered, filtered to reduce missing data, and using a whitelist
setwd("01_replicates")
proteus.error.wl <- read.PLINK("replicates_whitelist_1-3snps_r66.raw", parallel=FALSE)
source("C:/Users/hans_/Desktop/Projects/Proteus_phylogeny/04_analysis/08_replicates/PairsDiff.R")
source("C:/Users/hans_/Desktop/Projects/Proteus_phylogeny/04_analysis/08_replicates/SNPs_error.R") # read in the R functions, which also calls the needed packages
SNP_error(proteus.error.wl)
# transpose data to analyse in excel and find those loci that do not match
rep_256 <- read.table("256_replicates.txt", sep="\t", header=TRUE)
rep_592 <- read.table("592_replicates.txt", sep="\t", header=TRUE)
rep_592_t <- t(rep_592) 
#write.table(rep_256_t,"rep_256_t.txt")
#write.table(rep_592_t,"rep_592_t.txt")


########### POPULATION GENETICS ########### 
#### ADMIXTURE ####
setwd("02_admixture")
n <- 10
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))

colours.ten<-sample(col_vector, n)
par(mfrow=c(1,1))
data<-read.table("proteus_p1_r66_minmaf05_locus_chr.8.Q")
barplot(t(as.matrix(data)), col=colours.ten,
        xlab="Individual #", ylab="Ancestry", border=NA,las=2)
barplot(t(as.matrix(data)), col=rainbow(5),
        xlab="Individual #", ylab="Ancestry", border=NA,las=2)

## create admixture plots 
# for Q = 2
tbl_2<-read.table("proteus_p1_r66_minmaf05_locus_chr.2.Q")
fam <- read.table("proteus_p1_r66_minmaf05_locus_chr.fam")
rownames(tbl_2) <- fam$V2
colnames(tbl_2) <- c("V01", "V02")
# for Q = 3
tbl_3<-read.table("proteus_p1_r66_minmaf05_locus_chr.3.Q")
fam <- read.table("proteus_p1_r66_minmaf05_locus_chr.fam")
rownames(tbl_3) <- fam$V2
colnames(tbl_3) <- c("V01", "V02", "V03")
# for Q = 4
tbl_4<-read.table("proteus_p1_r66_minmaf05_locus_chr.4.Q")
fam <- read.table("proteus_p1_r66_minmaf05_locus_chr.fam")
rownames(tbl_4) <- fam$V2
colnames(tbl_4) <- c("V01", "V02", "V03", "V04")
# for Q = 5
tbl_5<-read.table("proteus_p1_r66_minmaf05_locus_chr.5.Q")
fam <- read.table("proteus_p1_r66_minmaf05_locus_chr.fam")
rownames(tbl_5) <- fam$V2
colnames(tbl_5) <- c("V01", "V02", "V03", "V04", "V05")
# for Q = 6
tbl_6<-read.table("proteus_p1_r66_minmaf05_locus_chr.6.Q")
fam <- read.table("proteus_p1_r66_minmaf05_locus_chr.fam")
rownames(tbl_6) <- fam$V2
colnames(tbl_6) <- c("V01", "V02", "V03", "V04", "V05", "V06")
# for Q = 7
tbl_7<-read.table("proteus_p1_r66_minmaf05_locus_chr.7.Q")
fam <- read.table("proteus_p1_r66_minmaf05_locus_chr.fam")
rownames(tbl_7) <- fam$V2
colnames(tbl_7) <- c("V01", "V02", "V03", "V04", "V05", "V06", "V07")
# for Q = 8
tbl_8<-read.table("proteus_p1_r66_minmaf05_locus_chr.8.Q")
fam <- read.table("proteus_p1_r66_minmaf05_locus_chr.fam")
rownames(tbl_8) <- fam$V2
colnames(tbl_8) <- c("V01", "V02", "V03", "V04", "V05", "V06", "V07", "V08")
# for Q = 9
tbl_9<-read.table("proteus_p1_r66_minmaf05_locus_chr.9.Q")
fam <- read.table("proteus_p1_r66_minmaf05_locus_chr.fam")
rownames(tbl_9) <- fam$V2
colnames(tbl_9) <- c("V01", "V02", "V03", "V04", "V05", "V06", "V07", "V08", "V09")
# for Q = 10
tbl_10<-read.table("proteus_p1_r66_minmaf05_locus_chr.10.Q")
fam <- read.table("proteus_p1_r66_minmaf05_locus_chr.fam")
rownames(tbl_10) <- fam$V2
colnames(tbl_10) <- c("V01", "V02", "V03", "V04", "V05", "V06", "V07", "V08", "V09","V10")
# using tidyverse
plot_data_2 <- tbl_2 %>% 
  mutate(id = rownames(tbl_2)) %>% 
  gather('pop', 'prob', V01:V02) %>% 
  group_by(id) %>% 
  mutate(likely_assignment = pop[which.max(prob)],
         assingment_prob = max(prob)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))
plot_data_2$id <- factor(plot_data_2$id,levels = c("PA331", "PC622", "PC621", "PC623", "PA329", "PA332", "PA334", "PC580", "PC581", "PC616", "PC617", "PC618", "PA040", "PA047", "PA049", "PA341", "PC587", "PC592", "PC614", "PC631-r", "PA354", "PC630", "PA255", "PA256", "PA343", "PA346", "PC582", "PA363", "PC578", "PC579", "PC586", "PC593"))

plot_data_3 <- tbl_3 %>% 
  mutate(id = rownames(tbl_3)) %>% 
  gather('pop', 'prob', V01:V03) %>% 
  group_by(id) %>% 
  mutate(likely_assignment = pop[which.max(prob)],
         assingment_prob = max(prob)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))
plot_data_3$id <- factor(plot_data_3$id,levels = c("PA331", "PC622", "PC621", "PC623", "PA329", "PA332", "PA334", "PC580", "PC581", "PC616", "PC617", "PC618", "PA040", "PA047", "PA049", "PA341", "PC587", "PC592", "PC614", "PC631-r", "PA354", "PC630", "PA255", "PA256", "PA343", "PA346", "PC582", "PA363", "PC578", "PC579", "PC586", "PC593"))

plot_data_4 <- tbl_4 %>% 
  mutate(id = rownames(tbl_4)) %>% 
  gather('pop', 'prob', V01:V04) %>% 
  group_by(id) %>% 
  mutate(likely_assignment = pop[which.max(prob)],
         assingment_prob = max(prob)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))
plot_data_4$id <- factor(plot_data_4$id,levels = c("PA331", "PC622", "PC621", "PC623", "PA329", "PA332", "PA334", "PC580", "PC581", "PC616", "PC617", "PC618", "PA040", "PA047", "PA049", "PA341", "PC587", "PC592", "PC614", "PC631-r", "PA354", "PC630", "PA255", "PA256", "PA343", "PA346", "PC582", "PA363", "PC578", "PC579", "PC586", "PC593"))

plot_data_5 <- tbl_5 %>% 
  mutate(id = rownames(tbl_5)) %>% 
  gather('pop', 'prob', V01:V05) %>% 
  group_by(id) %>% 
  mutate(likely_assignment = pop[which.max(prob)],
         assingment_prob = max(prob)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))
plot_data_5$id <- factor(plot_data_5$id,levels = c("PA331", "PC622", "PC621", "PC623", "PA329", "PA332", "PA334", "PC580", "PC581", "PC616", "PC617", "PC618", "PA040", "PA047", "PA049", "PA341", "PC587", "PC592", "PC614", "PC631-r", "PA354", "PC630", "PA255", "PA256", "PA343", "PA346", "PC582", "PA363", "PC578", "PC579", "PC586", "PC593"))

plot_data_6 <- tbl_6 %>% 
  mutate(id = rownames(tbl_6)) %>% 
  gather('pop', 'prob', V01:V06) %>% 
  group_by(id) %>% 
  mutate(likely_assignment = pop[which.max(prob)],
         assingment_prob = max(prob)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))
plot_data_6$id <- factor(plot_data_6$id,levels = c("PA331", "PC622", "PC621", "PC623", "PA329", "PA332", "PA334", "PC580", "PC581", "PC616", "PC617", "PC618", "PA040", "PA047", "PA049", "PA341", "PC587", "PC592", "PC614", "PC631-r", "PA354", "PC630", "PA255", "PA256", "PA343", "PA346", "PC582", "PA363", "PC578", "PC579", "PC586", "PC593"))

plot_data_7 <- tbl_7 %>% 
  mutate(id = rownames(tbl_7)) %>% 
  gather('pop', 'prob', V01:V07) %>% 
  group_by(id) %>% 
  mutate(likely_assignment = pop[which.max(prob)],
         assingment_prob = max(prob)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))
plot_data_7$id <- factor(plot_data_7$id,levels = c("PA331", "PC622", "PC621", "PC623", "PA329", "PA332", "PA334", "PC580", "PC581", "PC616", "PC617", "PC618", "PA040", "PA047", "PA049", "PA341", "PC587", "PC592", "PC614", "PC631-r", "PA354", "PC630", "PA255", "PA256", "PA343", "PA346", "PC582", "PA363", "PC578", "PC579", "PC586", "PC593"))

plot_data_8 <- tbl_8 %>% 
  mutate(id = rownames(tbl_8)) %>% 
  gather('pop', 'prob', V01:V08) %>% 
  group_by(id) %>% 
  mutate(likely_assignment = pop[which.max(prob)],
         assingment_prob = max(prob)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))
plot_data_8$id <- factor(plot_data_8$id,levels = c("PA331", "PC622", "PC621", "PC623", "PA329", "PA332", "PA334", "PC580", "PC581", "PC616", "PC617", "PC618", "PA040", "PA047", "PA049", "PA341", "PC587", "PC592", "PC614", "PC631-r", "PA354", "PC630", "PA255", "PA256", "PA343", "PA346", "PC582", "PA363", "PC578", "PC579", "PC586", "PC593"))

plot_data_9 <- tbl_9 %>% 
  mutate(id = rownames(tbl_9)) %>% 
  gather('pop', 'prob', V01:V09) %>% 
  group_by(id) %>% 
  mutate(likely_assignment = pop[which.max(prob)],
         assingment_prob = max(prob)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))
plot_data_9$id <- factor(plot_data_9$id,levels = c("PA331", "PC622", "PC621", "PC623", "PA329", "PA332", "PA334", "PC580", "PC581", "PC616", "PC617", "PC618", "PA040", "PA047", "PA049", "PA341", "PC587", "PC592", "PC614", "PC631-r", "PA354", "PC630", "PA255", "PA256", "PA343", "PA346", "PC582", "PA363", "PC578", "PC579", "PC586", "PC593"))

plot_data_10 <- tbl_10 %>% 
  mutate(id = rownames(tbl_10)) %>% 
  gather('pop', 'prob', V01:V10) %>% 
  group_by(id) %>% 
  mutate(likely_assignment = pop[which.max(prob)],
         assingment_prob = max(prob)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))
plot_data_10$id <- factor(plot_data_10$id,levels = c("PA331", "PC622", "PC621", "PC623", "PA329", "PA332", "PA334", "PC580", "PC581", "PC616", "PC617", "PC618", "PA040", "PA047", "PA049", "PA341", "PC587", "PC592", "PC614", "PC631-r", "PA354", "PC630", "PA255", "PA256", "PA343", "PA346", "PC582", "PA363", "PC578", "PC579", "PC586", "PC593"))

q2<- ggplot(plot_data_2, aes(id, prob, fill = pop)) +
  geom_col() +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_fill_brewer(palette='Paired') +
  theme_classic()
q3<- ggplot(plot_data_3, aes(id, prob, fill = pop)) +
  geom_col() +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_fill_brewer(palette='Paired') +
  theme_classic()
q4<- ggplot(plot_data_4, aes(id, prob, fill = pop)) +
  geom_col() +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_fill_brewer(palette='Paired') +
  theme_classic()
q5<- ggplot(plot_data_5, aes(id, prob, fill = pop)) +
  geom_col() +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_fill_brewer(palette='Paired') +
  theme_classic()
q6<- ggplot(plot_data_6, aes(id, prob, fill = pop)) +
  geom_col() +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_fill_brewer(palette='Paired') +
  theme_classic()
q7<- ggplot(plot_data_7, aes(id, prob, fill = pop)) +
  geom_col() +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_fill_brewer(palette='Paired') +
  theme_classic()
q8<- ggplot(plot_data_8, aes(id, prob, fill = pop)) +
  geom_col() +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_fill_brewer(palette='Paired') +
  theme_classic()
q9<- ggplot(plot_data_9, aes(id, prob, fill = pop)) +
  geom_col() +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_fill_brewer(palette='Paired') +
  theme_classic()
q10<- ggplot(plot_data_10, aes(id, prob, fill = pop)) +
  geom_col() +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_fill_brewer(palette='Paired') +
  theme_classic()
grid.arrange(q2,q3,q4,q5,q6,q7,q8,q9,q10,nrow=9)


#### FINERADSTRUCTURE ####
### 1) EDIT THE FOLLOWING THREE LINES TO PROVIDE PATHS TO THE fineRADstructure OUTPUT 
setwd("03_fineradstructure/") ## The directory where the files are located
chunkfile<-"proteus_haplotypes_filtered_reordered_chunks.out" ## RADpainter output file
mcmcfile<-"proteus_haplotypes_filtered_reordered_chunks.mcmc.xml" ## finestructure mcmc file
treefile<-"proteus_haplotypes_filtered_reordered_chunks.mcmcTree.xml" ## finestructure tree file
### 2) EDIT THIS PATH TO WHERE YOU WANT THE PLOTS:
plotsFolder <- "03_fineradstructure/"
### 3) SET VALUES FOR THESE VARIABLES: "analysisName" will be included in output plots
analysisName <- "Proteus population structure";  maxIndv <- 10000; maxPop<-10000
### 4) EDIT THE PATH TO YOUR COPY of FinestructureLibrary.R
source("C:/Users/hans_/Desktop/Projects/Asellus_cave_evolution/04_analysis/02_fineradstructure/FinestructureLibrary.R", chdir = TRUE) # read in the R functions, which also calls the needed packages
### 5) EXECUTE THE CODE ABOVE AND THE REST OF THE CODE BELOW
## make some colours
some.colors<-MakeColorYRP() # these are yellow-red-purple
some.colorsEnd<-MakeColorYRP(final=c(0.2,0.2,0.2)) # as above, but with a dark grey final for capped values
### READ IN THE CHUNKCOUNT FILE
dataraw<-as.matrix(read.table(chunkfile,row.names=1,header=T,skip=1)) # read in the pairwise coincidence 
### READ IN THE MCMC FILES
mcmcxml<-xmlTreeParse(mcmcfile) ## read into xml format
mcmcdata<-as.data.frame.myres(mcmcxml) ## convert this into a data frame
### READ IN THE TREE FILES
treexml<-xmlTreeParse(treefile) ## read the tree as xml format
ttree<-extractTree(treexml) ## extract the tree into ape's phylo format
## Reduce the amount of significant digits printed in the posteror assignment probabilities (numbers shown in the tree):
ttree$node.label[ttree$node.label!=""] <-format(as.numeric(ttree$node.label[ttree$node.label!=""]),digits=2)
# convert to dendrogram format
tdend<-myapetodend(ttree,factor=1)
## Now we work on the MAP state
mapstate<-extractValue(treexml,"Pop") # map state as a finestructure clustering
mapstatelist<-popAsList(mapstate) # .. and as a list of individuals in populations
popnames<-lapply(mapstatelist,NameSummary) # population names IN A REVERSIBLE FORMAT (I.E LOSSLESS)
## NOTE: if your population labels don't correspond to the format we used (NAME<number>) YOU MAY HAVE TROUBLE HERE. YOU MAY NEED TO RENAME THEM INTO THIS FORM AND DEFINE YOUR POPULATION NAMES IN popnamesplot BELOW
popnamesplot<-lapply(mapstatelist,NameMoreSummary) # a nicer summary of the populations
names(popnames)<-popnamesplot # for nicety only
names(popnamesplot)<-popnamesplot # for nicety only
popdend<-makemydend(tdend,mapstatelist) # use NameSummary to make popdend
popdend<-fixMidpointMembers(popdend) # needed for obscure dendrogram reasons
popdendclear<-makemydend(tdend,mapstatelist,"NameMoreSummary")# use NameMoreSummary to make popdend
popdendclear<-fixMidpointMembers(popdendclear) # needed for obscure dendrogram reasons
###
## Plot 1: COANCESTRY MATRIX
fullorder<-labels(tdend) # the order according to the tree
datamatrix<-dataraw[fullorder,fullorder] # reorder the data matrix
tmpmat<-datamatrix 
tmpmat[tmpmat>maxIndv]<-maxIndv #  # cap the heatmap
pdf(file=paste(plotsFolder,analysisName,"-SimpleCoancestry.pdf",sep=""),height=25,width=25)
plotFinestructure(tmpmat,dimnames(tmpmat)[[1]],dend=tdend,cols=some.colorsEnd,cex.axis=1.1,edgePar=list(p.lwd=0,t.srt=90,t.off=-0.1,t.cex=1.2))
dev.off()
###
## Plot 2: POPULATIONS AND COANCESTRY AVERAGES
popmeanmatrix<-getPopMeanMatrix(datamatrix,mapstatelist)
tmpmat<-popmeanmatrix
tmpmat[tmpmat>maxPop]<-maxPop # cap the heatmap
pdf(file=paste(plotsFolder,analysisName,"-PopAveragedCoancestry.pdf",sep=""),height=20,width=20)
plotFinestructure(tmpmat,dimnames(tmpmat)[[1]],dend=tdend,cols=some.colorsEnd,cex.axis=1.1,edgePar=list(p.lwd=0,t.srt=90,t.off=-0.1,t.cex=1.2))
dev.off()
###
## Plot 3: POPULATIONS AND COANCESTRY AVERAGES WITH PERHAPS MORE INFORMATIVE LABELS
mappopcorrectorder<-NameExpand(labels(popdend))
mappopsizes<-sapply(mappopcorrectorder,length)
labellocs<-PopCenters(mappopsizes)
xcrt=0
ycrt=45
pdf(file=paste(plotsFolder,analysisName,"-PopAveragedCoancestry2.pdf",sep=""),height=25,width=25)
plotFinestructure(tmpmat,dimnames(tmpmat)[[1]],labelsx=labels(popdendclear),labelsatx=labellocs,xcrt=xcrt,cols=some.colorsEnd,ycrt=ycrt,dend=tdend,cex.axis=1.1,edgePar=list(p.lwd=0,t.srt=90,t.off=-0.1,t.cex=1.2),hmmar=c(3,0,0,1))
dev.off()


#### PCA ####
## adegenet package
setwd("04_pca/")
all.proteus <- read.PLINK("p6_r66_single_snp_minmaf05_all.raw", parallel=FALSE)
istria.rem.proteus <- read.PLINK("p6_r66_single_snp_minmaf05_rem_istria.raw", parallel=FALSE)
slo.proteus <- read.PLINK("p4_r66_single_snp_minmaf05_slovenia.raw", parallel=FALSE)
slo.light.proteus <- read.PLINK("p2_r66_single_snp_minmaf05_slo_light.raw", parallel=FALSE)
cro.bos.proteus <- read.PLINK("p2_r66_single_snp_minmaf05_cro_bos.raw", parallel=FALSE)

pca1 <- glPca(all.proteus, nf = 11)
pca2 <- glPca(istria.rem.proteus, nf = 11)
pca3 <- glPca(slo.proteus, nf = 6)
pca4 <- glPca(slo.light.proteus, nf = 5)
pca4 <- glPca(cro.bos.proteus, nf = 5)
# simple plot to check
scatter(pca1, posi="bottomright")
add.scatter.eig(pca1$eig[1:10],2,3,4, posi="topright", inset=.05, ratio=.3)
s.label(pca1$scores,xax=3,yax=5)
add.scatter.eig(ca1$eig,nf=3,xax=2,yax=3,posi="topleft")
scatter(pca2, posi="bottomright")
add.scatter.eig(pca2$eig[1:10],2,1,2, posi="topright", inset=.05, ratio=.3)
scatter(pca3, posi="topleft")
add.scatter.eig(pca3$eig[1:10],2,1,2, posi="topright", inset=.05, ratio=.3)
scatter(pca4, posi="topright")
add.scatter.eig(pca4$eig[1:10],2,1,2, posi="topright", inset=.05, ratio=.3)
# colour plot
myCol <- colorplot(pca1$scores,pca1$scores, transp=TRUE, cex=4)
abline(h=0,v=0, col="grey")
add.scatter.eig(pca1$eig[1:10],2,1,2, posi="topright", inset=.08, ratio=.2)
myCol <- colorplot(pca1$scores,pca1$scores, axes=1,3,4, transp=TRUE, cex=4)
abline(h=0,v=0, col="grey")
add.scatter.eig(pca1$eig[1:10],2,3,4, posi="topright", inset=.08, ratio=.2)
myCol <- colorplot(pca2$scores,pca2$scores, transp=TRUE, cex=4)
abline(h=0,v=0, col="grey")
add.scatter.eig(pca$eig[1:10],2,1,2, posi="topleft", inset=.08, ratio=.2)
myCol <- colorplot(pca3$scores,pca3$scores, transp=TRUE, cex=4)
abline(h=0,v=0, col="grey")
add.scatter.eig(pca3$eig[1:10],2,1,2, posi="topleft", inset=.08, ratio=.2)
colorplot(pca4$scores,pca4$scores, transp=TRUE, cex=4)
abline(h=0,v=0, col="grey")
add.scatter.eig(pca4$eig[1:10],2,1,2, posi="topleft", inset=.08, ratio=.2)

### Fst correlation plots ###
setwd("05_fst/")
dat <- read.csv("Proteus_fst_vs_dxy_matrix.tsv", sep = "\t")
data <- read.csv("Proteus_fst_vs_dxy_matrix_2.tsv", header=FALSE, sep = "\t")
rows<-(dat$X)
cols<-(dat$X)
colnames(data)<-cols
rownames(data)<-rows
names(data) <- cols
fst<-as.matrix(data)
dxy<-as.matrix(data)
dxy[upper.tri(dxy)]=NA
fst[lower.tri(fst)]=NA
pheatmap(dxy, display_numbers = T, number_format='%.3f', color = colorRampPalette(c('goldenrod2','firebrick3'))(100), cluster_rows = F, cluster_cols = F, fontsize_number = 15)
pheatmap(fst, display_numbers = T, number_format="%.2f", color = colorRampPalette(c('goldenrod2','firebrick3'))(100), cluster_rows = F, cluster_cols = F, fontsize_number = 15)

#### TREEMIX ####
setwd("06_treemix/")
source("plotting_funcs.R")
plot_tree("proteus_treemix")
plot_tree("proteus_treemix_sz_corr")
plot_tree("proteus_treemix_m1")
plot_tree("proteus_treemix_m2")
plot_tree("proteus_treemix_m3")
plot_tree("proteus_treemix_m4")
plot_tree("proteus_treemix_m5")
plot_resid("proteus_treemix", "poporder.txt")
plot_resid("proteus_treemix_sz_corr", "poporder.txt")
plot_resid("proteus_treemix_m1", "poporder.txt")
plot_resid("proteus_treemix_m2", "poporder.txt")
plot_resid("proteus_treemix_m3", "poporder.txt")
treemix_res <- optM("C:/Users/hans_/Desktop/Projects/Proteus_phylogeny/04_analysis/05_treemix/02_results/", method = "Evanno", thresh = 0.05)
plot_optM(treemix_res, method = "Evanno", plot = TRUE, pdf = "evanno_treemix.pdf")

#### BPP ####
#create random sample of 100 loci from the 29145 loci
ran_samples_100 <- sample(x = 1:29145, size = 100)
# load mcmc file and convert to data matrix
mcmc <- read.table("proteus_100_loci_mcmc.txt",header=TRUE)
mcmc_proteus <- as.data.frame(mcmc)
# load tree
#tree <- readMCMCtree("FigTree.tre")
tree <- read.tree("species_tree.tre")

mcmc.summary(mcmc, prob = 0.95)
test1 <- curve(bppr::dslnorm(x, shift=14, meanlog=0, sdlog=0.01), from=14, to=20, n=1000)
str(test1)
max(test1$y)
root_calmsc <- msc2time.t(mcmc=mcmc_proteus, 
                     node="10Para_litKrajinaLikaKras_CarsoSW_SloSticnaSE_SloparkeljIstria", 
                     calf=rslnorm, shift=15.8, meanlog=0.01, sdlog=0.1)
mcmc.summary(root_calmsc)
root_RAD_BSC2 <- msc2time.t(mcmc=mcmc_proteus, 
                          node="10Para_litKrajinaLikaKras_CarsoSW_SloSticnaSE_SloparkeljIstria", 
                          calf=rslnorm, shift=15.4, meanlog=0.1, sdlog=0.1)
mcmc.summary(root_RAD_BSC2)
root_mtDNA_starbeast2 <- msc2time.t(mcmc=mcmc_proteus, 
                          node="10Para_litKrajinaLikaKras_CarsoSW_SloSticnaSE_SloparkeljIstria", 
                          calf=rslnorm, shift=14.7, meanlog=0.1, sdlog=0.1)
mcmc.summary(root_mtDNA_starbeast2)

mcmc2densitree(tree, root_calmsc,"t_", thin=0.05, alpha=0.01)

root_calmsc_uni <- msc2time.t(mcmc=mcmc_proteus, 
                          node="10Para_litKrajinaLikaKras_CarsoSW_SloSticnaSE_SloparkeljIstria", 
                          calf=runif, shift=8.4, max=20.2)

plot(density(root_calmsc$t_10Para_litKrajinaLikaKras_CarsoSW_SloSticnaSE_SloparkeljIstria, adj=.1), xlab="Time (Ma)",
     main="root age")
rug(root_calmsc$t_10Para_litKrajinaLikaKras_CarsoSW_SloSticnaSE_SloparkeljIstria)

mcmc2densitree(tree, root_calmsc,"t_", thin=0.05, alpha=0.01)
title(xlab="Divergence time (Ma)")


mcmc.summary(root_calmsc)
apply(root_calmsc, 2, mean)
coda::HPDinterval(coda::as.mcmc(root_calmsc))
curve(bppr::dslnorm(x, shift=3, meanlog=0.5, sdlog=0.5), from=0, to=15, n=1e3)
root_uplift_calmsc_ln <- msc2time.t(mcmc=mcmc_proteus, 
                     node="14LikaKras_Carso", 
                     calf=runif, shift=3,  meanlog=0.5, sdlog=0.5)
root_uplift_calmsc_nd <- msc2time.t(mcmc=mcmc_proteus, 
                                 node="14LikaKras_Carso", 
                                 calf=runif, min=3, max=5)
str(root_calmsc_uni)
mcmc.summary(root_uplift_calmsc_ln)
mcmc.summary(root_uplift_calmsc_nd)
apply(root_uplift_calmsc_nd, 2, mean)
mcmc2densitree(tree, root_uplift_calmsc,"t_", thin=0.05, alpha=0.01)


### MCMCtree ###

## shortest HPD sum tree
setwd("08_MCMCtree/01_root_LK_BSC2")
phy <- readMCMCtree("FigTree.tre")
mcmc <-read.table("mcmc.txt",header=TRUE)
MCMCtree.posterior <- as.data.frame(mcmc)
MCMC.tree.plot(phy, build.tree =TRUE, MCMC.chain = MCMCtree.posterior, cex.tips = 1.5, 
               time.correction = 100, plot.type = "distributions", cex.age = 1.5, 
               cex.labels = 1.5, col.tree = "grey40", 
               scale.res = c("Epoch"), add.time.scale	= TRUE, relative.height = 0.08, add.abs.time =TRUE,
               density.col = "#00000050", density.border.col = "#00000080")

## only Lika split and root
setwd("08_MCMCtree/01_LK")
phy <- readMCMCtree("FigTree.tre")
#phy <- read.tree("timetree_abs_ages.tre")
mcmc <-read.table("mcmc.txt",header=TRUE)
MCMCtree.posterior <- as.data.frame(mcmc)
node.post <-read.table("node_posteriors.txt",header=TRUE,check.names = FALSE)
node.posterior <- as.data.frame(node.post)
str(node.posterior)
MCMC.tree.plot(phy, build.tree =TRUE, MCMC.chain = MCMCtree.posterior, cex.tips = 1.5, 
               time.correction = 100, plot.type = "distributions", cex.age = 1.5, 
               cex.labels = 1.5, col.tree = "grey40", 
               scale.res = c("Epoch"), add.time.scale	= TRUE, relative.height = 0.08, add.abs.time =TRUE,
               density.col = "#00000050", density.border.col = "#00000080")

## mitochondrial DNA tree
setwd("08_MCMCtree/03_mtDNA")
phy <- readMCMCtree("FigTree.tre")
mcmc <-read.table("mcmc.txt",header=TRUE)
MCMCtree.posterior <- as.data.frame(mcmc)
MCMC.tree.plot(phy, build.tree =TRUE, MCMC.chain = MCMCtree.posterior, cex.tips = 1.5, 
               time.correction = 100, plot.type = "distributions", cex.age = 1.5, 
               cex.labels = 1.5, col.tree = "grey40", 
               scale.res = c("Epoch"), add.time.scale	= TRUE, relative.height = 0.08, add.abs.time =TRUE,
               density.col = "#00000050", density.border.col = "#00000080")


########### MtDNA analysis ########### 
### Plot different indices of diversity against samples and locations
setwd("08_arlequin")
moldiv <- read.csv("proteus_mol_diversity2.csv", header=TRUE)
str(moldiv)
plot(moldiv)
moldiv <- read.csv("area_vs_diversity.csv", header=TRUE)
str(moldiv)
plot(Pi~log.area., data=moldiv, pch=20, cex=3)
abline(lm(Pi~log.area., data=moldiv))
summary(lm(Pi~log.area., data=moldiv))
plot(Pi_genomics~log.area., data=moldiv, pch=20, cex=3)
abline(lm(Pi_genomics~log.area., data=moldiv))
summary(lm(Pi_genomics~log.area., data=moldiv))


### Phylogenetic analyses of mtDNA
### MAXIMUM PARSIMONY ###
## read in the trees
setwd("C:/Users/hans_/Desktop/Projects/Proteus_phylogeny/04_analysis/09_mtDNA/04_MP")
jack <- read.nexus("MP_tree_500_jackknife2")
boot<- read.nexus("MP_tree_500_bootstraps")
plot(jack)
plot(boot)
collapsed = collapse_tree_at_resolution(boot, resolution=380,shorten=FALSE,rename_collapsed_nodes = TRUE)$tree
plot(collapsed)
coll.tree<- collapse_tree_at_resolution(jack, 
                            resolution             = 0, 
                            by_edge_count          = FALSE,
                            shorten                = TRUE,
                            rename_collapsed_nodes = FALSE,
                            criterion              = 'max_tip_depth')
plot(coll.tree)


### ANCESTRAL TRAIT RECONSTRUCTION ###
## read in the trees
#corHMM
setwd("11_ancestral_state_reconstruction/")
morphs <- read.table("Proteus_states_ances.txt", header=TRUE)
str(morphs)
tree=read.tree("RAxML_bestTree.proteus_all_1-3snps_r66_MD_root_ancestor.tre")
# ancestral trait reconstruction
# null model
ARD_morph_fixed=rayDISC(tree,morphs,ntraits=1,charnum=1,model="ARD",node.states="marginal",root.p=c(1,0))
# equal rates
ARD_morph_equal=rayDISC(tree,morphs,ntraits=1,charnum=1,model="ER",node.states="marginal",root.p=c(1,0))
# no reversal
ARD_Dollo=rayDISC(tree,morphs,ntraits=1,charnum=1, p=c(0,1), model="ARD",node.states="marginal",root.p=c(1,0))
# no reversal equal rates
ARD_Dollo_equal=rayDISC(tree,morphs,ntraits=1,charnum=1, p=c(0,1), model="ER",node.states="marginal",root.p=c(1,0))
# no reversal equal rates 2
ARD_Dollo_equal=rayDISC(tree,morphs,ntraits=1,charnum=1, p=c(0,1), root.p=c(1,0))

plotRECON(tree,ARD_Dollo_equal$states,piecolors=col.list,
          cex=1.1, pie.cex=0.7)
ARD_Dollo_equal$tip.states

# null model rates
ASR_morph_rates=rayDISC(tree,morphs,ntraits=1,charnum=1, p=c(0.709, 1.892),model="ARD",node.states="marginal",root.p=c(1,0))

round(ARD_morph_fixed$states,2)
round(ARD_morph_fixed$states,2)

## plot the trees
col.list<-c("#52453A","#DDAF92")
plotRECON(tree,ARD_morph_equal$states,piecolors=col.list,
          cex=1.1, pie.cex=0.7,height=9, width=9, file="ER2.pdf")
plotRECON(tree,ARD_Dollo$states,piecolors=col.list,
          cex=1.1, pie.cex=0.7,height=9, width=9, file="Dollo-test2.pdf")
plotRECON(tree,ARD_morph_fixed$states,piecolors=col.list,
          cex=1.1, pie.cex=0.7)
plotRECON(tree,ARD_morph_equal$states,piecolors=col.list,
          cex=1.1, pie.cex=0.7)
plotRECON(tree,ASR_morph_rates$states,piecolors=col.list,
          cex=1.1, pie.cex=0.7)
plotRECON(tree,ARD_Dollo_equal$states,piecolors=col.list,
          cex=1.1, pie.cex=0.7)
plotRECON(tree,ARD_Dollo$states,piecolors=col.list,
          cex=1.1, pie.cex=0.7)


# calculate AICs
str(ARD_morph_fixed)
str(ARD_Dollo)
ARD_morph_fixed$AICc
ARD_morph_equal$AICc
ARD_Dollo$AICc
ARD_Dollo_equal$AICc
ASR_morph_rates$AICc
lr.test(ARD_morph_equal$loglik,ARD_Dollo$loglik)
ARD_morph_fixed$solution
ARD_Dollo$solution

### use only species as tips
setwd("11_ancestral_state_reconstruction/")
morphs <- read.table("Proteus_states_ances_species.txt", header=TRUE)
str(morphs)
tree=read.tree("species_tree.tre")
## start ancestral trait reconstruction
# null model
ARD_morph_fixed=rayDISC(tree,morphs,ntraits=1,charnum=1,model="ARD",node.states="marginal",root.p=c(1,0))
# equal rates
ARD_morph_equal=rayDISC(tree,morphs,ntraits=1,charnum=1,model="ER",node.states="marginal",root.p=c(1,0))
# no reversal
ARD_Dollo=rayDISC(tree,morphs,ntraits=1,charnum=1, p=c(0,1), model="ARD",node.states="marginal",root.p=c(1,0))
# no reversal equal rates
ARD_Dollo_equal=rayDISC(tree,morphs,ntraits=1,charnum=1, p=c(0,1), model="ER",node.states="marginal",root.p=c(1,0))
# no reversal equal rates 2
ARD_Dollo_equal=rayDISC(tree,morphs,ntraits=1,charnum=1, p=c(0,1), root.p=c(1,0))
plotRECON(tree,ARD_Dollo_equal$states,piecolors=col.list,
          cex=1.1, pie.cex=0.7)
ARD_Dollo_equal$tip.states
# null model rates
ASR_morph_rates=rayDISC(tree,morphs,ntraits=1,charnum=1, p=c(0.00001, 0.5),model="ARD",node.states="marginal",root.p=c(1,0))
round(ARD_morph_fixed$states,2)
round(ARD_morph_fixed$states,2)
plotRECON(tree,ASR_morph_rates$states,piecolors=col.list,
          cex=1.1, pie.cex=0.7)
ASR_morph_rates$AICc

## plot the trees
col.list<-c("#52453A","#DDAF92")
plotRECON(tree,ARD_morph_equal$states,piecolors=col.list,
          cex=1.1, pie.cex=0.7,height=9, width=9, file="species_tree_ER2.pdf")
plotRECON(tree,ARD_Dollo$states,piecolors=col.list,
          cex=1.1, pie.cex=0.7,height=9, width=9, file="species_tree_Dollo-test2.pdf")
plotRECON(tree,ARD_morph_fixed$states,piecolors=col.list,
          cex=1.1, pie.cex=0.7)
plotRECON(tree,ARD_morph_equal$states,piecolors=col.list,
          cex=1.1, pie.cex=0.7)
plotRECON(tree,ASR_morph_rates$states,piecolors=col.list,
          cex=1.1, pie.cex=0.7)
plotRECON(tree,ARD_Dollo_equal$states,piecolors=col.list,
          cex=1.1, pie.cex=0.7)
plotRECON(tree,ARD_Dollo$states,piecolors=col.list,
          cex=1.1, pie.cex=0.7)
# calculate AICs
str(ARD_morph_fixed)
str(ARD_Dollo)
ARD_morph_fixed$AICc
ARD_morph_equal$AICc
ARD_Dollo$AICc
ARD_Dollo_equal$AICc
ASR_morph_rates$AICc
lr.test(ARD_morph_equal$loglik,ARD_Dollo$loglik)
ARD_morph_fixed$solution
ARD_Dollo$solution



### SPECIES DELIMITATION PLOTTING ###

## read in the trees
library(rmutil)
setwd("12_species_delimitation")
twelveS.tree=read.tree("UPGMA_tree_12S.nwk")
dloop.tree=read.tree("UPGMA_tree_dloop.nwk")
CO1.tree=read.tree("UPGMA_tree_CO1.nwk")
#CO1 tree has a taxon which has been excluded later (PA134), -> exclude this tip
CO1.tree.new <- drop.tip(CO1.tree, c("PA134"), trim.internal = TRUE, subtree = FALSE,
         root.edge = 0, rooted = is.rooted(CO1.tree))

## read in species delimitations
species_12S <- read.csv("species_12S.csv", header = T,stringsAsFactors = TRUE)
species_12S_m <- species_12S[,-1]
rownames(species_12S_m) <- species_12S[,1]
species_12S_m
species_dloop <- read.csv("species_dloop.csv", header = T,stringsAsFactors = TRUE)
species_dloop_m <- species_dloop[,-1]
rownames(species_dloop_m) <- species_dloop[,1]
species_dloop_m
species_CO1 <- read.csv("species_CO1.csv", header = T,stringsAsFactors = TRUE)
species_CO1_m <- species_CO1[,-1]
rownames(species_CO1_m) <- species_CO1[,1]
species_CO1_m

# plot species delimitations on trees
trait.plot(twelveS.tree, species_12S_m, cols = list(ABGD = c("#E3191C", "#693D99","#2F9E2A","#F89896","darkolivegreen4","#FE9418","#E3D2FF","#1F77B3","#B2DE89"), 
                                                    ASAP = c("#E3191C", "#693D99","#2F9E2A","#1F77B3","#F89896","#B2DE89","#FE9418")))
trait.plot(dloop.tree, species_dloop_m, cols = list(ABGD = c("#E3191C", "#693D99","#2F9E2A","#F89896","#B2DE89","#FE9418","#E3D2FF","#1F77B3"), 
                                                    ASAP = c("#E3191C", "#693D99","#2F9E2A","#F89896","tan4","#B2DE89","#FDBE6F","#FE9418","#E3D2FF","#1F77B3")))
trait.plot(CO1.tree.new, species_CO1_m, cols = list(ABGD = c("#E3191C", "#693D99","#2F9E2A","#F89896","#B2DE89","#FE9418","#1F77B3"), 
                                                    ASAP = c("#E3191C", "#693D99","mediumorchid1","#2F9E2A","#F89896","darkolivegreen4","#B2DE89","seagreen1","#FDBE6F","#FE9418","brown2","yellow2","#E3D2FF","#1F77B3")))

### Add tree based analyses
# read in tree
beast.tree=read.tree("proteus_constant_pop.tre")
newtree <- drop.tip(beast.tree, c("PA766"), trim.internal = TRUE, subtree = FALSE,
                         root.edge = 0, rooted = is.rooted(beast.tree))
str(newtree)
newtree2 <- keep.tip(newtree, c("PA331","PC623","PA300","PA263","PC540", "PA816","PA260",
                                "PA315", "PC617","PA276","PA297","PA157","PA321","PA348","PA033","PA255",
                                "PB600","PA011","PA004"))
plot(newtree2)
## create input files separately for each gene
# GMYC & PTP results
results <- read.csv("results_GMYC_PTP_no_lin.csv", header=TRUE)
str(results)
results_m <- results[,-1]
rownames(results_m) <- results[,1]
results_m

#all.results <- merge(all.data, results, by = 'ID',all.y=TRUE)
#write.table(all.results,"all.results.txt",sep="\t")
# plot species delimitations on trees
trait.plot(newtree, results_m, cols = list(BPTP = c("#E3191C", "#693D99","mediumorchid1","#2F9E2A","#F89896","darkolivegreen4","darkolivegreen3","#B2DE89","seagreen1","#FDBE6F","#FE9418","brown2","yellow2","#E3D2FF","#1F77B3"), 
                                           PTP = c("#E3191C", "#693D99","mediumorchid1","#2F9E2A","#F89896","darkolivegreen4","darkolivegreen3","#B2DE89","seagreen1","#FDBE6F","#FE9418","brown2","yellow2","#E3D2FF","#1F77B3"),
                                           GMYC = c("#E3191C", "#693D99","mediumorchid1","#2F9E2A","#F89896","darkolivegreen4","darkolivegreen3","#B2DE89","seagreen1","#FDBE6F","#FE9418","brown2","yellow2","#E3D2FF","#1F77B3")))

### SPATIAL ANALYSIS ###
setwd("13_spatial_analysis")
## Read sequence data in fasta
gene.12S <- read.dna(file="proteus_12S_sorted_gaps_rem.fas", format="fasta")
gene.CO1 <- read.dna(file="proteus_CO1_sorted.fas", format="fasta")
gene.Dloop <- read.dna(file="proteus_DLoop_sorted_rem_gaps.fas", format="fasta")
class(gene.12S)
class(gene.CO1)
class(gene.Dloop)
# add population annotation and make genind object
proteus.12S.annot <- read.csv("proteus_12S_lineages_and_coordinates.csv", header=TRUE, row.names=1)
proteus.CO1.annot <- read.csv("proteus_CO1_lineages_and_coordinates.csv", header=TRUE, row.names=1)
proteus.Dloop.annot <- read.csv("proteus_Dloop_lineages_and_coordinates.csv", header=TRUE, row.names=1)
# combine data
proteus.12S.genind <- DNAbin2genind(x=gene.12S, pop=proteus.12S.annot[["Lineage"]])
proteus.CO1.genind <- DNAbin2genind(x=gene.CO1, pop=proteus.CO1.annot[["Lineage"]])
proteus.Dloop.genind <- DNAbin2genind(x=gene.Dloop, pop=proteus.Dloop.annot[["Lineage"]])
# add coordinates
proteus.12S.coord <- read.csv("proteus_12S_coordinates.csv", header=TRUE,quote="", dec=".", row.names=1)
proteus.CO1.coord <- read.csv("proteus_CO1_coordinates.csv", header=TRUE,quote="", dec=".", row.names=1)
proteus.Dloop.coord <- read.csv("proteus_Dloop_coordinates.csv", header=TRUE,quote="", dec=".", row.names=1)
proteus.Dloop.coord
proteus.CO1.coord
proteus.12S.coord
# Add coordinates - note identification of slots within object
proteus.12S.genind$other$xy <- proteus.12S.coord
proteus.CO1.genind$other$xy <- proteus.CO1.coord
proteus.Dloop.genind$other$xy <- proteus.Dloop.coord
## Mantel test across all lineages##
# Euclidean distance for individuals (plain ordinary distance matrix)
proteus.12S.dist <- dist(x=proteus.12S.genind, method="euclidean", diag=T, upper=T)
proteus.CO1.dist <- dist(x=proteus.CO1.genind, method="euclidean", diag=T, upper=T)
proteus.Dloop.dist <- dist(x=proteus.Dloop.genind, method="euclidean", diag=T, upper=T)
proteus.12S.dist
proteus.CO1.dist
proteus.Dloop.dist
?dist()
# Geographical distance
proteus.12S.gdist <- dist(x=proteus.12S.genind$other$xy, method="euclidean",diag=TRUE, upper=TRUE)
proteus.CO1.gdist <- dist(x=proteus.CO1.genind$other$xy, method="euclidean",diag=TRUE, upper=TRUE)
proteus.Dloop.gdist <- dist(x=proteus.Dloop.genind$other$xy, method="euclidean",diag=TRUE, upper=TRUE)
# Mantel test
proteus.12S.mantel <- mantel.randtest(m1=proteus.12S.dist, m2=proteus.12S.gdist,nrepet=1000)
proteus.CO1.mantel <- mantel.randtest(m1=proteus.CO1.dist, m2=proteus.CO1.gdist,nrepet=1000)
proteus.Dloop.mantel <- mantel.randtest(m1=proteus.Dloop.dist, m2=proteus.Dloop.gdist,nrepet=1000)
plot(proteus.12S.mantel, nclass=30)
plot(proteus.CO1.mantel, nclass=30)
plot(proteus.Dloop.mantel, nclass=30)

### perform analysis separately based on lineages
### CO1 gene
setwd("C:/Users/hans_/Desktop/Projects/Proteus_phylogeny/04_analysis/13_spatial_analysis")
## Read sequence data in fasta
Dolenjska.CO1 <- read.dna(file="Dolenjska_CO1.fas", format="fasta")
Ljubljanica.CO1 <- read.dna(file="Ljubljanica_CO1.fas", format="fasta")
Kras.CO1 <- read.dna(file="Kras_CO1.fas", format="fasta")
Paralittoral.CO1 <- read.dna(file="Paralittoral_CO1.fas", format="fasta")
Lika.CO1 <- read.dna(file="Lika_CO1.fas", format="fasta")
parkelj.CO1 <- read.dna(file="parkelj_CO1.fas", format="fasta")
Sticna.CO1 <- read.dna(file="Sticna_CO1.fas", format="fasta")
Krajina.CO1 <- read.dna(file="Krajina_CO1.fas", format="fasta")
Istria.CO1 <- read.dna(file="Istria_CO1.fas", format="fasta")
class(Dolenjska.CO1)
# combine data
Dolenjska.CO1.genind <- DNAbin2genind(x=Dolenjska.CO1)
Ljubljanica.CO1.genind <- DNAbin2genind(x=Ljubljanica.CO1)
Kras.CO1.genind <- DNAbin2genind(x=Kras.CO1)
Paralittoral.CO1.genind <- DNAbin2genind(x=Paralittoral.CO1)
Lika.CO1.genind <- DNAbin2genind(x=Lika.CO1)
parkelj.CO1.genind <- DNAbin2genind(x=parkelj.CO1)
Krajina.CO1.genind <- DNAbin2genind(x=Krajina.CO1)
Sticna.CO1.genind <- DNAbin2genind(x=Sticna.CO1)
Istria.CO1.genind <- DNAbin2genind(x=Istria.CO1)
# add coordinates
Dolenjska.CO1.coord <- read.csv("Dolenjska_CO1_coordinates.csv", header=TRUE,quote="", dec=".", row.names=1)
Ljubljanica.CO1.coord <- read.csv("Ljubljanica_CO1_coordinates.csv", header=TRUE,quote="", dec=".", row.names=1)
Kras.CO1.coord <- read.csv("Kras_CO1_coordinates.csv", header=TRUE,quote="", dec=".", row.names=1)
Paralittoral.CO1.coord <- read.csv("Paralittoral_CO1_coordinates.csv", header=TRUE,quote="", dec=".", row.names=1)
Lika.CO1.coord <- read.csv("Lika_CO1_coordinates.csv", header=TRUE,quote="", dec=".", row.names=1)
parkelj.CO1.coord <- read.csv("parkelj_CO1_coordinates.csv", header=TRUE,quote="", dec=".", row.names=1)
Sticna.CO1.coord <- read.csv("Sticna_CO1_coordinates.csv", header=TRUE,quote="", dec=".", row.names=1)
Krajina.CO1.coord <- read.csv("Krajina_CO1_coordinates.csv", header=TRUE,quote="", dec=".", row.names=1)
Istria.CO1.coord <- read.csv("Istria_CO1_coordinates.csv", header=TRUE,quote="", dec=".", row.names=1)
# Add coordinates - note identification of slots within object
Dolenjska.CO1.genind$other$xy <- Dolenjska.CO1.coord
Ljubljanica.CO1.genind$other$xy <- Ljubljanica.CO1.coord
Kras.CO1.genind$other$xy <- Kras.CO1.coord
Paralittoral.CO1.genind$other$xy <- Paralittoral.CO1.coord
Lika.CO1.genind$other$xy <- Lika.CO1.coord
parkelj.CO1.genind$other$xy <- parkelj.CO1.coord
Krajina.CO1.genind$other$xy <- Krajina.CO1.coord
Sticna.CO1.genind$other$xy <- Sticna.CO1.coord
Istria.CO1.genind$other$xy <- Istria.CO1.coord
# Geographical distance
Dolenjska.CO1.gdist <- dist(x=Dolenjska.CO1.genind$other$xy, method="euclidean",diag=TRUE, upper=TRUE)
Ljubljanica.CO1.gdist <- dist(x=Ljubljanica.CO1.genind$other$xy, method="euclidean",diag=TRUE, upper=TRUE)
Kras.CO1.gdist <- dist(x=Kras.CO1.genind$other$xy, method="euclidean",diag=TRUE, upper=TRUE)
Paralittoral.CO1.gdist <- dist(x=Paralittoral.CO1.genind$other$xy, method="euclidean",diag=TRUE, upper=TRUE)
Lika.CO1.gdist <- dist(x=Lika.CO1.genind$other$xy, method="euclidean",diag=TRUE, upper=TRUE)
parkelj.CO1.gdist <- dist(x=parkelj.CO1.genind$other$xy, method="euclidean",diag=TRUE, upper=TRUE)
Krajina.CO1.gdist <- dist(x=Krajina.CO1.genind$other$xy, method="euclidean",diag=TRUE, upper=TRUE)
Sticna.CO1.gdist <- dist(x=Sticna.CO1.genind$other$xy, method="euclidean",diag=TRUE, upper=TRUE)
Istria.CO1.gdist <- dist(x=Istria.CO1.genind$other$xy, method="euclidean",diag=TRUE, upper=TRUE)
# Euclidean distance for individuals (plain ordinary distance matrix)
Dolenjska.CO1.dist <- dist(x=Dolenjska.CO1.genind, method="euclidean", diag=T, upper=T)
Ljubljanica.CO1.dist <- dist(x=Ljubljanica.CO1.genind, method="euclidean", diag=T, upper=T)
Kras.CO1.dist <- dist(x=Kras.CO1.genind, method="euclidean", diag=T, upper=T)
Paralittoral.CO1.dist <- dist(x=Paralittoral.CO1.genind, method="euclidean", diag=T, upper=T)
Lika.CO1.dist <- dist(x=Lika.CO1.genind, method="euclidean", diag=T, upper=T)
parkelj.CO1.dist <- dist(x=parkelj.CO1.genind, method="euclidean", diag=T, upper=T)
Krajina.CO1.dist <- dist(x=Krajina.CO1.genind, method="euclidean", diag=T, upper=T)
Sticna.CO1.dist <- dist(x=Sticna.CO1.genind, method="euclidean", diag=T, upper=T)
Istria.CO1.dist <- dist(x=Istria.CO1.genind, method="euclidean", diag=T, upper=T)
# Mantel test
Dolenjska.CO1.mantel <- mantel.randtest(m1=Dolenjska.CO1.dist, m2=Dolenjska.CO1.gdist,nrepet=1000)
Ljubljanica.CO1.mantel <- mantel.randtest(m1=Ljubljanica.CO1.dist, m2=Ljubljanica.CO1.gdist,nrepet=1000)
Kras.CO1.mantel <- mantel.randtest(m1=Kras.CO1.dist, m2=Kras.CO1.gdist,nrepet=1000)
Paralittoral.CO1.mantel <- mantel.randtest(m1=Paralittoral.CO1.dist, m2=Paralittoral.CO1.gdist,nrepet=1000)
Lika.CO1.mantel <- mantel.randtest(m1=Lika.CO1.dist, m2=Lika.CO1.gdist,nrepet=1000)
parkelj.CO1.mantel <- mantel.randtest(m1=parkelj.CO1.dist, m2=parkelj.CO1.gdist,nrepet=1000)
Krajina.CO1.mantel <- mantel.randtest(m1=Krajina.CO1.dist, m2=Krajina.CO1.gdist,nrepet=1000)
Sticna.CO1.mantel <- mantel.randtest(m1=Sticna.CO1.dist, m2=Sticna.CO1.gdist,nrepet=1000)
Istria.CO1.mantel <- mantel.randtest(m1=Istria.CO1.dist, m2=Istria.CO1.gdist,nrepet=1000)
plot(Dolenjska.CO1.mantel, nclass=30)
plot(Ljubljanica.CO1.mantel, nclass=30)
plot(Kras.CO1.mantel, nclass=30)
plot(Paralittoral.CO1.mantel, nclass=30)
plot(Lika.CO1.mantel, nclass=30)
plot(parkelj.CO1.mantel, nclass=30)
plot(Krajina.CO1.mantel, nclass=30)
plot(Sticna.CO1.mantel, nclass=30)
plot(Istria.CO1.mantel, nclass=30)

### Dloop gene
setwd("C:/Users/hans_/Desktop/Projects/Proteus_phylogeny/04_analysis/13_spatial_analysis")
## Read sequence data in fasta
Dolenjska.Dloop <- read.dna(file="Dolenjska_Dloop.fas", format="fasta")
Ljubljanica.Dloop <- read.dna(file="Ljubljanica_Dloop.fas", format="fasta")
Kras.Dloop <- read.dna(file="Kras_Dloop.fas", format="fasta")
Paralittoral.Dloop <- read.dna(file="Paralittoral_Dloop.fas", format="fasta")
Lika.Dloop <- read.dna(file="Lika_Dloop.fas", format="fasta")
parkelj.Dloop <- read.dna(file="parkelj_Dloop.fas", format="fasta")
Sticna.Dloop <- read.dna(file="Sticna_Dloop.fas", format="fasta")
Krajina.Dloop <- read.dna(file="Krajina_Dloop.fas", format="fasta")
Istria.Dloop <- read.dna(file="Istria_Dloop.fas", format="fasta")
class(Dolenjska.Dloop)
# combine data
Dolenjska.Dloop.genind <- DNAbin2genind(x=Dolenjska.Dloop)
Ljubljanica.Dloop.genind <- DNAbin2genind(x=Ljubljanica.Dloop)
Kras.Dloop.genind <- DNAbin2genind(x=Kras.Dloop)
Paralittoral.Dloop.genind <- DNAbin2genind(x=Paralittoral.Dloop)
Lika.Dloop.genind <- DNAbin2genind(x=Lika.Dloop)
parkelj.Dloop.genind <- DNAbin2genind(x=parkelj.Dloop)
Krajina.Dloop.genind <- DNAbin2genind(x=Krajina.Dloop)
Sticna.Dloop.genind <- DNAbin2genind(x=Sticna.Dloop)
Istria.Dloop.genind <- DNAbin2genind(x=Istria.Dloop)
# add coordinates
Dolenjska.Dloop.coord <- read.csv("Dolenjska_Dloop_coordinates.csv", header=TRUE,quote="", dec=".", row.names=1)
Ljubljanica.Dloop.coord <- read.csv("Ljubljanica_Dloop_coordinates.csv", header=TRUE,quote="", dec=".", row.names=1)
Kras.Dloop.coord <- read.csv("Kras_Dloop_coordinates.csv", header=TRUE,quote="", dec=".", row.names=1)
Paralittoral.Dloop.coord <- read.csv("Paralittoral_Dloop_coordinates.csv", header=TRUE,quote="", dec=".", row.names=1)
Lika.Dloop.coord <- read.csv("Lika_Dloop_coordinates.csv", header=TRUE,quote="", dec=".", row.names=1)
parkelj.Dloop.coord <- read.csv("parkelj_Dloop_coordinates.csv", header=TRUE,quote="", dec=".", row.names=1)
Sticna.Dloop.coord <- read.csv("Sticna_Dloop_coordinates.csv", header=TRUE,quote="", dec=".", row.names=1)
Krajina.Dloop.coord <- read.csv("Krajina_Dloop_coordinates.csv", header=TRUE,quote="", dec=".", row.names=1)
Istria.Dloop.coord <- read.csv("Istria_Dloop_coordinates.csv", header=TRUE,quote="", dec=".", row.names=1)
# Add coordinates - note identification of slots within object
Dolenjska.Dloop.genind$other$xy <- Dolenjska.Dloop.coord
Ljubljanica.Dloop.genind$other$xy <- Ljubljanica.Dloop.coord
Kras.Dloop.genind$other$xy <- Kras.Dloop.coord
Paralittoral.Dloop.genind$other$xy <- Paralittoral.Dloop.coord
Lika.Dloop.genind$other$xy <- Lika.Dloop.coord
parkelj.Dloop.genind$other$xy <- parkelj.Dloop.coord
Krajina.Dloop.genind$other$xy <- Krajina.Dloop.coord
Sticna.Dloop.genind$other$xy <- Sticna.Dloop.coord
Istria.Dloop.genind$other$xy <- Istria.Dloop.coord
# Geographical distance
Dolenjska.Dloop.gdist <- dist(x=Dolenjska.Dloop.genind$other$xy, method="euclidean",diag=TRUE, upper=TRUE)
Ljubljanica.Dloop.gdist <- dist(x=Ljubljanica.Dloop.genind$other$xy, method="euclidean",diag=TRUE, upper=TRUE)
Kras.Dloop.gdist <- dist(x=Kras.Dloop.genind$other$xy, method="euclidean",diag=TRUE, upper=TRUE)
Paralittoral.Dloop.gdist <- dist(x=Paralittoral.Dloop.genind$other$xy, method="euclidean",diag=TRUE, upper=TRUE)
Lika.Dloop.gdist <- dist(x=Lika.Dloop.genind$other$xy, method="euclidean",diag=TRUE, upper=TRUE)
parkelj.Dloop.gdist <- dist(x=parkelj.Dloop.genind$other$xy, method="euclidean",diag=TRUE, upper=TRUE)
Krajina.Dloop.gdist <- dist(x=Krajina.Dloop.genind$other$xy, method="euclidean",diag=TRUE, upper=TRUE)
Sticna.Dloop.gdist <- dist(x=Sticna.Dloop.genind$other$xy, method="euclidean",diag=TRUE, upper=TRUE)
Istria.Dloop.gdist <- dist(x=Istria.Dloop.genind$other$xy, method="euclidean",diag=TRUE, upper=TRUE)
# Euclidean distance for individuals (plain ordinary distance matrix)
Dolenjska.Dloop.dist <- dist(x=Dolenjska.Dloop.genind, method="euclidean", diag=T, upper=T)
Ljubljanica.Dloop.dist <- dist(x=Ljubljanica.Dloop.genind, method="euclidean", diag=T, upper=T)
Kras.Dloop.dist <- dist(x=Kras.Dloop.genind, method="euclidean", diag=T, upper=T)
Paralittoral.Dloop.dist <- dist(x=Paralittoral.Dloop.genind, method="euclidean", diag=T, upper=T)
Lika.Dloop.dist <- dist(x=Lika.Dloop.genind, method="euclidean", diag=T, upper=T)
parkelj.Dloop.dist <- dist(x=parkelj.Dloop.genind, method="euclidean", diag=T, upper=T)
Krajina.Dloop.dist <- dist(x=Krajina.Dloop.genind, method="euclidean", diag=T, upper=T)
Sticna.Dloop.dist <- dist(x=Sticna.Dloop.genind, method="euclidean", diag=T, upper=T)
Istria.Dloop.dist <- dist(x=Istria.Dloop.genind, method="euclidean", diag=T, upper=T)
# Mantel test
Dolenjska.Dloop.mantel <- mantel.randtest(m1=Dolenjska.Dloop.dist, m2=Dolenjska.Dloop.gdist,nrepet=1000)
Ljubljanica.Dloop.mantel <- mantel.randtest(m1=Ljubljanica.Dloop.dist, m2=Ljubljanica.Dloop.gdist,nrepet=1000)
Kras.Dloop.mantel <- mantel.randtest(m1=Kras.Dloop.dist, m2=Kras.Dloop.gdist,nrepet=1000)
Paralittoral.Dloop.mantel <- mantel.randtest(m1=Paralittoral.Dloop.dist, m2=Paralittoral.Dloop.gdist,nrepet=1000)
Lika.Dloop.mantel <- mantel.randtest(m1=Lika.Dloop.dist, m2=Lika.Dloop.gdist,nrepet=1000)
parkelj.Dloop.mantel <- mantel.randtest(m1=parkelj.Dloop.dist, m2=parkelj.Dloop.gdist,nrepet=1000)
Krajina.Dloop.mantel <- mantel.randtest(m1=Krajina.Dloop.dist, m2=Krajina.Dloop.gdist,nrepet=1000)
Sticna.Dloop.mantel <- mantel.randtest(m1=Sticna.Dloop.dist, m2=Sticna.Dloop.gdist,nrepet=1000)
Istria.Dloop.mantel <- mantel.randtest(m1=Istria.Dloop.dist, m2=Istria.Dloop.gdist,nrepet=1000)
plot(Dolenjska.Dloop.mantel, nclass=30)
plot(Ljubljanica.Dloop.mantel, nclass=30)
plot(Kras.Dloop.mantel, nclass=30)
plot(Paralittoral.Dloop.mantel, nclass=30)
plot(Lika.Dloop.mantel, nclass=30)
plot(parkelj.Dloop.mantel, nclass=30)
plot(Krajina.Dloop.mantel, nclass=30)
plot(Sticna.Dloop.mantel, nclass=30)
plot(Istria.Dloop.mantel, nclass=30)
# test if dimensions are the same in case of errors
#dim(as.matrix(Ljubljanica.Dloop.dist))
#dim(as.matrix(Ljubljanica.Dloop.gdist))
## All Mantel tests with sufficient sample sizes:
Dolenjska.CO1.mantel
Ljubljanica.CO1.mantel
Kras.CO1.mantel
Paralittoral.CO1.mantel
Lika.CO1.mantel
Krajina.CO1.mantel
Sticna.CO1.mantel 

Dolenjska.Dloop.mantel
Ljubljanica.Dloop.mantel
Kras.Dloop.mantel
Paralittoral.Dloop.mantel
Lika.Dloop.mantel
Istria.Dloop.mantel 
>>>>>>> 0eeaa07fdd09f5b40a347605c2b3755f580cb76a
