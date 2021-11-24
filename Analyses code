#Script to repeat analyses for the article: Gonzalez-Voyer et al. In press. Sex roles
#in birds: phylogenetic analyses of the influence of climate, life histories and social
#environment. Ecology Letters.
#authors: Alejandro Gonzalez Voyer & Gavin H. Thomas


library(ape)
library(phytools)
library(caper)
library(geiger)
library(phylolm)



#Determine the relationship among different estimates of the four major components of sex roles, i. e. mate competition role (traits related with sexual size dimorphism), mate courtship role (estimates of plumage dimorphism), pair bonding role (mating system and related traits) and, parental care role (male-female difference in care provisioning). 
#Some curation of the database that is publically available in Dryad may be necessary to run this script.

#This is the phylogeny that was used in all analyses

phy<-read.tree("BirdzillaNivosusNew.tre")

#Code to add Charadrius nivosus to the Jetz et al. phylogeny:
#which(phyT$tip.label=="Charadrius_alexandrinus")
#phy<-bind.tip(phyT, "Charadrius_nivosus", where=2388, position=0.5*phyT$edge.length[which(phyT$edge[,2]==2388)])


#Filter the dataset to make sure we have the same number of species as those in the tree, as the dataset has 10377 species whereas the tree has 9994 (having added C nivosus, see Methods in the accompanying publication).
#There are minor differences in spelling for two species in the dataset compared to the tree, therefore Brachypteracias_squamigera will change to Brachypteracias_squamiger and Nectarinia_neergaardi will change to Nectarinia_neergardi.  

dat$species<-as.character(dat$species)
dat$species[dat$species=="Brachypteracias_squamigera"]<-"Brachypteracias_squamiger"
dat$species[dat$species=="Nectarinia_neergaardi"]<-"Nectarinia_neergardi"
row.names(dat)<-dat$species
drop.data<-setdiff(dat$species, phy$tip.label)
exclude <- drop.data
dat<-dat[is.na(match(row.names(dat), exclude)), ] 
setdiff(phy$tip.label, dat$species)
setdiff(dat$species,phy$tip.label)


## 1.Mate competition role

#The mate competition role is mainly related to sexual size dimorphism. Compare dimorphism in body mass with dimorphism in tail length, tarsus length and wing length. 

mass.ssd<-log10(dat$m.mass/dat$f.mass)
wing.ssd<-log10(dat$m.wing/dat$f.wing)
tarsus.ssd<-log10(dat$m.tarsus/dat$f.tarsus)
species<-as.character(dat$species)
ssd<-as.data.frame(cbind(mass.ssd, wing.ssd, tarsus.ssd))
ssd$species<-species
ssd<-ssd[,c(4,1,2,3)]


#Look at distribution of the different estimates of SSD: 
  
hist(ssd$mass.ssd, main="Histogram of Mass SSD", xlab="Mass SSD")

hist(ssd$wing.ssd, main="Histogram of Wing SSD", xlab="Wing SSD")

hist(ssd$tarsus.ssd, main="Histogram of Tarsus SSD", xlab="Tarsus SSD")


#The correlation between them (to speed up the analyses without controlling for phylogeny):

cor1<-cor.test(ssd$mass.ssd, ssd$wing.ssd)

cor1

plot(ssd$wing.ssd~ssd$mass.ssd, xlab="Mass SSD", ylab="Wing SSD")

cor2<-cor.test(ssd$mass.ssd, ssd$tarsus.ssd)

cor2

plot(ssd$tarsus.ssd~ssd$mass.ssd, xlab="Mass SSD", ylab="Tarsus SSD")

cor3<-cor.test(ssd$wing.ssd, ssd$tarsus.ssd)

cor3

plot(ssd$wing.ssd~ssd$tarsus.ssd, xlab="Tarsus SSD", ylab="Wing SSD")


#Given the high correlations between the three estimates of SSD, we can use them in combination to estimate SSD.   

ssdMean<-rowMeans(ssd[,2:4], na.rm=TRUE)
ssdMean<-as.data.frame(cbind(ssd$species, ssdMean),stringsAsFactors=FALSE)
colnames(ssdMean)<-c("species", "ssdM")
ssdMean$ssdM<-as.numeric(as.character(ssdMean$ssdM))
hist(na.omit(ssdMean$ssdM), xlab="Mean SSD", main="Histogram of Averaged SSD")


## 2.Mate courting role

#Mean of the different estimates of plumage dimorphism, namely that of the back, belly, head, tail and wings.   

ornament<-cbind(dat$pl.head, dat$pl.back, dat$pl.belly, dat$pl.wings, dat$pl.tail)
dat$pl.mean2<-rowMeans(ornament, na.rm=TRUE)

#Histogram of distribution of sexual dichromatism values:

hist(dat$pl.mean2, xlab="Plumage dimorphism", main="")


## 3. Pair bonding role

#For this sex role component the most straight forward proxie is mating system, based on male polygamy - female polygamy. 

mating.sys<-dat$mpg.scr-dat$fpg.scr
#Sample size for mating system:
mate.sys<-mating.sys[!is.na(mating.sys)]
length(mate.sys)

hist(mating.sys, xlab="Mating system")

## 4.Parental role

#There are 8 variables scored for the male investment in each of them. However, there are differences in the sample size among the variables. One way to maximize sample size could be to pick those variables for which we have data for most species, if these correlate well with a mean calculated from all variables. 
#Degree of coverage for each of the different care categories:   

nest<-length(na.omit(dat$nest.bld)) 
inc<-length(na.omit(dat$inc)) 
nestgr<-length(na.omit(dat$nest.grd)) 
brood<-length(na.omit(dat$chick.brd)) 
feed<-length(na.omit(dat$chick.feed))
def<-length(na.omit(dat$chick.dfc)) 
postf<-length(na.omit(dat$postf.feed))
postfg<-length(na.omit(dat$postf.grd)) 


#Differences in the number of species for which data is available for the different care categories:   

setNames(c(nest,inc, nestgr, brood, feed, def, postf, postfg), c("NestBuilding","Incubation","NestGuarding","ChickBrooding","ChickFeeding","ChickDefence","PostFledgFeeding","PostFledgeGuarding"))

#Correlation between the mean for the 4 variables for which there is data for more than 1 000 species with the mean obtained for all care categories, for species for which we have data on all variables.   

care<-with(dat, data.frame(species,nest.bld,inc,chick.brd,chick.feed,chick.dfc,nest.grd,postf.feed,postf.grd))
careALL<-na.omit(care)
care4<-with(careALL, data.frame(nest.bld, inc, chick.brd, chick.feed))
mean4<-rowMeans(care4)
meanA<-rowMeans(careALL[,2:9])

#Given the reasonably high correlation between the two means we could use the mean of all variables to represent male parental investment.   

corel<-cor.test(mean4, meanA)
corel

#To center the care data on 0, which should reflect equal invesment by males and females, substract the value of 2 to all cases.

meanCare<-rowMeans(care[,2:9], na.rm=TRUE)
meanCare<-meanCare-2
mean(meanCare,na.rm=TRUE)

#To ensure this variable follows the order of all others, where positive values indicate conventional sex roles and negative values reversed sex roles, invert the sign of the values:

meanCareN<--meanCare
mean(meanCareN,na.rm=TRUE)
hist(meanCareN)


#Create a dataframe with the 4 components of sex roles:
  
Sexrole<-as.data.frame(cbind(ssdMean, dat$pl.mean, mating.sys, meanCareN), stringsAsFactors=FALSE)
colnames(Sexrole)<-c("species", "ssd", "dichro", "mating.sys", "care")

###Do sex role component means differ from zero?

#Run a PGLS analysis with each sex role component as the response an a vector of 1 as the predictor and use the intercept.
SSD<-Sexrole[,1:2]
SSD<-na.omit(SSD)
row.names(SSD)<-SSD$species
to.dropSSD<-setdiff(phy$tip.label,as.character(SSD$species))
phy.SSD<-drop.tip(phy, to.dropSSD)
ssd1<-phylolm(ssd~1, data=SSD, phy=phy.SSD, model="lambda")
summary(ssd1)

#Now for plumage dichromatism:

Dichro<-Sexrole[,c(1,3)]
Dichro<-na.omit(Dichro)
row.names(Dichro)<-Dichro$species
to.dropDichro<-setdiff(phy$tip.label,as.character(Dichro$species))
phy.Dichro<-drop.tip(phy, to.dropDichro)
dichro1<-phylolm(dichro~1, data=Dichro, phy=phy.Dichro, model="lambda")
summary(dichro1)

#Now on to mating system:

Matesys<-Sexrole[,c(1,4)]
Matesys<-na.omit(Matesys)
row.names(Matesys)<-Matesys$species
to.dropMatesys<-setdiff(phy$tip.label,as.character(Matesys$species))
phy.Matesys<-drop.tip(phy, to.dropMatesys)
mate1<-phylolm(mating.sys~1, data=Matesys, phy=phy.Matesys, model="lambda")
summary(mate1)

#Finally parental care:
  
PCare<-Sexrole[,c(1,5)]
PCare<-na.omit(PCare)
row.names(PCare)<-PCare$species
to.dropPCare<-setdiff(phy$tip.label,as.character(PCare$species))
phy.PCare<-drop.tip(phy, to.dropPCare)
care1<-phylolm(care~1, data=PCare, phy=phy.PCare, model="lambda")
summary(care1)

#Based on the results above, all the mean values of the Sex role components are not significantly different from 0. 


#Relationships among different components of sex roles

#For the phylogenetic PCA we need to standardize the variance of traits, so that traits with more variance don't have undue influence in the PCA.

SSDMean<-ssdMean$ssdM/sd(na.omit(ssdMean$ssdM))
pl.mean2<-dat$pl.mean/sd(na.omit(dat$pl.mean))
Mating.sys<-mating.sys/sd(na.omit(mating.sys))
MeanCare1<-meanCareN/sd(na.omit(meanCareN))
SexroleSD<-as.data.frame(cbind(SSDMean, pl.mean2, Mating.sys,MeanCare1))
SexroleSD<-cbind(dat$species, SexroleSD)
colnames(SexroleSD)<-c("species","mean.ssd","pl.mean","mating.sys","care")

#Sample size for which we have data on all 4 traits

summary(SexroleSD)
nrow(na.omit(SexroleSD))

PCAdata<-na.omit(Sexrole)
PCAdata<-droplevels(PCAdata)
row.names(PCAdata)<-PCAdata$species
PCAdata<-PCAdata[,-1]
#Cut tree to match dataframe:
to.drop<-setdiff(phy$tip.label,row.names(PCAdata))
phycut<-drop.tip(phy, to.drop)
setdiff(phycut$tip.label,row.names(PCAdata))



#Phylogenetic PCA, with a  lambda model of trait evolution 

phylo.pca<-phyl.pca(phycut,PCAdata, method="lambda", mode="corr")   

phylo.pca

biplot(phylo.pca, xlabs=rep("o",1861))


#Histogram of variance explained for the phylo PCA:   

Eval<-diag(phylo.pca$Eval)
p.var<-Eval/sum(Eval)
barplot(100*p.var, las=2, xlab='', ylab="% Variance Explained")

#To estimate the sex roles distance variable:

pc_means <- colMeans(phylo.pca$S)

sr_dist <- matrix(NA, ncol=1, nrow=dim(phylo.pca$S)[1])
rownames(sr_dist) <- rownames(phylo.pca$S)

for (i in 1:dim(sr_dist)[1]) {
  sr_dist[i,] <- dist(rbind(phylo.pca$S[i,], pc_means))
}



proles <-  data.frame(spp=prolesAll$species, matsys=prolesAll$mating.sys, ssd=prolesAll$ssd, dichro=prolesAll$dichro, care=prolesAll$care)
proles_dist <- cbind(proles[match(rownames(sr_dist), proles$spp),], sexroledist=log10(sr_dist[,1]))

proles_dist_all <- data.frame(spp=prolesAll$species, sexroledist=NA)
proles_dist_all[match(rownames(sr_dist), proles$spp), "sexroledist"] <- sr_dist[,1]

#Add sex role distance to the dataframe
dat<-merge(dat, proles_dist_all, by.x="species", all.x=TRUE)

#Bivariate correlations between sex role components

#Select sex role components

Sexrole<-subset(dat, select=c("species", "ssd", "dichro", "mating.sys", "care"))
Sexrole<-na.omit(Sexrole)

#Cut tree to match dataframe

to.drop<-setdiff(phy$tip.label, Sexrole$species)
phy.cut<-drop.tip(phy, to.drop)


#Use the estimated multivariate-lambda from the phylogenetic PCA to rescale the covariance matrix

library(corrplot)

invC<-solve(vcv.phylo(rescale(phy.cut, "lambda",0.72877)))

Sexrole<-Sexrole[match(phy.cut$tip.label,Sexrole$species),]
#Use dataframe of sex role components to estimate bivariate correlations as below (Changed code to create square matrix):

nsps<-length(phy.cut$tip.label)

mycorvec<-c()

i1<-c("ssd","dichro","mating.sys","care")
j2<-c("ssd","dichro","mating.sys","care")

for (i in 1:4) {
  for (j in 1:4){
    
    x<-Sexrole[,i1[i]]
    y<-Sexrole[,j2[j]]
    
    mean.x<-colSums(invC%*%x)/sum(invC)
    mean.y<-colSums(invC%*%y)/sum(invC)
    vector.ones<-as.matrix(rep(1,nsps))
    var.x<-t(x-vector.ones%*%mean.x)%*%invC%*%(x-vector.ones%*%mean.x)/(nsps-1)
    var.y<-t(y-vector.ones%*%mean.x)%*%invC%*%(y-vector.ones%*%mean.x)/(nsps-1)
    cor.xy<-(t(x-vector.ones%*%mean.x)%*%invC%*%(y-vector.ones%*%mean.x)/(nsps-1))/sqrt(var.x*var.y)
    mycorvec<-c(mycorvec, cor.xy)
    
  }
}

mycorvec2<-mycorvec[c(1,2,3,5,6,9)]
names(mycorvec2)<-c("ssd-dichro","ssd-mating.sys","ssd-care","dichro-mating.sys","dichro-care","mating.sys-care")
mycorvec2
#Create a square matrix to use for correlation table plotting
mycor.mat<-matrix(mycorvec, nrow=4, ncol=4)

#Add row and column names for the correlation table
colnames(mycor.mat)<-c("SSD","Sexual dichromatism","Mating system","Parental care")
rownames(mycor.mat)<-c("SSD","Sexual dichromatism","Mating system","Parental care")

#Define color palette for background plotting and plot
pdf(file="CorrelationPlot.pdf", height=8, width = 8)
col1 <- colorRampPalette(c("white","yellow","orange","red","#7F0000"))
corrplot(mycor.mat, type="upper", is.corr=FALSE, col=col1(6), addCoef.col = "grey", diag=FALSE, tl.col="black", cl.pos="n")
dev.off

#c("white","yellow", "#FF7F00", "red","#7F0000")


## Adding Extra pair copulation

#Include epb as one more potential "predictor" of sex roles, within the social environment hypothesis. Test whether epb influences sex roles. 

#Association between SSD and extra-pair copulations, there are 269 species with data for both:

SSD.epb<-as.data.frame(cbind(Sexrole$ssd, epb.data))
names(SSD.epb)<-c("SSD", "species", "epb")
SSD.epb<-na.omit(SSD.epb)
to.drop4<-setdiff(phy$tip.label, SSD.epb$species)
phy.cutSSD<-drop.tip(phy, to.drop4)
comp.dat.SSD<-comparative.data(phy.cutSSD, SSD.epb, names.col=species, vcv = TRUE, warn.dropped = TRUE)
m5<-pgls(SSD~epb, data=comp.dat.SSD, lambda = 'ML')
summary(m5)

#Association between dichromatism and extra-pair copulations, there are 284 species with data for both variables:

dichro.epb<-as.data.frame(cbind(Sexrole$dichro, epb.data))
names(dichro.epb)<-c("dichro", "species", "epb")
dichro.epb<-na.omit(dichro.epb)
to.drop5<-setdiff(phy$tip.label, dichro.epb$species)
phy.cutdichro<-drop.tip(phy, to.drop5)
comp.dat.dichro<-comparative.data(phy.cutdichro, dichro.epb, names.col = species, vcv=TRUE, warn.dropped = TRUE)
m6<-pgls(dichro~epb, data=comp.dat.dichro, lambda = 'ML')
summary(m6)

#Check residuals of the model.   
par(mfrow=c(2,2))
plot(m6)

#Relationship between parental care and extra-pair broods.   

care.epb<-as.data.frame(cbind(Sexrole$care, epb.data))
names(care.epb)<-c("care", "species", "epb")
care.epb<-na.omit(care.epb)
to.drop6.1<-setdiff(phy$tip.label, care.epb$species)
phy.cutcare<-drop.tip(phy, to.drop6.1)
comp.dat.care.epb<-comparative.data(phy.cutcare, care.epb, names.col = species, vcv=TRUE, warn.dropped = TRUE)
m6.1<-pgls(care~epb, data=comp.dat.care.epb, lambda = 'ML')
summary(m6.1)


#Tests of the different hypotheses.
#First life-history hypotheses

#Phylogenetic multiple regression models for life-history 
#traits, including clutch size, incubation duration, female mass, developmental mode
#and the interactions between female mass and clutch size, and female mass and incubation.

#First SSD

ssd.lh<-subset(dat, select=c(species, clutch.size, inc.dur, dev.mod.new, f.mass, ssd))
ssd.lh<-na.omit(ssd.lh)
to.drop<-setdiff(phy$tip.label, ssd.lh$species)
phy.ssd.lh<-drop.tip(phy, to.drop)
rownames(ssd.lh)<-ssd.lh$species
m.ssd.lh<-phylolm(ssd~clutch.size+inc.dur+dev.mod.new+log10(f.mass)+clutch.size*log10(f.mass)+inc.dur*log10(f.mass),
                  data=ssd.lh, phy=phy.ssd.lh, model="lambda")
summary(m.ssd.lh)

#Now dichromatism
dichro.lh<-subset(dat, select=c(species, clutch.size, inc.dur, dev.mod.new, f.mass, dichro))
dichro.lh<-na.omit(dichro.lh)
to.drop<-setdiff(phy$tip.label,dichro.lh$species)
phy.dichro.lh<-drop.tip(phy, to.drop)
rownames(dichro.lh)<-dichro.lh$species
m.dichro.lh<-phylolm(dichro~clutch.size+inc.dur+dev.mod.new+log10(f.mass)+clutch.size*log10(f.mass)+inc.dur*log10(f.mass),
                     data=dichro.lh, phy=phy.dichro.lh, model="lambda")
summary(m.dichro.lh)

#Now mating system
mate.lh<-subset(dat, select=c(species, clutch.size, inc.dur, dev.mod.new, f.mass, mating.sys))
mate.lh<-na.omit(mate.lh)
to.drop<-setdiff(phy$tip.label, mate.lh$species)
phy.mate.lh<-drop.tip(phy, to.drop)
rownames(mate.lh)<-mate.lh$species
m.mate.lh<-phylolm(mating.sys~clutch.size+inc.dur+dev.mod.new+log10(f.mass)+clutch.size*log10(f.mass)+inc.dur*log10(f.mass),
                   data=mate.lh, phy=phy.mate.lh, model="lambda")
summary(m.mate.lh)

#Now parental care
care.lh<-subset(dat, select=c(species, clutch.size, inc.dur, dev.mod.new, f.mass, care))
care.lh<-na.omit(care.lh)
to.drop<-setdiff(phy$tip.label, care.lh$species)
phy.care.lh<-drop.tip(phy, to.drop)
rownames(care.lh)<-care.lh$species
m.care.lh<-phylolm(care~clutch.size+inc.dur+dev.mod.new+log10(f.mass)+clutch.size*log10(f.mass)+inc.dur*log10(f.mass),
                   data=care.lh, phy=phy.care.lh, model="lambda")
summary(m.care.lh)

#Finally extent of sex role bias:   
bias.lh<-subset(dat, select=c(species, clutch.size, inc.dur, dev.mod.new, f.mass, sexroledist))
bias.lh<-na.omit(bias.lh)
to.drop<-setdiff(phy$tip.label,bias.lh$species)
phy.bias.lh<-drop.tip(phy, to.drop)
rownames(bias.lh)<-bias.lh$species
m.bias.lh<-phylolm(sexroledist~clutch.size+inc.dur+dev.mod.new+log10(f.mass)+clutch.size*log10(f.mass)+inc.dur*log10(f.mass),
                   data=bias.lh, phy=phy.bias.lh, model="lambda")
summary(m.bias.lh)

#Now the phylogenetic linear models testing the remaining hypotheses
#Predictors: Adult survival, ASR, Coloniality (Grouping), Extra-pair broods (epb)

#First SSD

ssd.surv<-subset(dat, select=c(species, survival, ssd))
ssd.surv<-na.omit(ssd.surv)
rownames(ssd.surv)<-ssd.surv$species
to.drop<-setdiff(phy$tip.label, ssd.surv$species)
phy.ssd.surv<-drop.tip(phy, to.drop)
m.ssd.surv<-phylolm(ssd~survival, data=ssd.surv, phy=phy.ssd.surv, model="lambda")
summary(m.ssd.surv)

ssd.asr<-subset(dat, select=c(species, asr, ssd))
ssd.asr<-na.omit(ssd.asr)
rownames(ssd.asr)<-ssd.asr$species
to.drop<-setdiff(phy$tip.label,ssd.asr$species)
phy.ssd.asr<-drop.tip(phy, to.drop)
m.ssd.asr<-phylolm(ssd~asr, data=ssd.asr, phy=phy.ssd.asr, model="lambda")
summary(m.ssd.asr)

ssd.col<-subset(dat, select=c(species, Grouping, ssd))
ssd.col<-na.omit(ssd.col)
rownames(ssd.col)<-ssd.col$species
to.drop<-setdiff(phy$tip.label,ssd.col$species)
phy.ssd.col<-drop.tip(phy, to.drop)
m.ssd.col<-phylolm(ssd~Grouping, data=ssd.col, phy=phy.ssd.col, model="lambda")
summary(m.ssd.col)

ssd.epb<-subset(dat, select=c(species, epb, ssd))
ssd.epb<-na.omit(ssd.epb)
rownames(ssd.epb)<-ssd.epb$species
to.drop<-setdiff(phy$tip.label, ssd.epb$species)
phy.ssd.epb<-drop.tip(phy, to.drop)
m.ssd.epb<-phylolm(ssd~epb, data=ssd.epb, phy=phy.ssd.epb, model="lambda")
summary(m.ssd.epb)

#Now dichromatism

dichro.surv<-subset(dat, select=c(species, survival, dichro))
dichro.surv<-na.omit(dichro.surv)
rownames(dichro.surv)<-dichro.surv$species
to.drop<-setdiff(phy$tip.label, dichro.surv$species)
phy.dichro.surv<-drop.tip(phy, to.drop)
m.dichro.surv<-phylolm(dichro~survival, data=dichro.surv, phy=phy.dichro.surv, model="lambda")
summary(m.dichro.surv)

dichro.asr<-subset(dat, select=c(species, asr, dichro))
dichro.asr<-na.omit(dichro.asr)
rownames(dichro.asr)<-dichro.asr$species
to.drop<-setdiff(phy$tip.label,dichro.asr$species)
phy.dichro.asr<-drop.tip(phy,to.drop)
m.dichro.asr<-phylolm(dichro~asr, data=dichro.asr, phy=phy.dichro.asr, model="lambda")
summary(m.dichro.asr)

dichro.group<-subset(dat, select=c(species, Grouping, dichro))
dichro.group<-na.omit(dichro.group)
rownames(dichro.group)<-dichro.group$species
to.drop<-setdiff(phy$tip.label,dichro.group$species)
phy.dichro.group<-drop.tip(phy, to.drop)
m.dichro.group<-phylolm(dichro~Grouping, data=dichro.group, phy=phy.dichro.group, model="lambda")
summary(m.dichro.group)

dichro.epb<-subset(dat, select=c(species, epb, dichro))
dichro.epb<-na.omit(dichro.epb)
rownames(dichro.epb)<-dichro.epb$species
to.drop<-setdiff(phy$tip.label, dichro.epb$species)
phy.dichro.epb<-drop.tip(phy,to.drop)
m.dichro.epb<-phylolm(dichro~epb, data=dichro.epb, phy=phy.dichro.epb, model="lambda")
summary(m.dichro.epb)

#Now Mating system

mate.surv<-subset(dat, select=c(species, survival, mating.sys))
mate.surv<-na.omit(mate.surv)
rownames(mate.surv)<-mate.surv$species
to.drop<-setdiff(phy$tip.label,mate.surv$species)
phy.mate.surv<-drop.tip(phy,to.drop)
m.mate.surv<-phylolm(mating.sys~survival, data=mate.surv, phy=phy.mate.surv, model="lambda")
summary(m.mate.surv)

mate.asr<-subset(dat, select=c(species, asr, mating.sys))
mate.asr<-na.omit(mate.asr)
rownames(mate.asr)<-mate.asr$species
to.drop<-setdiff(phy$tip.label, mate.asr$species)
phy.mate.asr<-drop.tip(phy,to.drop)
m.mate.asr<-phylolm(mating.sys~asr, data=mate.asr, phy=phy.mate.asr, model="lambda")
summary(m.mate.asr)

mate.group<-subset(dat, select=c(species, Grouping, mating.sys))
mate.group<-na.omit(mate.group)
rownames(mate.group)<-mate.group$species
to.drop<-setdiff(phy$tip.label,mate.group$species)
phy.mate.group<-drop.tip(phy,to.drop)
m.mate.group<-phylolm(mating.sys~Grouping, data=mate.group, phy=phy.mate.group, model="lambda")
summary(m.mate.group)

mate.epb<-subset(dat, select=c(species, epb, mating.sys))
mate.epb<-na.omit(mate.epb)
rownames(mate.epb)<-mate.epb$species
to.drop<-setdiff(phy$tip.label,mate.epb$species)
phy.mate.epb<-drop.tip(phy, to.drop)
m.mate.epb<-phylolm(mating.sys~epb, data=mate.epb, phy=phy.mate.epb, model="lambda")
summary(m.mate.epb)

#Parental care

care.surv<-subset(dat, select=c(species, survival, care))
care.surv<-na.omit(care.surv)
rownames(care.surv)<-care.surv$species
to.drop<-setdiff(phy$tip.label,care.surv$species)
phy.care.surv<-drop.tip(phy,to.drop)
m.care.surv<-phylolm(care~survival, data=care.surv, phy=phy.care.surv, model="lambda")
summary(m.care.surv)

care.asr<-subset(dat, select=c(species, asr, care))
care.asr<-na.omit(care.asr)
rownames(care.asr)<-care.asr$species
to.drop<-setdiff(phy$tip.label,care.asr$species)
phy.care.asr<-drop.tip(phy,to.drop)
m.care.asr<-phylolm(care~asr, data=care.asr, phy=phy.care.asr, model="lambda")
summary(m.care.asr)

care.group<-subset(dat, select=c(species, Grouping, care))
care.group<-na.omit(care.group)
rownames(care.group)<-care.group$species
to.drop<-setdiff(phy$tip.label,care.group$species)
phy.care.group<-drop.tip(phy,to.drop)
m.care.group<-phylolm(care~Grouping, data=care.group, phy=phy.care.group, model="lambda")
summary(m.care.group)

care.epb<-subset(dat, select=c(species, epb, care))
care.epb<-na.omit(care.epb)
rownames(care.epb)<-care.epb$species
to.drop<-setdiff(phy$tip.label,care.epb$species)
phy.care.epb<-drop.tip(phy,to.drop)
m.care.epb<-phylolm(care~epb, data=care.epb, phy=phy.care.epb, model="lambda")
summary(m.care.epb)

#Finally sex role bias   

bias.surv<-subset(dat, select=c(species, survival, sexroledist))
bias.surv<-na.omit(bias.surv)
rownames(bias.surv)<-bias.surv$species
to.drop<-setdiff(phy$tip.label,bias.surv$species)
phy.bias.surv<-drop.tip(phy,to.drop)
m.bias.surv<-phylolm(sexroledist~survival, data=bias.surv, phy=phy.bias.surv, model="lambda")
summary(m.bias.surv)

bias.asr<-subset(dat, select=c(species, asr, sexroledist))
bias.asr<-na.omit(bias.asr)
rownames(bias.asr)<-bias.asr$species
to.drop<-setdiff(phy$tip.label,bias.asr$species)
phy.bias.asr<-drop.tip(phy,to.drop)
m.bias.asr<-phylolm(sexroledist~asr, data=bias.asr, phy=phy.bias.asr, model="lambda")
summary(m.bias.asr)

bias.epb<-subset(dat, select=c(species, epb, sexroledist))
bias.epb<-na.omit(bias.epb)
rownames(bias.epb)<-bias.epb$species
to.drop<-setdiff(phy$tip.label,bias.epb$species)
phy.bias.epb<-drop.tip(phy,to.drop)
m.bias.epb<-phylolm(sexroledist~epb, data=bias.epb, phy=phy.bias.epb, model="lambda")
summary(m.bias.epb)

#Analyses of the relationship between climate and sex roles

#Mating system and ecology

#Mating system and its relation with mean temperature during the breeding season:   
#These analyses were run with the caper package, for transparency it is presented 
#as such, it is much quicker if run with phylolm.

mate.eco<-subset(dat, select=c(species,mating.sys,mean_breed_temp,breed_temp_var, mean_breed_prec))
mate.eco<-na.omit(mate.eco)
to.drop<-setdiff(phy$tip.label,mate.eco$species)
mate.eco.phy<-drop.tip(phy,to.drop)
mate.eco.comp<-comparative.data(mate.eco.phy, mate.eco, names.col="species", vcv=TRUE, warn.dropped=TRUE)

m42<-pgls(mating.sys~mean_breed_temp, data=mate.eco.comp, lambda='ML')

plot(m42)

summary(m42)

plot(mate.eco$mean_breed_temp, mate.eco$mating.sys, xlab="Mean temperature breeding season", ylab="Mating system")



#Now mating system and its realtionship with variation in temperature during the breeding season:   

m43<-pgls(mating.sys~breed_temp_var, data=mate.eco.comp, lambda='ML')

plot(m43)

summary(m43)

plot(mate.eco$breed_temp_var, mate.eco$mating.sys, xlab="Variation in temperature during breeding season", ylab = "Mating system")


#Finally, the relationship between mating system and mean precipitation during the breeding season:   

m44<-pgls(Mating.sys~mean_breed_prec, data=mate.eco.comp, lambda='ML')

plot(m44)

summary(m44)

plot(mate.eco$mean_breed_prec, mate.eco$mating.sys, xlab="Mean precipitation breeding season", ylab="Mating system")

#Now parental care and climate
#First the relationship between parental care and mean breeding temperature:   

care.eco<-subset(dat, select=c(species,care,mean_breed_temp,coldest_breed_month_temp,hottest_breed_month_temp,breed_temp_var, mean_breed_prec))
care.eco<-na.omit(care.eco)
to.drop<-setdiff(phy$tip.label,care.eco$species)
care.eco.phy<-drop.tip(phy,to.drop)
row.names(care.eco)<-care.eco$species
#care.eco.comp<-comparative.data(care.eco.phy, care.eco, names.col="species", vcv=TRUE, warn.dropped=TRUE)

m45<-pgls(care~mean_breed_temp, data=care.eco.comp, lambda='ML')

summary(m45)

plot(m45)

m45.1<-phylolm(care~mean_breed_temp, data=care.eco, phy=care.eco.phy, model='lambda')

summary(m45.1)

plot(care.eco$mean_breed_temp, care.eco$care, xlab="Mean temperature breeding season", ylab = "Parental care")

#Now parental care and variation in temperature during the breeding season:   

m46<-pgls(care~breed_temp_var, data=care.eco.comp, lambda='ML')

plot(m46)

summary(m46)

m46.1<-phylolm(care~breed_temp_var, data=care.eco, phy=care.eco.phy, model='lambda')

summary(m46.1)

plot(care.eco$breed_temp_var, care.eco$care, xlab="Variation in temperature breeding season", ylab="Parental care")


#Finally the relationship between parental care and mean precipitation during the breeding season:   

m47<-pgls(care~mean_breed_prec, data=care.eco.comp, lambda='ML')

plot(m47)

summary(m47)

m47.1<-phylolm(care~mean_breed_prec, data=care.eco, phy=care.eco.phy, model='lambda')

summary(m47.1)

plot(care.eco$mean_breed_prec, care.eco$care, xlab="Mean precipitation breeding season", ylab="Parental Care")

#Now sexual dichromatism and climate
#First dichromatism versus mean breeding temperature:   

dichro.eco<-subset(dat, select=c(species,dichro,mean_breed_temp,coldest_breed_month_temp,hottest_breed_month_temp,breed_temp_var, mean_breed_prec))
dichro.eco<-na.omit(dichro.eco)
to.drop<-setdiff(phy$tip.label,dichro.eco$species)
dichro.eco.phy<-drop.tip(phy,to.drop)
dichro.eco.comp<-comparative.data(dichro.eco.phy, dichro.eco, names.col="species", vcv=TRUE, warn.dropped=TRUE)
m39<-pgls(dichro~mean_breed_temp, data=dichro.eco.comp, lambda='ML')

plot(m39)

summary(m39)

plot(dichro.eco$mean_breed_temp, dichro.eco$dichro, xlab="Mean breeding temperature", ylab="Sexual dichromatism")

#Now dichromatism versus variation in temperature during the breeding season:   

m40<-pgls(dichro~breed_temp_var, data=dichro.eco.comp, lambda='ML')

plot(m40)

summary(m40)

plot(dichro.eco$breed_temp_var, dichro.eco$dichro, xlab="Variation in temperature during breeding season", ylab = "Sexual dichromatism")


#Finally, dichromatism versus mean precipitation during the breeding season:   

m41<-pgls(dichro~mean_breed_prec, data=dichro.eco.comp, lambda='ML')

plot(m41)

summary(m41)

plot(dichro.eco$mean_breed_prec, dichro.eco$dichro, xlab="Mean precipitation during breeding season", ylab = "Sexual dichromatism")

#Finally sexual size dimorphism and climate
#SSD versus mean breeding temperature


ssd.eco<-subset(dat, select=c(species,ssd,mean_breed_temp,coldest_breed_month_temp,hottest_breed_month_temp,breed_temp_var, mean_breed_prec))
ssd.eco<-na.omit(ssd.eco)
to.drop<-setdiff(phy$tip.label,ssd.eco$species)
ssd.eco.phy<-drop.tip(phy,to.drop)
ssd.eco.comp<-comparative.data(ssd.eco.phy, ssd.eco, names.col="species", vcv=TRUE, warn.dropped=TRUE)
m36<-pgls(ssd~mean_breed_temp, data=ssd.eco.comp, lambda='ML')

plot(m36)

summary(m36)

plot(ssd.eco$mean_breed_temp, ssd.eco$ssd, xlab="Mean breeding temperature", ylab="Sexual size dimor")


# Now SSD versus variation in temperature during the breeding season:   

m37<-pgls(ssd~breed_temp_var, data=ssd.eco.comp, lambda='ML')

plot(m37)

summary(m37)

plot(ssd.eco$breed_temp_var, ssd.eco$ssd, xlab="Variation in breeding season temperature", ylab="Sexual size dimor")



#Now SSD versus mean precipitation during the breeding season:   

m38<-pgls(ssd~mean_breed_prec, data=ssd.eco.comp, lambda='ML')

plot(m38)

summary(m38)

plot(ssd.eco$mean_breed_prec, ssd.eco$ssd, xlab="Mean precipitation during breeding season", ylab="Sexual size dimor")


############################################
#Plots for Supplementary materials
############################################

pdf(file="Histograms plot")

par(mfrow=c(2,2))

hist(dat$ssd, xlab="Sexual size dimorphism", main=NULL)

hist(dat$dichro, xlab="Sexual dichromatism", ylab=NULL, main=NULL)

hist(dat$mating.sys, xlab="Mating system", main=NULL)

hist(dat$care, xlab="Parental care", ylab=NULL, main=NULL)

dev.off()

pdf(file="ASR models plot", width = 9, height = 9)

par(mfrow=c(3,2))

plot(dat$ssd~dat$asr, xlab="Adult sex ratio", ylab="SSD")
abline(a=m.ssd.asr$coefficients[1], b=m.ssd.asr$coefficients[2], col="red",lwd=2)

plot(dat$dichro~dat$asr, xlab="Adult sex ratio", ylab="Sexual dichromatism")
abline(a=m.dichro.asr$coefficients[1], b=m.dichro.asr$coefficients[2], col="red",lwd=2)

plot(dat$mating.sys~dat$asr, xlab="Adult sex ratio", ylab="Mating system")
abline(a=m.mate.asr$coefficients[1], b=m.mate.asr$coefficients[2], col="red", lwd=2)

plot(dat$care~dat$asr, xlab="Adult sex ratio", ylab="Parental care")
abline(a=m.care.asr$coefficients[1], b=m.care.asr$coefficients[2], col="red", lwd=2)

plot(dat$sexroledist~dat$asr, xlab="Adult sex ratio", ylab="Extent of sex role bias", ylim=c(0,80))
abline(a=m.bias.asr$coefficients[1], b=m.bias.asr$coefficients[2], col="red",lwd=2)

dev.off()


