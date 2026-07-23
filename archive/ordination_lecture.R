mtcars.pca <- prcomp(mtcars[,c(1:7,10,11)],center=TRUE,scale.=TRUE)
summary(mtcars.pca)
str(mtcars.pca)
x.var <- mtcars.pca$sdev ^ 2
x.pvar <- x.var/sum(x.var)
plot(x.pvar,xlab="Principal component", ylab="Proportion of variance explained", ylim=c(0,1), type='b')
mtcars.pca
library(devtools)
install_github('vqv/ggbiplot')
library(ggbiplot)
ggbiplot(mtcars.pca)
ggbiplot(mtcars.pca,labels=rownames(mtcars))
mtcars.country <- c(rep("Japan",3),rep("US",4),rep("Europe",7),rep("US",3),"Europe",rep("Japan",3),rep("US",4),rep("Europe",3),"US",rep("Europe",3))
ggbiplot(mtcars.pca,ellipse=TRUE,labels=rownames(mtcars),groups=mtcars.country)
ggbiplot(mtcars.pca,choices=c(3,4),ellipse=TRUE,labels=rownames(mtcars),groups=mtcars.country)
TeslaMarsCar <- c(1000,60,50,500,0,0.5,2.5,0,1,0,0)

mtcarsplus <- rbind(mtcars, TeslaMarsCar)
mtcars.countryplus <- c(mtcars.country, "Mars")

mtcarsplus.pca <- prcomp(mtcarsplus[,c(1:7,10,11)], center = TRUE,scale. = TRUE)
summary(mtcarsplus.pca)

ggbiplot(mtcarsplus.pca, obs.scale = 1, var.scale = 1, ellipse = TRUE, circle = FALSE, var.axes=TRUE, labels=c(rownames(mtcars), "TeslaMarsCar"), groups=mtcars.countryplus)+
  scale_colour_manual(name="Origin", values= c("forest green", "red", "violet", "dark blue"))+
  ggtitle("PCA of mtcars dataset, with extra sample added")+
  theme_minimal()+
  theme(legend.position = "bottom")
tmc.sc <- scale(t(TeslaMarsCar[c(1:7,10,11)]), center= mtcars.pca$center)
tmc.pred <- tmc.sc %*% mtcars.pca$rotation

mtcars.plusproj.pca <- mtcars.pca
mtcars.plusproj.pca$x <- rbind(mtcars.plusproj.pca$x, tmc.pred)

ggbiplot(mtcars.plusproj.pca, obs.scale = 1, var.scale = 1, ellipse = TRUE, circle = FALSE, var.axes=TRUE, labels=c(rownames(mtcars), "TeslaMarsCar"), groups=mtcars.countryplus)+
  scale_colour_manual(name="Origin", values= c("forest green", "red3", "violet", "dark blue"))+
  ggtitle("PCA of mtcars dataset, with extra sample projected")+
  theme_minimal()+
  theme(legend.position = "bottom")
fa.data <- read.csv('EFA.csv',header=TRUE)
head(fa.data)

library(psych)
#install.packages('psych')
library(GPArotation)
#install.packages('GPArotation')
parallel <- fa.parallel(fa.data, fm='minres',fa='fa')
fourfactor <- fa(fa.data,nfactors=4,rotate='oblimin',fm='minres')
print(fourfactor$loadings,cutoff=0.35)
fa.diagram(fourfactor)
print(fourfactor)
plot(0:10,0:10,type="n",axes=F,xlab="Abundance of Species 1",ylab="") 
axis(1)
points(5,0); text(5.5,0.5,labels="community A")
points(3,0); text(3.2,0.5,labels="community B")
points(0,0); text(0.8,0.5,labels="community C")
plot(0:10,0:10,type="n",xlab="Abundance of Species 1",
     ylab="Abundance of Species 2")
points(5,5); text(5,4.5,labels="community A")
points(3,3); text(3,3.5,labels="community B")
points(0,5); text(0.8,5.5,labels="community C")
# install.packages("scatterplot3d")
library(scatterplot3d)
d=scatterplot3d(0:10,0:10,0:10,type="n",xlab="Abundance of Species 1", ylab="Abundance of Species 2",zlab="Abundance of Species 3"); 
d$points3d(5,5,0); text(d$xyz.convert(5,5,0.5),labels="community A")
d$points3d(3,3,3); text(d$xyz.convert(3,3,3.5),labels="community B")
d$points3d(0,5,5); text(d$xyz.convert(0,5,5.5),labels="community C")
#install.packages("vegan")
library(vegan)
set.seed(2)
community_matrix <- matrix(
   sample(1:100,300,replace=T),nrow=10,
   dimnames=list(paste("community",1:10,sep=""),paste("sp",1:30,sep="")))

example_NMDS <- metaMDS(community_matrix,k=2)
example_NMDS <- metaMDS(community_matrix,k=3,trymax=100)
stressplot(example_NMDS)
plot(example_NMDS)
ordiplot(example_NMDS,type="n")
orditorp(example_NMDS,display="species",col="red",air=0.01)
orditorp(example_NMDS,display="sites",cex=1.25,air=0.01)
treat <- c(rep("Treatment1",5),rep("Treatment2",5))
ordiplot(example_NMDS,type="n")
ordiellipse(example_NMDS,groups=treat,draw="polygon",col="grey90",label=F)
orditorp(example_NMDS,display="species",col="red",air=0.01)
orditorp(example_NMDS,display="sites",col=c(rep("green",5),rep("blue",5)),
   air=0.01,cex=1.25)
colors=c(rep("red",5),rep("blue",5))
ordiplot(example_NMDS,type="n")

for(i in unique(treat)) {
  ordiellipse(example_NMDS$point[grep(i,treat),],draw="polygon",
   groups=treat[treat==i],col=colors[grep(i,treat)],label=F) } 
orditorp(example_NMDS,display="species",col="red",air=0.01)
orditorp(example_NMDS,display="sites",col=c(rep("green",5),
  rep("blue",5)),air=0.01,cex=1.25)
ordiplot(example_NMDS,type="n")
orditorp(example_NMDS,display="species",col="red",air=0.01)
orditorp(example_NMDS,display="sites",col=c(rep("green",5),rep("blue",5)),
         air=0.01,cex=1.25)
ordicluster(example_NMDS,hclust(vegdist(community_matrix,"bray"))) 
# Define random elevations for previous example
elevation=runif(10,0.5,1.5)
# Use the function ordisurf to plot contour lines
ordisurf(example_NMDS,elevation,main="",col="forestgreen")
# Finally, display species on plot
orditorp(example_NMDS,display="species",col="grey30",air=0.1,
   cex=1)
