
############################################################
####                                                    ####  
####  NRES 746, Student-led topic #6                    ####
####                                                    ####
############################################################


############################################################
####  Ordination!                                       ####
############################################################



suppressWarnings(library(vegan))
suppressWarnings(library(MASS))
data(varespec)
head(varespec)
sum(varespec==0)/(nrow(varespec)*ncol(varespec)) # 42% of the data are zeroes
dim(varespec) # 44 different species, not distributed evenly among sites


data(mtcars)
head(mtcars)
dim(mtcars)
cars<-mtcars[,-c(8:9)] # remove categorical variables

mtcars.pca<-prcomp(mtcars[,c(1:7,10,11)], center = TRUE,scale. = TRUE) 
summary(mtcars.pca)


str(mtcars.pca) # Look at PCA object


biplot(mtcars.pca) # This is a little messy with sample names


suppressWarnings(library(devtools))
#install_github("vqv/ggbiplot")
suppressWarnings(library(ggbiplot))
ggbiplot(mtcars.pca)+
  theme_bw()+ # removes gray background
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+ # removes gridlines
  ggtitle("Cars PCA")


ggbiplot(mtcars.pca, labels=rownames(mtcars))+
  theme_bw()+ # removes gray background
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+ # removes gridlines
  ggtitle("Cars PCA")+
  xlab("PC1")+
  ylab("PC2")


mtcars.country <- c(rep("Japan", 3), rep("US",4), rep("Europe", 7),rep("US",3), "Europe", rep("Japan", 3), rep("US",4), rep("Europe", 3), "US", rep("Europe", 3))

ggbiplot(mtcars.pca,ellipse=TRUE,  labels=rownames(mtcars), groups=mtcars.country, obs.scale=1, var.scale=1)+ # ellipse=T puts ellipses around groups, obs and var scale help space out the plot a bit
  scale_colour_manual(name="Origin", values= c("black", "red3", "blue"))+ # change default grouping colors
  theme_bw()+ # removes gray background
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+ # removes gridlines
  ggtitle("Cars PCA")+
  xlab("PC1")+
  ylab("PC2")

# Simplify the plot again to remove car names.
ggbiplot(mtcars.pca,ellipse=TRUE, groups=mtcars.country, obs.scale=1, var.scale=1)+ # ellipse=T puts ellipses around groups, obs and var scale help space out the plot a bit
  scale_colour_manual(name="Origin", values= c("black", "red3", "blue"))+ # change default grouping colors
  theme_bw()+ # removes gray background
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+ # removes gridlines
  ggtitle("Cars PCA")+
  xlab("PC1")+
  ylab("PC2")


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


library(scatterplot3d)
d=scatterplot3d(0:10,0:10,0:10,type="n",xlab="Abundance of Species 1",
                ylab="Abundance of Species 2",
                zlab="Abundance of Species 3")
d$points3d(5,5,0); text(d$xyz.convert(5,5,0.5),labels="community A")
d$points3d(3,3,3); text(d$xyz.convert(3,3,3.5),labels="community B")
d$points3d(0,5,5); text(d$xyz.convert(0,5,5.5),labels="community C")


data(varespec)
head(varespec)
# Create a matrix of dissimilarity indices
vare.dis<-vegdist(varespec) # default is Bray-Curtis dissimilarity index
head(vare.dis)


vare.mds<-metaMDS(varespec,k=2,trymax=100,trace=F) #k = # of dimensions you want for the solution, trymax is the max number of iterations to try to converge on a solution with the lowest stress
vare.mds


stressplot(vare.mds,vare.dis)


ordiplot(vare.mds,type="t") #shows site ordination scores in species space based on species averages


data(varechem)
head(varechem)
rankindex(scale(varechem), varespec, c("euc","man","bray","jac","kul"))


ordisurf(vare.mds~N,varechem)


ef<-envfit(vare.mds,varechem,permu=999)
ef


plot(vare.mds,display="sites")
plot(ef,p.max=0.1) #limit to most significant variables

