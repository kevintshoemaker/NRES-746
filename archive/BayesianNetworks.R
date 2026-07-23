
############################################################
####                                                    ####  
####  NRES 746, Student-led topic #4                    ####
####                                                    ####
############################################################


############################################################
####  Bayesian networks                                 ####
############################################################



#install.packages("bnlearn")
library(bnlearn)


data(coronary)
head(coronary)


bn_df <- data.frame(coronary)
res <- hc(bn_df)
plot(res)


res$arcs <- res$arcs[-which((res$arcs[,'from'] == "M..Work" & res$arcs[,'to'] == "Family")),]


fittedbn <- bn.fit(res, data = bn_df)
print(fittedbn$Proteins)


cpquery(fittedbn, event = (Proteins=="<3"), evidence = ( Smoking=="no") )

cpquery(fittedbn, event = (Proteins=="<3"), evidence = ( Smoking=="no" & Pressure==">140" ) )

cpquery(fittedbn, event = (Pressure==">140"), evidence = ( Proteins=="<3" ) )


#install.packages("sparsebn")
library(sparsebn)


# source("http://bioconductor.org/biocLite.R")
# biocLite("Rgraphviz")   
library("Rgraphviz")


setPlotPackage("graph")   
plot(cytometryContinuous$dag)


setPlotPackage("igraph")
plot(cytometryContinuous$dag,
layout = igraph::layout_(to_igraph(cytometryContinuous$dag),
                         igraph::in_circle()),
vertex.label = names(cytometryContinuous$dag),
vertex.size = 30,
vertex.label.color = gray(0),
vertex.color = colors(),
edge.color = "red",
edge.arrow.size = 0.3)

