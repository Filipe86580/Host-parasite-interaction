################################################################################
##### SET THE STAGE 
################################################################################


#Set the working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Delete all previous objects
#rm(list= ls())

#Load the required packages
library(igraph)
library(bipartite)
library(Rmisc)
library(vegan)
library(gdata)
library(ggplot2)
library(gridExtra)
library(grid)

#Load some custom-made functions
source("RestNullModel.R")
source("PosteriorProb.R")
source("MyDiamond.R")


################################################################################
##### PROCESS THE NETWORK PEP
################################################################################


#Import the network
data_all <- as.matrix(read.csv2("data/Pasta2.csv", header = T, row.names=1))

#Inspect the network
class(data_all)
data_all
dim(data_all)
min(data_all)
max(data_all)

#Plot the matrix
visweb(data_all)

#Plot and export the graph
tiff(filename= "figures/data_all.tif", res= 300, height= 1500, width= 1700)
par(mfrow=c(1,1),mar=c(1,1,1,5))
visweb(data_all)
par(mfrow=c(1,1))
dev.off()

#Convert the network to igraph format
data2_all <- graph_from_incidence_matrix(data_all, directed = F, weighted = TRUE) 

#Inspect object
class(data2_all)
data2_all
E(data2_all)
V(data2_all)$name


#Inform which nodes represent which taxonomic groups
V(data2_all)$set[1:nrow(data_all)] = c(rep("Host", 80))
V(data2_all)$set[(nrow(data_all)+1):(nrow(data_all)+ncol(data_all))] = "Para"

#Set seed
set.seed(123)


################################################################################
##### DRAW THE NETWORK PEP
################################################################################


#Set layout
lay1_all <- layout_nicely(data2_all)

#Set edge mode and width
E(data2_all)$arrow.mode = 0
E(data2_all)$width = E(data2_all)$weight/5+1

#Import "diamond" vertex shape
source("MyDiamond.R")
add_shape("diamond", clip=shapes("circle")$clip,
          plot=MyDiamond)

#Set vertex shapes
V(data2_all)$shape = V(data2_all)$set 

V(data2_all)$shape = gsub("Host","square",V(data2_all)$shape)

V(data2_all)$shape = gsub("Para","circle",V(data2_all)$shape)

#Calculate DIRTLPAwb+ modularity and save the output as a data frame and a list
data.mod_all <- computeModules(data_all, method = "DormannStrauss")

data.mod_all
data.modules_all <- module2constraints(data.mod_all)
data.modules_all

data.df_all <- data.frame(c(rownames(data_all), colnames(data_all)), data.modules_all) 
data.df_all
colnames(data.df_all) <- c("vertices", "modules")

data.df_all
data.list_all <- split(data.df_all$vertices, data.df_all$modules)
data.list_all

##Set node and cloud colors by modularity
colors_all <- rainbow(length(data.list_all), alpha = 0.3, s = 1, v = 0.8)
V(data2_all)$color <- colors_all[data.df_all$modules]
clouds_all = colors_all

#Plot and export the graph
tiff(filename= "figures/all_FIG1.tif", res= 300, height= 3000, width= 5100)
par(mfrow=c(1,1),mar=c(1,1,1,5))
plot(data2_all,
     col = V(data2_all)$color,
     mark.groups = data.list_all,
     mark.border = "lightgrey", 
     mark.col = adjustcolor(clouds_all, alpha = 0.5),
     vertex.size = 10,
     vertex.label = V(data2_all)$name,
     vertex.label.font= 3,
     vertex.label.color = "black",
     vertex.label.cex = 1,
     vertex.frame.color = NA,
     edge.color = adjustcolor("grey", alpha.f = .5),
     edge.curved = 0.3,
     edge.width = 3,
     layout=lay1_all)
#legend(x = 0.9,y = 1.0, legend = c("Host", "Parasite"),
#      pch = c(18,15),  title = "Taxon",
#     text.col = "gray20", title.col = "black",
#     box.lwd = 0, cex = 2, col = c("grey", "grey", "grey"))
par(mfrow=c(1,1))
dev.off()


################################################################################
##### NETWORK LEVEL ANALYSIS (TOPOLOGY) PEP
################################################################################

#Set the number of permutations to be used in all null model analyses
#Here we set a low value just for testing the script. In our paper we've set it
# to 1,000.
permutations <- 1000

#Generate randomized matrices
nulls_all <- nullmodel(data_all, N=permutations, method="vaznull")


##### MODULARITY



Mod_all<-data.mod_all

#Extract module membership
Part_all  <- bipartite::module2constraints(Mod_all)
row.Part_all <- Part_all [1:nrow(data_all)]
col.Part_all <- Part_all [(nrow(data_all)+1):(nrow(data_all)+ncol(data_all))]

#Calculate metric for the randomized networks
nullmod_all <- sapply(nulls_all, computeModules, method = "DormannStrauss")
modnull_all <- sapply(nullmod_all, function(x) x@likelihood)
(Mod_all@likelihood - mean(modnull_all))/sd(modnull_all) # Z value
Mod.sig_all <- sum(modnull_all>(Mod_all@likelihood)) / length(modnull_all) # p value
Mod.sig_all

#Plot the observed value against the distribution of randomized values
plot(density(modnull_all ), main="Observed vs. randomized",
     xlim=c(min((Mod_all@likelihood), min(modnull_all )), 
            max((Mod_all@likelihood), max(modnull_all ))))
abline(v=Mod_all@likelihood, col="red", lwd=2, xlab="")


Mod_all@likelihood #observed
mean(modnull_all ) #randomized mean
sd(modnull_all ) #randomized SD
(Mod_all@likelihood - mean(modnull_all ))/sd(modnull_all ) # Z-value
sum(modnull_all >(Mod_all@likelihood)) / length(modnull_all ) #randomized > observed
sum(modnull_all <(Mod_all@likelihood)) / length(modnull_all ) #randomized < observed


##### SPECIALIZATION 

#Calculate metric for the original network
Spec <- networklevel(data_all, index="H2")
class(Spec)

#Calculate metric for the randomized networks
randomized.Spec <- unlist(sapply(nulls_all, networklevel, index="H2"))
(Spec - mean(randomized.Spec))/sd(randomized.Spec) # Z value
Spec.sig <- sum(randomized.Spec>Spec)/length(randomized.Spec) # p value
Spec.sig

#Plot the observed value against the distribution of randomized values
plot(density(randomized.Spec), main="Observed vs. randomized",
     xlim=c(min((Spec), min(randomized.Spec)), 
            max((Spec), max(randomized.Spec))))
abline(v=Spec, col="red", lwd=2, xlab="")

Spec #observed
mean(randomized.Spec) #randomized mean
sd(randomized.Spec) #randomized SD
(Spec - mean(randomized.Spec))/sd(randomized.Spec) # Z-value
sum(randomized.Spec>(Spec)) / length(randomized.Spec) #randomized > observed
sum(randomized.Spec<(Spec)) / length(randomized.Spec) #randomized < observed


##### NESTEDNESS

#Calculate metric for the original network
Nest_all<- networklevel(data_all, index="weighted NODF")

#Calculate metric for the randomized networks
randomized.Nest_all <- unlist(sapply(nulls_all, networklevel, index="weighted NODF"))
(Nest_all - mean(randomized.Nest_all))/sd(randomized.Nest_all) # Z value
Nest.sig_all <- sum(randomized.Nest_all>Nest_all)/length(randomized.Nest_all) # p value
Nest.sig_all

#Plot the observed value against the distribution of randomized values
plot(density(randomized.Nest_all), main="Observed vs. randomized",
     xlim=c(min((Nest_all), min(randomized.Nest_all)), 
            max((Nest_all), max(randomized.Nest_all))))
abline(v=Nest_all, col="red", lwd=2, xlab="")

Nest_all #observed
mean(randomized.Nest_all) #randomized mean
sd(randomized.Nest_all) #randomized SD
(Nest_all - mean(randomized.Nest_all))/sd(randomized.Nest_all) # Z-value
sum(randomized.Nest_all>(Nest_all)) / length(randomized.Nest_all) #randomized > observed
sum(randomized.Nest_all<(Nest_all)) / length(randomized.Nest_all) #randomized < observed


##### COMPOUND TOPOLOGY 

#Calculate the desired nestedness metric (here WNODA) for the original network.
obs.com_all <- unlist(bipartite::nest.smdm(x = data_all, 
                                       constraints = Part_all, #Input the modular structured recovered from step 2
                                       weighted = T, #By considering the edge weights, you are choosing WNODA
                                       decreasing = "abund"))

#Check the scores
obs.com_all

#Calculate constrained interaction probabilities considering the network's modular structure
Pij_all <- PosteriorProb(M = data_all, 
                     R.partitions = row.Part_all, #Input the modular structured recovered from step 2
                     C.partitions = col.Part_all, #Input the modular structured recovered from step 2
                     Prior.Pij = "degreeprob", #Choose the null model
                     Conditional.level = "modules") #Choose the kind of constraints

#Check what those probabilities look like
Pij_all
dim(Pij_all)

#Generate randomized networks with the null model of your choice, considering the interaction probabilities calculated before. 
nulls.com_all <- RestNullModel(M = data_all, 
                           Pij.Prob = Pij_all, #Recover the probabilities calculated in the previous command
                           Numbernulls = permutations, #This step may take long, so start experimenting with low values
                           Print.null = T, 
                           allow.degeneration = F, #Choose whether you allow orphan rows and columns to be removed or not
                           return.nonrm.species = F, 
                           connectance = T, byarea = T, 
                           R.partitions = row.Part_all, 
                           C.partitions = col.Part_all)

#Calculate the nestedness within and between modules
rest.nest_all <- nest.smdm(data_all, constraints = Part_all, 
                       weighted = T, 
                       decreasing = "abund", 
                       sort = T)

unlist(rest.nest_all)

null.com_all <- sapply(nulls.com_all, 
                   function(x) bipartite::nest.smdm(x = x,
                                                    constraints = Part_all, 
                                                    weighted = T, 
                                                    decreasing = "abund"))
WNODA.null.com_all <- unlist(null.com_all[3,])
WNODAsm.null.com_all <- unlist(null.com_all[8,])
WNODAdm.null.com_all <- unlist(null.com_all[9,])

#Plot the observed nestedness value against the distribution of randomized values
quartz()

par(mfrow = c(1,3))

plot(density(WNODA.null.com_all), xlim=c(min(obs.com_all[3], min(WNODA.null.com_all)),
                                     max(obs.com_all[3], max(WNODA.null.com_all))), 
     main="observed vs. randomized", xlab = "WNODA matrix")
abline(v=obs.com_all[3], col="red", lwd=2)

plot(density(WNODAsm.null.com_all), xlim=c(min(obs.com_all[8], min(WNODAsm.null.com_all)),
                                       max(obs.com_all[8], max(WNODAsm.null.com_all))), 
     main="observed vs. randomized", xlab = "WNODAsm matrix")
abline(v=obs.com_all[8], col="red", lwd=2)    

plot(density(WNODAdm.null.com_all), xlim=c(min(obs.com_all[9], min(WNODAdm.null.com_all)),
                                       max(obs.com_all[9], max(WNODAdm.null.com_all))), 
     main="observed vs. randomized", xlab = "WNODAdm matrix")
abline(v=obs.com_all[9], col="red", lwd=2)    

par(mfrow = c(1,1))

#Estimate the p-values

#Nestedness in th entire network
praw.WNODA_all <- sum(WNODA.null.com_all>obs.com_all[3]) / length(WNODA.null.com_all)
p.WNODA_all <- ifelse(praw.WNODA_all > 0.5, 1- praw.WNODA_all, praw.WNODA_all)    # P-value
p.WNODA_all

#Nestedness within the modules
praw.WNODAsm_all <- sum(WNODAsm.null.com_all>obs.com_all[8]) / length(WNODAsm.null.com_all)
p.WNODAsm_all <- ifelse(praw.WNODAsm_all > 0.5, 1- praw.WNODAsm_all, praw.WNODAsm_all)    # P-value
p.WNODAsm_all

#Nestedness between the modules
praw.WNODAdm_all <- sum(WNODAdm.null.com_all>obs.com_all[9]) / length(WNODAdm.null.com_all)
p.WNODAdm_all <- ifelse(praw.WNODAdm_all > 0.5, 1- praw.WNODAdm_all, praw.WNODAdm_all)    # P-value
p.WNODAdm_all


##### PLOT THE COMPOUND TOPOLOGY

par(mfrow = c(1,1))

#Sort the matrix in a way that facilitates visualizing the compound topology
data.comp_all <- bipartite::sortmatrix(matrix = data_all, topology = "compound", 
                                   sort_by = "weights", 
                                   row_partitions = row.Part_all, 
                                   col_partitions = col.Part_all)

#Assign colors for the modules
modcol_all <- rainbow((length(unique(Part_all))), alpha=1, s = 1, v = 1)

#Plot the matrix
png("figures/compound_all.png", width = 3000, height = 3000, res = 300)
plotmatrix(data.comp_all$matrix, 
           row_partitions = data.comp_all$row_partitions, 
           col_partitions = data.comp_all$col_partitions, 
           border = T,
           binary = F,
           modules_colors = modcol_all,
           within_color = modco_all, 
           between_color = "lightgrey")
dev.off()


##### EXPORT A SUMMARY OF THE TOPOLOGICAL RESULTS

sink(file = "results/results_topology_all.txt")


paste("The network has", nrow(data_all), "rows and", ncol(data_all), "columns.")
cat("\n")
paste("The network's modularity (DIRT_LPA+) is", round(Mod_all@likelihood, 2), ",", "P =", round(Mod.sig_all, 2), ",", "and it contains", length(unique(Part_all)),"modules.")
cat("\n")
paste("The network's nestedness (WNODF) is", round(Nest_all/100, 2),",", "P =", round(Nest.sig_all, 2))
cat("\n")
paste("The network shows the following scores of nestedness (WNODA):")
cat("\n")
paste("Entire network =", round(rest.nest_all$WNODAmatrix/100, 2), ",", "P =", round(p.WNODA_all, 2))
cat("\n")
paste("Between the modules =", round(rest.nest_all$WNODA_DM_matrix/100, 2), ",", "P =", round(p.WNODAdm_all, 2))
cat("\n")
paste("Within the modules =", round(rest.nest_all$WNODA_SM_matrix/100, 2), ",", "P =", round(p.WNODAsm_all, 2))
cat("\n")

sink(file = NULL, )


sink(file = "results/results_specialization.txt")


paste("The network's specialization (H2) is", round(Spec, 2),",", "P =", round(Spec.sig, 2))
cat("\n")


sink(file = NULL, )


################################################################################
##### SPECIES ROLES PEP
################################################################################
#Agora vamos calcular essas duas centralidades (c e z)#
papeis.Host.full_all <- czvalues(data.mod_all,weighted = F, level = "lower")
papeis.Parasite.full_all <- czvalues(data.mod_all, weighted = F,level = "higher")      

#Podemos combinar esses resultados, que foram calculados no bipartite e, assim, respeitam a 
#estrutura bipartida da rede, em ?nico data frame, inluindo a identidade de classe com um fator:
papeis.Host.full2_all <- as.data.frame(papeis.Host.full_all)
papeis.Parasite.full2_all <- as.data.frame(papeis.Parasite.full_all)
papeis.full_all <- rbind(papeis.Host.full2_all, papeis.Parasite.full2_all)
papeis.full_all$classe <- as.factor(c(rep("host", length(papeis.Host.full_all$c)),
                                  rep("parasite", length(papeis.Parasite.full_all$c))))

write.csv(papeis.full_all, "results/papeis_full_all.csv")
papeis.bip<-papeis.full_all


#Cor Species Roles##
papeis.bip2 <- papeis.bip

#i) kinless hubs ( [z ≥ 2.5 and c > 0.75])
papeis.bip2$papel <- ifelse(papeis.bip$c > 0.75 & papeis.bip$z >= 2.5, "darkblue", "#D1E5F0")

#ii) connector hubs ( [z ≥ 2.5 and 0.30 < c ≤ 0.75])
papeis.bip2$papel <- ifelse(papeis.bip$c > 0.30 & papeis.bip$c <= 0.75 & papeis.bip$z >= 2.5, "#F4A582", papeis.bip2$papel)

#iii) provincial hubs ( [z ≥ 2.5 and c ≤  0.30])
papeis.bip2$papel <- ifelse(papeis.bip$c <= 0.30 & papeis.bip$z >= 2.5, "#B2182B", papeis.bip2$papel) #B2182B

#iv) non-hub kinless vertices ( [z  < 2.5 and c > 0.80])
papeis.bip2$papel <- ifelse(papeis.bip$c > 0.80 & papeis.bip$z  < 2.5, "#92C5DE", papeis.bip2$papel)

#v) non-hub connector vertices ( [z  < 2.5 and 0.62 < c ≤ 0.80])
papeis.bip2$papel <- ifelse(papeis.bip$c > 0.62 & papeis.bip$c <= 0.80 & papeis.bip$z  < 2.5, "#4393C3", papeis.bip2$papel)

#vi)peripheral vertices ( [z  < 2.5 and 0.05 < c ≤ 0.62])
papeis.bip2$papel <- ifelse(papeis.bip$c > 0.05 & papeis.bip$c <= 0.62 & papeis.bip$z  < 2.5,"#FDDBC7" , papeis.bip2$papel) 

#vii)ultraperipheral vertices [z  < 2.5 and c ≤ 0.05])
papeis.bip2$papel <- ifelse(papeis.bip$c <= 0.05 & papeis.bip$z  < 2.5, "#D1E5F0", papeis.bip2$papel)

papeis.bip2$papel <- ifelse(is.na(papeis.bip2$pape), "#D1E5F0", papeis.bip2$papel)

#Figure 2
data3_all<-data2_all


tiff(filename= "figures/all_FIG2.tif", res= 300, height= 3000, width= 4100)
par(mfrow=c(1,1),mar=c(1,1,1,5))
plot(data3_all,
     vertex.color = papeis.bip2$papel, 
     vertex.size = 8,
     vertex.label = V(data2_all)$name,
     vertex.label.font= 2,
     vertex.label.color = "black",
     vertex.label.cex = 0.65,
     vertex.frame.color = "black",
     edge.color = adjustcolor("grey", alpha.f = .5),
     edge.curved = 0.3,
     edge.width = 3,
     layout=lay1_all 
)
legend(x = 1,y = -0.7, legend = c("Provincial Hubs","Peripheral Vertices","Ultraperipheral Vertices"),
       pch = c(19),
       col = c("#B2182B","#FDDBC7","#D1E5F0"),  
       title = "Species Roles",
       text.col = "gray20", title.col = "black",
       box.lwd = 0, cex = 1)
par(mfrow=c(1,1))
dev.off()

#Figure3
tiff(filename= "figures/all_FIG3.tif", res= 300, height= 3000, width= 5100)
par(mfrow=c(1,1),mar=c(5,5,5,5))
plot(data = papeis.bip2,
     c ~ z,
     pch = 20,
     cex = 4,
     col = papeis.bip2$papel,cex.axis=1.2,cex.names=1.2,ylab="Pi", xlab="zi")
text(papeis.bip2$z,                        
     papeis.bip2$c, 
     labels = rownames(papeis.bip2), 
     cex = 0.8, col = "black", pos = 4, srt=60)

legend(x = 1,y = -0.7, legend = c("Provincial Hubs","Peripheral Vertices","Ultraperipheral Vertices"),
       pch = c(19),
       col = c("#B2182B","#FDDBC7","#D1E5F0"),  
       title = "Species Roles",
       text.col = "gray20", title.col = "black",
       box.lwd = 0, cex = 1)
par(mfrow=c(1,1))
dev.off()
