#Setwd, load libraries
setwd("~/Gentianella_publication/final_run/")

#Libraries
library(ggplot2)
library(vcfR)
library(adegenet)
library(StAMPP)
library(dartR)

#Colour palette
calcis_colours <- c("#EAB64D","#AFC494","#0BB24E","#5c6f9d","deepskyblue","#DA8051","#857E61"
                    ,"#90C8AD","orchid4","yellowgreen","#C97E8C","#CBD4C2")
#read in VCF file
calcis <- read.vcfR("noreps_R0.8_minmaf_0.05_max_obshet0.75_randomSNP_noparas.vcf")

#read in pop data - this will allow us to see if groups match up with a prior expectations
pop.data <- read.table("popmap_noreps.txt", sep = "\t", header = T)


#Create Genlight from the vcfR object and add information to it
c_genlight <- vcfR2genlight(calcis)
pop(c_genlight) <- pop.data$Population

####First Analysis: PCA###
#PCA using adgenet and ggplot2
c_PCA <- glPca(c_genlight, center = TRUE, scale = TRUE, nf = 5, loadings = TRUE,
               alleleAsUnit = FALSE, useC = TRUE, parallel = FALSE,
               n.cores = NULL, returnDotProd=FALSE, matDotProd=NULL)

barplot(100*c_PCA$eig/sum(c_PCA$eig), col = heat.colors(90), main="PCA Eigenvalues")
title(ylab="Percent of variance\nexplained", line = 2)
title(xlab="Eigenvalues", line = 1)

#Calculate variation explained by axes
eig <- round((100*c_PCA$eig/sum(c_PCA$eig)), digits = 1)

#Extract PCA scores and merge with the population data
pca_scores <- as.data.frame(c_PCA$scores)
pca_scores["AccessID"] <- rownames(pca_scores)
row.names(pca_scores) <- NULL
pca_scores <- merge(pca_scores, pop.data, by ="AccessID")

#PLOT the PCA - the sprintf function allows us to use a variable in a string 

PCA_axes1_2 <- ggplot(pca_scores, aes(x=PC1, y = PC2)) + 
  geom_hline(yintercept = 0, color = "grey", linewidth=0.5, linetype = 2) + 
  geom_vline(xintercept = 0, color = "grey", linewidth =0.5, linetype = 2) +
  geom_point(aes(color=Population, shape = Taxon), size =2) +
  stat_ellipse(aes(group = Population, color = Population)) +
  xlab(sprintf("PC1 %s %%", eig[1])) +
  ylab(sprintf("PC2 %s %%", eig[2])) +
  ggtitle("PCA axes 1 & 2") +
  scale_colour_manual(values = calcis_colours) +
  theme_classic() +
  theme(legend.position = "none") 


PCA_axes3_4 <- ggplot(pca_scores, aes(x=PC3, y = PC4)) + 
  geom_hline(yintercept = 0, color = "grey", linewidth=0.5, linetype = 2) + 
  geom_vline(xintercept = 0, color = "grey", linewidth =0.5, linetype = 2) +
  geom_point(aes(color=Population, shape = Taxon), size =2) +
  stat_ellipse(aes(group = Population, color = Population)) +
  xlab(sprintf("PC3 %s %%", eig[3])) +
  ylab(sprintf("PC4 %s %%", eig[4])) +
  ggtitle("PCA axes 3 & 4") +
  scale_colour_manual(values = calcis_colours) +
  theme_classic() +
  theme(legend.position = "right") +
  theme(legend.text=element_text(size=9))
  

PCA_axes1_2 + PCA_axes3_4

###Write a .phylip file for splitstree###

#Create a matrix of pairwise Nei's distances between samples
gen_dist <- stamppNeisD(c_genlight, pop = FALSE)

#Write this out as a .phy file to the working directory
stamppPhylip(gen_dist, file = "input_for_splitstree")

###AMOVA###
#Libraries
library(ade4)
library(poppr)

#Add the AMOVA strata file
strata(c_genlight) <- pop.data[2]
nameStrata(c_genlight) <- ~Population

#Perform an AMOVA that tests at the population level
amova_obj <- poppr.amova(c_genlight, 
                         hier = ~Population,
                         threads = 3,
                         method = "ade4",
                         within = FALSE)

#Check your results
amova_obj
amova_obj$varcomp/sum(amova_obj$varcomp)

#Test if pops significantly different using randomization test
set.seed(2847)
c_sig <- randtest(amova_obj, nrepet = 10000)


#Visualise the results - black line shouldn't fall within observed 
plot(c_sig)
c_sig

write.table(amova_obj$componentsofcovariance, file = "AMOVA_covariance.txt", quote = FALSE, sep = "\t")

###Pairwise Fst###
#Fst with P values from the StAMMP package
pvaluefst <- stamppFst(c_genlight, nboots = 50000, percent = 95.000, nclusters = 3)
#If you use too few bootstraps, your p values will be small (i.e., 2dp)
fstmat <- pvaluefst$Fsts

#even with over 100,000 boostraps we didn't get P values of greater than 0
# based on trials using less differentiated groups, we can conclude everything has a P value of at least ~0.00001
pvals <- pvaluefst$Pvalues
pvals <- format((pvals[1-12,1-12]+0.00001), scientific = FALSE)

#apply BF correction for multiple comparisons
adjusted <- p.adjust(pvals, method = "bonferroni")
adjusted

write.table(fstmat, file = "Pairwise_Fst.txt", quote = FALSE, sep = "\t")


###IBD###
#Input a new popmap consisting of the sample names, population, taxon, and sampling coords (lat lon format)

popmap_samplingcoords <- read.table("popmap_and_sampcoords.txt", sep = "\t", header = T)

#Separate out only the lat lon 
latlon <- popmap_samplingcoords[,4:5]

#Perform IBD test using euclidean distances
ibd_result <- gl.ibd(c_genlight,
                     distance = 'euclidean',
                     coordinates = latlon,
                     permutations = 99999)

