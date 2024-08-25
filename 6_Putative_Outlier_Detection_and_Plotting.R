#####Outlier Detection######
#Following this tutorial https://bcm-uga.github.io/pcadapt/articles/pcadapt.html
#Set working directory, read in data
setwd("~/Gentianella_publication/final_run/dataset3/")
library(vcfR)
library(adegenet)
library(dartR)

#colours
calcis_colours <- c("#EAB64D","#AFC494","#0BB24E","#5c6f9d","deepskyblue","#DA8051","#857E61"
                    ,"#90C8AD","orchid4","yellowgreen","#C97E8C","#CBD4C2")

##Read the files
ecology <- read.vcfR("populations.snps.vcf")

#Population data
pop.data <- read.table("popmap_lessthan50percentofmean.txt", sep = "\t", header = TRUE)

v_genlight <- vcfR2genlight(ecology)
pop(v_genlight) <- pop.data$Population


##Using PCAdapt####
library(dartR)
library(pcadapt)

#Set the path to the VCF file
path <- "~/Gentianella_publication/final_run/dataset3/populations.snps.vcf"

#Read in the file
v_adapt <- read.pcadapt(path, type = ("vcf"))

##Perform the pcadapt - specify a number of axes to retain initially
pc_v <- pcadapt(input=v_adapt,K=20)

#Check the axes using a screeplot. We are looking for where the points deviate from the straight
#Line (i.e., Cattell's rule)
plot(pc_v, option ="screeplot")
#Based on this it looks like K=10 is most appropriate

#Repeat with K=10
pc_v <- pcadapt(input=v_adapt,K=10)
plot(pc_v,option="manhattan")
plot(pc_v, option = "qqplot")

#Plot the p values, those with really small ones are outliers
hist(pc_v$pvalues, xlab = "p-values", main = "Frequency of p-values", breaks = 50, col = "orchid3")

#Add the P values to their own dataframe, rename the column
pval <- as.data.frame(pc_v$pvalues)
colnames(pval)[1] <- "p_value"
#Add the SNP number to the next column (assumes the pvalues kept the order of the SNPS)
pval$SNP <- v_genlight$loc.names
#Perform bonferroni method which is a conservative way of finding SNP outliers
pval$bf <- p.adjust(pval$p_value, method="bonferroni")
#add 1s to identify which SNPs are outliers and which aren't (0)
pval$outlier <- ifelse(pval$bf < 0.1, '1',
                       ifelse(pval$p_value > 0.1, '0', '0'))
#Check how many there are
length(which(pval$outlier==1))

#Separate the Two dataframes
SNP_ouliers <- pval[which(pval$outlier==1),]
SNP_normals <- pval[which(pval$outlier==0),]

#Create a new GL using dartR that contains only the outlier loci
gl_outliers <- gl.keep.loc(v_genlight, SNP_ouliers$SNP)


# Create two PCA for the outlier dataset 
o_PCA <- glPca(gl_outliers, center = TRUE, scale = TRUE, nf = 8, loadings = TRUE,
               alleleAsUnit = FALSE, useC = TRUE, parallel = FALSE,
               n.cores = NULL, returnDotProd=FALSE, matDotProd=NULL)

#Variation explained by axes
eig <- round((100*o_PCA$eig/sum(o_PCA$eig)), digits = 1)

#Extract PCA scores and merge with the population data
pca_scores <- as.data.frame(o_PCA$scores)
pca_scores["AccessID"] <- rownames(pca_scores)
row.names(pca_scores) <- NULL
pca_scores <- merge(pca_scores, pop.data, by ="AccessID")

#PLOT - the sprintf function allows us to use a variable in a string 
outlierPCA1_2 <- ggplot(pca_scores, aes(x=PC1, y = PC2)) + 
  geom_hline(yintercept = 0, color = "grey", linewidth=0.5, linetype = 2) + 
  geom_vline(xintercept = 0, color = "grey", linewidth =0.5, linetype = 2) +
  geom_point(aes(color=Population, shape=Taxon), size =2) +
  stat_ellipse(aes(group = Population, color = Population)) +
  xlab(sprintf("PC1 %s %%", eig[1])) +
  ylab(sprintf("PC2 %s %%", eig[2])) +
  ggtitle("PCA of outliers axes 1 & 2") +
  scale_colour_manual(values=calcis_colours) +
  theme_classic() + theme(legend.position = "none")
  
  

outlierPCA3_4 <- ggplot(pca_scores, aes(x=PC3, y = PC4)) + 
  geom_hline(yintercept = 0, color = "grey", linewidth=0.5, linetype = 2) + 
  geom_vline(xintercept = 0, color = "grey", linewidth =0.5, linetype = 2) +
  geom_point(aes(color=Population, shape=Taxon), size =2) +
  stat_ellipse(aes(group = Population, color = Population)) +
  xlab(sprintf("PC3 %s %%", eig[3])) +
  ylab(sprintf("PC4 %s %%", eig[4])) +
  ggtitle("PCA of outliers axes 3 & 4") +
  scale_colour_manual(values=calcis_colours) +
  theme_classic() + theme(legend.position = "right") + scale_shape(guide = FALSE) 
  
  

##Repeat this but with neutral loci
gl_normals <- gl.keep.loc(v_genlight, SNP_normals$SNP)
##Now run it back through the PCA code 
n_PCA <- glPca(gl_normals, center = TRUE, scale = FALSE, nf = 5, loadings = TRUE,
               alleleAsUnit = FALSE, useC = TRUE, parallel = FALSE,
               n.cores = NULL, returnDotProd=FALSE, matDotProd=NULL)


#Variation explained by axes
eig <- round((100*n_PCA$eig/sum(n_PCA$eig)), digits = 1)

#Extract PCA scores and merge with the population data
pca_scores <- as.data.frame(n_PCA$scores)
pca_scores["AccessID"] <- rownames(pca_scores)
row.names(pca_scores) <- NULL
pca_scores <- merge(pca_scores, pop.data, by ="AccessID")

#PLOT - the sprintf function allows us to use a variable in a string 
normalPCA1_2 <- ggplot(pca_scores, aes(x=PC1, y = PC2)) + 
  geom_hline(yintercept = 0, color = "grey", linewidth=0.5, linetype = 2) + 
  geom_vline(xintercept = 0, color = "grey", linewidth =0.5, linetype = 2) +
  geom_point(aes(color=Population, shape = Taxon), size =2) +
  stat_ellipse(aes(group = Population, color = Population)) +
  xlab(sprintf("PC1 %s %%", eig[1])) +
  ylab(sprintf("PC2 %s %%", eig[2])) +
  ggtitle("PCA of neutral SNPS axes 1 & 2") +
  scale_colour_manual(values = calcis_colours) +
  theme_classic() +
  theme(legend.position = "none") 
  


normalPCA3_4 <- ggplot(pca_scores, aes(x=PC3, y = PC4)) + 
  geom_hline(yintercept = 0, color = "grey", linewidth=0.5, linetype = 2) + 
  geom_vline(xintercept = 0, color = "grey", linewidth =0.5, linetype = 2) +
  geom_point(aes(color=Population, shape = Taxon), size =2) +
  stat_ellipse(aes(group = Population, color = Population)) +
  xlab(sprintf("PC3 %s %%", eig[3])) +
  ylab(sprintf("PC4 %s %%", eig[4])) +
  ggtitle("PCA of neutral SNPS axes 3 & 4") +
  scale_colour_manual(values = calcis_colours, guide = FALSE) +
  theme_classic() +
  theme(legend.position = "right")
  

combined <- outlierPCA1_2+outlierPCA3_4+normalPCA1_2+normalPCA3_4 

combined 

## Investigate potential for overdominance
outlier_matrix <- as.matrix(gl_outliers) 
normal_matrix <- as.matrix(gl_normals)

# Recode genotypes and add pop data
outlier_fre <- data.frame(row.names = gl_outliers$ind.names)
outlier_fre$ref <- rowSums(outlier_matrix == '0', na.rm = TRUE)
outlier_fre$alt <- rowSums(outlier_matrix == '2', na.rm = TRUE)
outlier_fre$het <- rowSums(outlier_matrix == '1', na.rm = TRUE)
outlier_fre$pops <- pop.data$Population

outlier_avgs <- aggregate(ref ~ pops, data = outlier_fre, mean)
outlier_avgs$alt <- (aggregate(alt ~ pops, data = outlier_fre, mean))$alt
outlier_avgs$het <- (aggregate(het ~ pops, data = outlier_fre, mean))$het

outlier_avgs <- t(outlier_avgs)
colnames(outlier_avgs) <- outlier_avgs[1,]
outlier_avgs <- outlier_avgs[-1,]

# Plot
par(mar = c(8, 4.1, 4.1, 7))
barplot(outlier_avgs, las = 2, cex.names = 0.75,
        col = c("orchid2","orchid3", "orchid4"),
        main = "Average Allele Frequencies in Outlier Dataset",
        ylab = "Proportion of Alleles",
        axes = TRUE,
        legend.text = rownames(outlier_avgs),
        args.legend = list(x = "topright",
                           inset = c(-0.2, 0)))

# Recode genotypes 
normal_fre <- data.frame(row.names = gl_outliers$ind.names)
normal_fre$ref <- rowSums(normal_matrix == '0', na.rm = TRUE)
normal_fre$alt <- rowSums(normal_matrix == '2', na.rm = TRUE)
normal_fre$het <- rowSums(normal_matrix == '1', na.rm = TRUE)
normal_fre$pops <- pop.data$Population

normal_avgs <- aggregate(ref ~ pops, data = normal_fre, mean)
normal_avgs$alt <- (aggregate(alt ~ pops, data = normal_fre, mean))$alt
normal_avgs$het <- (aggregate(het ~ pops, data = normal_fre, mean))$het

normal_avgs <- t(normal_avgs)
colnames(normal_avgs) <- normal_avgs[1,]
normal_avgs <- normal_avgs[-1,]

# Plot
par(mar = c(8, 4.1, 4.1, 7))
barplot(normal_avgs, las = 2, cex.names = 0.75,
        col = c("orchid2","orchid3", "orchid4"),
        main = "Average Allele Frequencies in Neutral Dataset",
        ylab = "Proportion of Alleles",
        axes = TRUE,
        legend.text = rownames(outlier_avgs),
        args.legend = list(x = "topright",
                           inset = c(-0.2, 0)))

# Investigate Fst of the twi datasets
library(StAMPP)
set.seed(28529)

pvaluefst_outliers <- stamppFst(gl_outliers, nboots = 50000, percent = 95.000, nclusters = 3)
#If you use too few bootstraps, your p values will be small (i.e., 2dp)
fstmat_outliers <- pvaluefst_outliers$Fsts

#even with over 100,000 boostraps we didn't get P values of greater than 0
# based on trials using less differentiated groups, we can conclude everything has a P value of at least ~0.00001
pvals_outliers <- pvaluefst_outliers$Pvalues
pvals_outliers <- format((pvals_outliers[1-12,1-12]+0.00001), scientific = FALSE)

#apply BF correction for multiple comparisons
adjusted_outliers <- p.adjust(pvals_outliers, method = "bonferroni")
adjusted_outliers

pvaluefst_neutrals <- stamppFst(gl_normals, nboots = 50000, percent = 95.000, nclusters = 3)
#If you use too few bootstraps, your p values will be small (i.e., 2dp)
fstmat_neutrals <- pvaluefst_neutrals$Fsts

#even with over 100,000 boostraps we didn't get P values of greater than 0
# based on trials using less differentiated groups, we can conclude everything has a P value of at least ~0.00001
pvals_neutrals <- pvaluefst_neutrals$Pvalues
pvals_neutrals <- format((pvals_neutrals[1-12,1-12]+0.00001), scientific = FALSE)

#apply BF correction for multiple comparisons
adjusted_neutrals <- p.adjust(pvals_neutrals, method = "bonferroni")
adjusted_neutrals

write.table(fstmat_neutrals, file = "Pairwise_Fst_neutrals.txt", quote = FALSE, sep = "\t")
write.table(fstmat_outliers, file = "Pairwise_Fst_outliers.txt", quote = FALSE, sep = "\t")


