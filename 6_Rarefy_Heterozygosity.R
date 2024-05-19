##Calculate Observed, Expected heterozygosity, 
setwd("~/Gentianella_publication/final_run/")
library(vcfR)
library(adegenet)
library(dartR)
library(dbplyr)

calcis <- read.vcfR("noreps_maf_maxhet_noparas_randomSNP.vcf")

calcis_genlight <- vcfR2genlight(calcis)

popdata <- read.table("popmap_noreps_samplesize.txt", sep = "\t", header = TRUE)
pop(calcis_genlight) <- popdata$Population

###### Function #####
#Genlight: A genlight object
#Popdata: Tab-separated file with two columns (AccessID and Population) which contain the sample name and populations
#num: number of times to re-calculate for each population
#resamp: size of population you want to resample to
SampleSizeDiversity2 <- function(Genlight, pop.data, num, resamp){
  SNPS=Genlight
  pop.data=pop.data
  num=num
  resamp=resamp
  #Compliance check and make sure there are enough individuals for each population
  SNPS <- gl.compliance.check(SNPS)
  
  #Create empty dataframe to store output
  test_stats <- data.frame(matrix(ncol=38, nrow=0))
  x <- c("pop", "n.Ind", "n.Loc", "n.Loc.adj", "polyLoc", "monoLoc", "all_NALoc", "Ho", "HoSD", "HoSE",     
          "HoLCI", "HoHCI", "Ho.adj", "Ho.adjSD", "Ho.adjSE", "Ho.adjLCI", "Ho.adjHCI", "He", "HeSD", "HeSE",     
          "HeLCI", "HeHCI", "uHe", "uHeSD", "uHeSE", "uHeLCI", "uHeHCI", "He.adj", "He.adjSD", "He.adjSE", 
          "He.adjLCI", "He.adjHCI", "FIS", "FISSD", "FISSE", "FISLCI", "FISHCI", "nInd" ) 
  colnames(test_stats) <- x
  samplesizes <- popdata %>% count(popdata$Population, sort = TRUE) #Count the number of samples per pop
  colnames(samplesizes) <- c("Population", "n")
  
  #Omit any populations which have sample sizes equal to or less than the amount chosen to subset
  if ((sum(samplesizes$n <= resamp))>=1)  {
    smallpops <- subset(samplesizes, n <=resamp)
    bigpops <- subset(samplesizes, n > resamp)
    smallpopsremoved <- gl.drop.pop(SNPS, pop.list = c(smallpops$Population), 
                                    mono.rm = TRUE, recalc = TRUE, verbose = 0)
    smallpopsremoved$other$loc.metrics <- as.data.frame(smallpopsremoved$other$loc.metrics)
    #mono.rm is set to false here because we want to simulate some of the loci being polymorphic in other populations, therefore being monomorphic in our target population
    
    #Make a list of populations
    pops <- bigpops$Population
    #For each population calculate genetic diversity based on the desired subset
    for (pop in pops) {
      print(pop)
      #Create list of samples in that pop
      indlist <- pop.data$AccessID[pop.data$Population==pop]
      counter <- 1
      #Downsample from population
      for(i in 1:num){
      
        #Create a random list of samples to remove
        print(counter)
        drop <- sample(indlist, (length(indlist)-resamp), replace = FALSE)
        
        #List which samples calculation was based on
        subsetlist <- indlist[!(indlist %in% drop)]
        
        #Reduce the dataset down to the desired number of samples
        subset_genlight <-gl.drop.ind(smallpopsremoved, ind.list = c(drop), recalc = TRUE, mono.rm = TRUE, verbose = 0 )
        subset_genlight$other$loc.metrics <- as.data.frame(subset_genlight$other$loc.metrics)
        
        #Calculate heterozygosity
        stats <- gl.report.heterozygosity(subset_genlight, plot.display = FALSE, verbose = 0)
        statsforpop <- stats[stats$pop == pop,]
        statsforpop$nInd <- list(toString(subsetlist))
        test_stats[nrow(test_stats) +1,] <- statsforpop
        counter <- (counter+1)
      }
    }
    
  } else { #if none of the populations are the same size or smaller than the subset, proceed normally here
    pops <- unique(pop.data$Population)
    for (pop in pops) {
      print(pop)
      #Create list of samples in that pop
      indlist <- pop.data$AccessID[pop.data$Population==pop]
      counter <- 1
      #Downsample from population
      for(i in 1:num){
        
        #Create a random list of samples to remove
        print(counter)
        drop <- sample(indlist, (length(indlist)-resamp), replace = FALSE)
        
        #List which samples calculation was based on
        subsetlist <- indlist[!(indlist %in% drop)]
        
        #Reduce the dataset down to the desired number of samples
        subset_genlight <-gl.drop.ind(smallpopsremoved, ind.list = c(drop), recalc = TRUE, mono.rm = TRUE, verbose = 0 )
        subset_genlight$other$loc.metrics <- as.data.frame(subset_genlight$other$loc.metrics)
        
        #Calculate heterozygosity
        stats <- gl.report.heterozygosity(subset_genlight, plot.display = FALSE, verbose = 0)
        statsforpop <- stats[stats$pop == pop,]
        statsforpop$nInd <- list(toString(subsetlist))
        test_stats[nrow(test_stats) +1,] <- statsforpop
        counter <- (counter+1)
      }
    }
  }
    
  return(test_stats)
}

#Calculate
set.seed(1837)
Trial <- SampleSizeDiversity2(calcis_genlight, popdata, 1000, 3)
boxplot(Trial$Ho~Trial$pop)

#Calculate mean and SD for each population and append to DF
resampledhet <- data.frame(matrix(0, ncol = 12, nrow = 2))
for (pop in (unique(Trial$pop))) {
  print(pop)
  resampledhet[[pop]] <- c(mean(Trial$Ho[Trial$pop == pop]), sd(Trial$Ho[Trial$pop == pop]))
}

#Write your files to txt 
T2 <- apply(Trial,2,as.character)
write.table(T2, file = "rarefiedgendiversity.txt", sep = "\t", row.names = F, col.names = T, quote = FALSE)
write.table(resampledhet, file = "rarefiedgeneticdiversitysummary.txt", sep = "\t", row.names = F, col.names = T, quote = FALSE)
