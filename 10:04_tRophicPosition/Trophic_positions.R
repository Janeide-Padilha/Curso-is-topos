######################################################################################################################
#TROPHIC POSITION ESTIMATION
######################################################################################################################

#Most of this script is based on Claudio Quezada-Romegialli's example scripts for his package tRophicPosition. 
#If you are looking for more example scripts and guidance to use the package, 
#please visit https://github.com/clquezada/tRophicPosition and 
#https://cran.r-project.org/web/packages/tRophicPosition/vignettes/.

#Prerequisites: Install the current version of R: https://cran.r-project.org/, 
#of RStudio: https://www.rstudio.com/products/rstudio/download/, 
#and of JAGS: http://mcmc-jags.sourceforge.net/.
#Create a Rstudio project associated with your working directory, and make sure your work will be saved automatically.

#Now let's get started!
install.packages("devtools")
devtools::install_github("https://github.com/clquezada/tRophicPosition")
getwd()
setwd("")
######################################################################################################################
#PART 1 - WORKING ON A SINGLE SPECIES
######################################################################################################################

#First, we will install the tRophicPosition package, and load it into our session environment.
install.packages("tRophicPosition", dependencies = TRUE)

library(tRophicPosition)

#Import the data. This file contains isotopic ratios for the consumers and the baseline items that we will use into
#our examples (5 consumer species, 2 baselines).
Fulldata <- read.csv("Data_tRophicPosition.csv", header=TRUE)

#We will also need TEF values. Here, we will use values of the McCutchan et al. paper, that can be extracted from the
#package directly. Alternatively, if we had more suitable TEFs for the species we are working with, we could import
#them from a file.
TEFs <- TDF(author = "McCutchan", element = "both", type = "all")

#To start with this first example, we will focus on only one species, the pycnogonid Decolopoda australis. To do so,
#we will have to select our species in the dataset, and format the data about Decolopoda so that they can be used 
#properly by the package's functions. This is done by mapping our data to relevant corresponding fields.
Decolopoda <- loadIsotopeData(Fulldata,
                              group = "Decolopoda",
                              groupsColumn = "Species",
                              consumer = "Decolopoda", 
                              consumersColumn = "Item",
                              b1 = "Pelagic_BL", 
                              b2 = "Sympagic_BL",
                              baselineColumn = "Item",
                              deltaN = TEFs$deltaN,
                              deltaC = TEFs$deltaC)

#Now that we created this "IsotopeData" class object, tRophicPosition can easily produce a neat plot to check that
#everything is in order.
plot(Decolopoda, b1 = "Pelagic baseline", b2 = "Sympagic baseline")

#First, we can use tRophicPosition to calculate the trophic level of Decolopoda, using a formula partly based on the
#classic paper by Post (2002). Pay attention to the "lambda" parameter. Here, it is set to 1 because we will estimate
#our baselines using isotopic ratios of primary producers directly. If we were using primary consumers instead, lambda
#would be set to 2.
Param.TP.Decolopoda <- parametricTP(Decolopoda, lambda=1, print=TRUE)

#Now, let's use a Bayesian model to estimate trophic position. First we will set a few parameters. Here, we will build
#a model using two baselines, as well as trophic enrichment factors for both isotopes. Lambda is still set to 1. If we 
#wanted to add a prior (e.g. based on gut contents), we would also specify it here.
model.params <- jagsBayesianModel(model = "twoBaselinesFull", lambda=1)


#Time to run our model. Many parameters can be changed, but this goes beyond the scope of this introduction. Those 
#"default" values are a trade-off between short run time and reliable results.
model.decolopoda <- TPmodel(data = Decolopoda, model.string = model.params,
                                 n.adapt = 20000, n.chains = 2)
posterior.decolopoda <- posteriorTP(model = model.decolopoda, n.iter = 20000,
                                 variable.names = "TP")

#Now, we can have a look at the model's output.
summary(posterior.decolopoda)

#As well as at the mode of the posterior distribution of solutions (value most commonly found in model solutions)
posterior.decolopoda.combined <- coda::mcmc(do.call(rbind, posterior.decolopoda))
getPosteriorMode(posterior.decolopoda.combined)

#Let's plot our model's results using a boxplot of credibility intervals. On this plot, by default, the dark, median, 
#and light grey are the 50, 75 and 95% credibility intervals, and the black dot is the mode.
plotTP(as.data.frame(posterior.decolopoda.combined), 
       xlab = NULL, xticklabels="Decolopoda", 
       ylims=c(2.2, 3.2),mgp = c(3, 0.7, 0))

#Finally, we can add the calculated value of trophic position to this plot, just to check whether the two methods
#yield similar results.
points(Param.TP.Decolopoda[[4]], col = "red", pch = 16)

#Now, let's move to a more complex example, with five different species.


######################################################################################################################
#PART 2 - WORKING ON MULTIPLE SPECIES
######################################################################################################################

#We will use the same data object as previously. However, this time we will not select a consumer but extract data for
#all 5 species.
Allspecies <- extractIsotopeData(Fulldata,
                                 b1 = "Pelagic_BL", 
                                 b2 = "Sympagic_BL",
                                 baselineColumn = "Item",
                                 consumersColumn = "Item",
                                 groupsColumn = "Species",
                                 deltaN = TEFs$deltaN,
                                 deltaC = TEFs$deltaC,
                                 d13C = "d13C", 
                                 d15N = "d15N")
                                 
#We can now look at data summaries, and plot the data for each species. We could do that in 5 times. 
#Alternatively, we can use a recursive function ("loop").
for (Species in Allspecies) {
  print(summary(Species))
  plot(Species)
} 

#Time to run our model. Since we have multiple species, this can take a while.
model.allspecies <- multiSpeciesTP(Allspecies,
                                   model = "twoBaselinesFull",
                                   lambda = 1,
                                   n.adapt = 20000, n.iter = 20000,
                                   burnin = 2000, n.chains = 2, print = FALSE)

#If we want to produce boxplots like before, we have to extract model outputs for each species
Decolopoda.TP <- model.allspecies$TPs$Decolopoda.Decolopoda.2bf
Flabegraviera.TP <- model.allspecies$TPs$Flabegraviera.Flabegraviera.2bf
Odontaster.TP <- model.allspecies$TPs$Odontaster.Odontaster.2bf
Parborlasia.TP <- model.allspecies$TPs$Parborlasia.Parborlasia.2bf
Sterechinus.TP <- model.allspecies$TPs$Sterechinus.Sterechinus.2bf

#The combine them all in a single data frame
AllSpecies.TP <- as.data.frame(cbind(Decolopoda.TP,Flabegraviera.TP,Odontaster.TP,Parborlasia.TP,Sterechinus.TP))

#That we can plot.
plotTP(AllSpecies.TP, xlab = "Species", 
       xticklabels=c("Decolopoda","Flabegraviera","Odontaster","Parborlasia","Sterechinus"), 
       ylims=c(1, 3.5),
       mgp = c(3, 0.7, 0))

#If we wish to add the parametric TP estimates as red dots, we need to calculate them for each species.
#Note that we will not do it for Decolopoda: it was already done in the first part of the script.
Flabegraviera <- loadIsotopeData(Fulldata,
                              group = "Flabegraviera",
                              groupsColumn = "Species", 
                              consumer = "Flabegraviera", 
                              consumersColumn = "Item",
                              b1 = "Pelagic_BL", 
                              b2 = "Sympagic_BL",
                              baselineColumn = "Item",
                              deltaN = TEFs$deltaN,
                              deltaC = TEFs$deltaC)
Param.TP.Flabegraviera <- parametricTP(Flabegraviera, lambda=1, print=TRUE)

Odontaster <- loadIsotopeData(Fulldata,
                              group = "Odontaster",
                              groupsColumn = "Species", 
                              consumer = "Odontaster", 
                              consumersColumn = "Item",
                              b1 = "Pelagic_BL", 
                              b2 = "Sympagic_BL",
                              baselineColumn = "Item",
                              deltaN = TEFs$deltaN,
                              deltaC = TEFs$deltaC)
Param.TP.Odontaster <- parametricTP(Odontaster, lambda=1, print=TRUE)

Parborlasia <- loadIsotopeData(Fulldata,
                              group = "Parborlasia",
                              groupsColumn = "Species", 
                              consumer = "Parborlasia", 
                              consumersColumn = "Item",
                              b1 = "Pelagic_BL", 
                              b2 = "Sympagic_BL",
                              baselineColumn = "Item",
                              deltaN = TEFs$deltaN,
                              deltaC = TEFs$deltaC)
Param.TP.Parborlasia <- parametricTP(Parborlasia, lambda=1, print=TRUE)

Sterechinus <- loadIsotopeData(Fulldata,
                              group = "Sterechinus",
                              groupsColumn = "Species",  
                              consumer = "Sterechinus", 
                              consumersColumn = "Item",
                              b1 = "Pelagic_BL", 
                              b2 = "Sympagic_BL",
                              baselineColumn = "Item",
                              deltaN = TEFs$deltaN,
                              deltaC = TEFs$deltaC)
Param.TP.Sterechinus <- parametricTP(Sterechinus, lambda=1, print=TRUE)

#Now we can combine them all in a frame
Param.TPs <-cbind(Param.TP.Decolopoda[[4]],Param.TP.Flabegraviera[[4]],Param.TP.Odontaster[[4]],
                  Param.TP.Parborlasia[[4]],Param.TP.Sterechinus[[4]])

#And add them to the plot
points(1:5,Param.TPs, col = "red", pch = 16)

#Finally, we can use our model to compute probabilities that species have different trophic positions. This command
#will return a double entry table, where each values answer the question "What is the probability that consumer in 
#row Y has a trophic position inferior or equal to the consumer in column X?". Those probabilities can be considered 
#meaningful when superior to 0.95.
pairwiseTP <- pairwiseComparisons(model.allspecies$TPs, print = TRUE)

#That will be it for today! Any questions?


######################################################################################################################
#PART 3 - INDEPENDENT WORK
######################################################################################################################

#Now it's your turn to work! Using what we learned today, try to complete the following tasks / answer the following 
#questions.


#1) Create a full model (two baselines, trophic fractionation for both isotopes) for Odontaster validus, and represent
#its output as a boxplot.

#2) Create a one baseline model based on the pelagic baseline for Odontaster validus, and represent its output 
#as a boxplot.

#3) Create another one baseline model for Odontaster validus, but based on sympagic algae. Represent its output 
#as a boxplot.
#Tip: you can interchange baselines 1 and 2 in the loadIsotopeData() function

#4) Plot the output of the three 3 models on a single graph.

#5) Interpret the result. What can you tell from this graph?


######################################################################################################################
#END OF SCRIPT
######################################################################################################################