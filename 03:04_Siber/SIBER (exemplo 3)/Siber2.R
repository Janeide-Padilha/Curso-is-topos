######################################################################################################################
#ISOTOPIC NICHE METRICS - CONVEX HULLS AND STANDARD ELLIPSES
######################################################################################################################

#Most of this script is based on Andrew Jackson's example scripts for his package SIBER. 
#More info at http://www.tcd.ie/Zoology/research/groups/jackson/
#If you are looking for more example scripts and guidance to use the package, 
#please visit https://github.com/AndrewLJackson/SIBER and https://cran.r-project.org/web/packages/SIBER/vignettes/.

#Prerequisites: Install the current version of R: https://cran.r-project.org/, 
#of RStudio: https://www.rstudio.com/products/rstudio/download/, 
#and of JAGS: http://mcmc-jags.sourceforge.net/ or https://sourceforge.net/projects/mcmc-jags/


#Now let's get started!
current_dir <- getwd()
print(current_dir)
setwd("/Users/janeidepadilha/Desktop/Multi-Crash/Aulas/SIA/Prática SIA/03_SIBER/")
######################################################################################################################
#PART 1 - CONVEX HULL METRICS ON COMMUNITIES
######################################################################################################################

#Let's install the SIBER, dplyr and ggplot packages from the CRAN repository, and load them into our session's environment.
install.packages("SIBER", dependencies = TRUE)
install.packages("dplyr", dependencies = TRUE)
install.packages("ggplot2", dependencies = TRUE)
library(SIBER)
library(dplyr)
library(ggplot2)

#Let's import the data, that are stored in the provided CSV file.
fulldata <- read.csv("SIBER_full.csv", header=T)

#Create a key for the community and group codes: it will be useful for quick reference later
key <- read.csv("Key.csv", header=F)

#Now, let's create a SIBER object that can be manipulated by the package.
siber.full <- createSiberObject(fulldata)

#Use SIBER to take a look at the data (colours = groups, symbols = communities)
par(mar=c(5.1, 5.1, 4.1, 2.1))
plotSiberObject(siber.full,
                ax.pad = 2, 
                hulls = F, 
                ellipses = F,
                group.hulls = F,
                bty = "o",
                iso.order = c(1,2),
                xlab = expression({delta}^13*C~'(\u2030)'),
                ylab = expression({delta}^15*N~'(\u2030)')
)

#Add community hulls, computed using the means of every consumer group
plotCommunityHulls(siber.full, plot.args = list(col = "black", lty = 2))

#Those graphs are good to take a look at your data. It's probably not your best option to produce nice figures.
#To do that, we could use ggplot. There are many online resources that will help you customize just about any
#aspect of your graph.

#First, let's calculate the means for each species in each community. To do that, we will use dplyr.
comm_means <- fulldata %>% group_by(community, group) %>%
  summarise(mC = mean(iso1, na.rm = TRUE),
            mN = mean(iso2, na.rm = TRUE))

#Now we can find the limits of our convex hulls for each community.
hulls_comm <- comm_means %>%
  group_by(community) %>%
  slice(chull(mC,mN))

#And plot all that using ggplot.
(hull_plot <- ggplot(comm_means, aes (x = mC, y = mN, colour = factor(community))) +
    geom_point(size = 2) +
    ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
    xlab(expression(paste(delta^{13}, "C (\u2030)"))) +
    scale_colour_manual(name= "Community", labels = c("Cap des Elephants", "Anse du Lion"), 
                        values = c("#377EB8","#E41A1C")) +
    scale_fill_manual(name = "Community", labels = c("Cap des Elephants", "Anse du Lion"), 
                      values = c("#377EB8","#E41A1C")) +
    aes(fill = factor(community)) +
    geom_polygon(data = hulls_comm, alpha = 0.3) +
    theme_bw())

#Now you have a nice-looking graph that you can export to illustrate documents. Time to calculate some metrics.

#Calculate the various Layman metrics on each of the communities.
Layman.full.calc <- communityMetricsML(siber.full) 

#See those metrics in the R console
print(Layman.full.calc)

#As explained during the theoretical course, we can also use a Bayesian model to estimate these values.
#To do that, we first need to set the options for running JAGS. Details are outside the scope of this course.
parms <- list()
parms$n.iter <- 2 * 10^4   # Number of model iterations. Here, we will keep it low to limit computation time.
parms$n.burnin <- 1 * 10^3 # Number of initial discarded values
parms$n.thin <- 10     # Interval to thin the posterior
parms$n.chains <- 2        # NUmber of chains to run

# We also need to define the priors. Once again, details are outside the scope of this course.
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3

#Now we will run the model, using the parameters and priors we just specified.
#This might take a while, and the run length will depend on the computer you use.
ellipses.post.full <- siberMVN(siber.full, parms, priors)

#Now we can extract the posterior means from our model, and calculate the corresponding distribution of layman 
#metrics. Once again, this takes a while.
mu.post <- extractPosteriorMeans(siber.full, ellipses.post.full)
layman.full.bayes <- bayesianLayman(mu.post)

#Let's produce density plots of the metrics to compare them across communities. Given the ranges that can be very
#different, it's more efficient to do it separately for each metric. Let's try it with TA (total area of the 
#convex hull).
par(mgp=c(3,0.6,0))
siberDensityPlot(cbind(layman.full.bayes[[1]][,"TA"], layman.full.bayes[[2]][,"TA"]),
                 xticklabels = c("Anse du Lion", "Cap des Elephants"), 
                 bty="L", ylim = c(0,60),
                 las = 1,
                 ylab = "Total area of the convex hull",
                 xlab = "Community")

#On this plot, by default, the dark, median, and light grey are the 50, 75 and 95% credibility intervals, and the
#black dot is the mode. We can also add the geometrically calculated value to compare our two methods.
points(1:2, Layman.full.calc[3,], col = "red", pch = 16)

#If we want, we can do that with other metrics. Here is an example with the standard deviation of nearest neighbour
#distance ("SDNND").
siberDensityPlot(cbind(layman.full.bayes[[1]][,"SDNND"], layman.full.bayes[[2]][,"SDNND"]),
                 xticklabels = c("Anse du Lion", "Cap des Elephants"), 
                 bty="L", ylim = c(0,1),
                 las = 1,
                 ylab = "SDNND",
                 xlab = "Community")
points(1:2, Layman.full.calc[6,], col = "red", pch = 16)

#Finally, we can use our model to compute the probability that Layman metrics' values are different between 
#communities. For TA: 
Prob.diff.TA <- sum(layman.full.bayes[[1]][,"TA"]>layman.full.bayes[[2]][,"TA"])/ NROW(layman.full.bayes[[2]][,"TA"])

#Here, this probability is inferior to 95%. It means that there are less than 95% of model runs where the TA of 
#community 1 is bigger than the one of community 2. If we do an analogy with traditional hypothesis test, it means
#that TA should be considered similar in the two communities.

#Now let's compare another metric, SDNND.
Prob.diff.SDNND <- sum(layman.full.bayes[[1]][,"SDNND"]>layman.full.bayes[[2]][,"SDNND"])/ 
  NROW(layman.full.bayes[[2]][,"SDNND"])

#The probability is around 66%, which is much lower than for TA. Once again, values of this metric can be considered
#similar in the two communities.



######################################################################################################################
#PART 2 - STANDARD ELLIPSE METRICS ON POPULATIONS
######################################################################################################################

#For this part, and to keep it simple, we will only work with two species: the sea stars Diplasterias brucei and

species_groups <- c(2, 3, 5, 14, 15, 16)
DipOdo <- fulldata[fulldata$group %in% species_groups & fulldata$community == 1, ]


#Now we can create our SIBER object.
siber.dipodo <- createSiberObject(DipOdo)

#And quickly visualize it.
par(mar=c(5.1, 5.1, 4.1, 2.1))
group.ellipses.args  <- list(n = 100, lty = 1, lwd = 2)
plotSiberObject(siber.dipodo,
                ax.pad = 2, 
                hulls = F, 
                ellipses = T,
                group.hulls = F,
                bty = "o",
                iso.order = c(1,2),
                xlab = expression({delta}^13*C~'(\u2030)'),
                ylab = expression({delta}^15*N~'(\u2030)')
)

#Like before, SIBER itself is not the best option to customize graphs. Instead, we can do it using ggplot.

#First, let's create a function that will extract ellipse values out of the SIBER object
pullEllipseStd = function(data, x, y) {
  return(as.data.frame(addEllipse(data$ML.mu[[x]][ , , y],
                                  data$ML.cov[[x]][ , , y],
                                  m = NULL,
                                  n = 100,
                                  p.interval = NULL,
                                  ci.mean = FALSE, do.plot = FALSE)))
}




#Now we'll use that function on each of our groups. The x argument is the community number, and the y argument
#is the group number, taken sequentially
# Crie as elipses para as espécies adicionais da comunidade 1
Ellipse_Char_Lions = pullEllipseStd(data=siber.dipodo, x=1, y=1)  # Para Charcotia obesa
Ellipse_Dec_Lions = pullEllipseStd(data=siber.dipodo, x=1, y=2)   # Para Decolopoda australis
Ellipse_Dip_Lions = pullEllipseStd(data=siber.dipodo, x=1, y=3)   # Para Diplasterias brucei (já estava no seu código)
Ellipse_Odo_Lions = pullEllipseStd(data=siber.dipodo, x=1, y=4)   # Para Odontaster validus (já estava no seu código)
Ellipse_Oph_Lions = pullEllipseStd(data=siber.dipodo, x=1, y=5)   # Para Ophiura sp.
Ellipse_Par_Lions = pullEllipseStd(data=siber.dipodo, x=1, y=6)   # Para Parborlasia corrugatus


#Now we can plot that using ggplot. Here as well as different symbols, I used different line types for each
#taxon (dashed lines for D. brucei, solid lines for O. validus)
# Definindo as cores para cada espécie
species_colors <- c("Charcotia obesa" = "#377EB8", 
                    "Decolopoda australis" = "#E41A1C", 
                    "Diplasterias brucei" = "#4DAF4A",
                    "Odontaster validus" = "#984EA3",
                    "Ophiura sp." = "green", 
                    "Parborlasia corrugatus" = "orange")

# Criar uma nova coluna em DipOdo com os nomes das espécies para usar na legenda
DipOdo$species <- factor(DipOdo$group, labels = c("Charcotia obesa", "Decolopoda australis", 
                                                  "Diplasterias brucei", "Odontaster validus", 
                                                  "Ophiura sp.", "Parborlasia corrugatus"))

# Agora usando a coluna species para a estética colour
ggplot(DipOdo, aes(x=iso1, y=iso2)) +
  geom_point(aes(colour=species, shape=species), size=2) +
  scale_colour_manual(values=species_colors) +
  scale_shape_manual(values=1:20) +
  geom_path(data=Ellipse_Dip_Lions, aes(x=V1, y=V2, colour="Diplasterias brucei"), size=1.25, linetype="solid") +
  geom_path(data=Ellipse_Char_Lions, aes(x=V1, y=V2, colour="Charcotia obesa"), size=1.25, linetype="solid") +
  geom_path(data=Ellipse_Dec_Lions, aes(x=V1, y=V2, colour="Decolopoda australis"), size=1.25, linetype="solid") +
  geom_path(data=Ellipse_Odo_Lions, aes(x=V1, y=V2, colour="Odontaster validus"), size=1.25, linetype="solid") +
  geom_path(data=Ellipse_Oph_Lions, aes(x=V1, y=V2, colour="Ophiura sp."), size=1.25, linetype="solid") +
  geom_path(data=Ellipse_Par_Lions, aes(x=V1, y=V2, colour="Parborlasia corrugatus"), size=1.25, linetype="solid") +
  scale_x_continuous(name=expression({delta}^13*C~('\u2030'))) +
  scale_y_continuous(name=expression({delta}^15*N~('\u2030'))) +
  theme_bw() +
  theme(legend.position="right") +
  guides(colour=guide_legend("Taxon"))
#Now that we have a pretty graph, let's have a look at the standard ellipse areas.
group.ML <- groupMetricsML(siber.dipodo)
print(group.ML)

#We will now use a Bayesian model to estimate SEA. To do that, we need parameters and priors like in the
#part 1 of the script. We will use the ones we defined earlier and run the model
ellipses.post.dipodo <- siberMVN(siber.dipodo, parms, priors)

#Now we can sample our posterior distribution.
SEA.B <- siberEllipses(ellipses.post.dipodo)

#Let's plot our model's output
siberDensityPlot(SEA.B, xticklabels = colnames(group.ML), 
                 xlab = c("Community | Group"),
                 ylab = expression("Standard Ellipse Area  " ('\u2030' ^2) ),
                 bty = "L",
                 las = 1,
                 main = "Model estimations of SEA",
                 ylims = c(0,9)
)

#And add red dots for the SEA computed earlier from the covariance matrix
points(1:6, group.ML[3,], col = "red", pch = 16)


######################################################################################################################
#END OF SCRIPT
######################################################################################################################
