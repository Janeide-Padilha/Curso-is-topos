######################################################################################################################
#PART 1 - CONVEX HULL METRICS ON COMMUNITIES
######################################################################################################################

library(SIBER)
library(dplyr)
library(ggplot2)
library(rjags)


data <- read.table(text = " iso1	iso2	group	community
-23.89	8.83	1	1
-24.29	9.53	1	1
-24.47	9.72	1	1
-24.65	9.37	1	1
-25.54	9.85	1	1
-24.48	9.16	1	1
-24.78	9.89	1	1
-24.59	9.99	1	1
-26.41	9.49	1	1
-24.84	9.18	1	1
-24.73	9.64	1	1
-24.74	9.34	1	1
-25.51	9.55	1	1
-24.44	9.25	2	1
-26.44	9.27	2	1
-25.26	9.54	2	1
-25.71	9.41	2	1
-25.30	9.77	2	1
-24.74	9.52	2	1
-25.91	9.02	2	1
-25.36	8.65	2	1
-24.80	8.29	2	1
-25.18	9.43	2	1
-24.48	9.83	2	1
-26.84	8.44	2	1
-25.23	9.3	2	1
-25.83	8.3	2	1
-24.81	9.99	2	1
-22.98	9.29	3	1
-25.24	9.76	3	1
-24.02	10.22	3	1
-24.50	10.06	3	1
-25.83	10.66	3	1
-25.28	9.9	3	1
-17.02	10.25	4	2
-19.54	9.78	4	2
-22.35	10.82	4	2
-20.87	11.32	4	2
-20.87	11.32	4	2
-16.76	15.42	5	2
-26.17	8.08	5	2
-25.27	9.17	5	2
-20.89	11.22	5	2
-21.86	10.10	5	2
-25.49	8.32	5	2
-24.18	9.15	5	2
-25.05	8.43	5	2
-22.96	8.56	5	2
-24.59	8.37	5	2
-25.22	9.12	5	2
-24.53	7.85	5	2
-26.17	8.08	5	2
-22.96	8.56	5	2
-24.07	9.35	5	2
-24.06	8.73	5	2
-20.92	13.23	6	2
-22.02	12.37	6	2
-21.34	12.69	6	2
-20.96	13.67	6	2
-21.79	12.82	6	2
-21.22	13.08	6	2
-21.12	13.46	6	2
-19.52	13.43	6	2
-23.58	11.41	6	2
-21.76	12.88	6	2
-21.74	12.89	6	2
-21.09	13.17	6	2
-21.37	12.89	6	2
-20.84	13.73	6	2
-21.34	12.69	6	2
-15.10	11.37	6	2
-21.24	12.83	6	2
-21.97	12.61	6	2
-23.58	11.41	6	2
-21.22	13.08	6	2
-21.62	12.99	6	2
-20.95	13.25	6	2
-16.63	12.62	7	2
-17.42	12.13	7	2
-18.84	12.34	7	2
-17.24	12.75	7	2
-17.78	14.10	7	2
-17.53	12.61	7	2
-17.88	13.35	7	2
-17.83	13.68	7	2
-19.00	11.25	7	2
-17.40	12.67	7	2
-18.07	11.15	8	2
-17.88	11.77	8	2
-20.08	11.56	8	2
-18.85	10.58	8	2
-17.39	10.43	8	2", header = TRUE, sep = "\t")


# Remover colunas desnecessárias
fulldata <- data %>% select(iso1, iso2, group, community)
key <- read.table(text = "community	1	Residente
	2	Migratória
group	1	Ppa
	2	Pan
	3	Pad
	4	Calb
	5	Svi
	6	Mgi
	7	Sma
	8	Ldo", header = TRUE, sep = "\t")

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
    scale_colour_manual(name= "Comunidade", labels = c("Residente", "Migratória"), 
                        values = c("#377EB8","#E41A1C")) +
    scale_fill_manual(name = "Comunidade", labels = c("Residente", "Migratória"), 
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
                 xticklabels = c("Residente", "Migratória"), 
                 bty="L", ylim = c(0,60),
                 las = 1,
                 ylab = "Área total da envoltória convexa",
                 xlab = "Comunidade")

#On this plot, by default, the dark, median, and light grey are the 50, 75 and 95% credibility intervals, and the
#black dot is the mode. We can also add the geometrically calculated value to compare our two methods.
points(1:2, Layman.full.calc[3,], col = "red", pch = 16)

#If we want, we can do that with other metrics. Here is an example with the standard deviation of nearest neighbour
#distance ("SDNND").
siberDensityPlot(cbind(layman.full.bayes[[1]][,"SDNND"], layman.full.bayes[[2]][,"SDNND"]),
                 xticklabels = c("Residente", "Migratória"), 
                 bty="L", ylim = c(0,1),
                 las = 1,
                 ylab = "SDNND",
                 xlab = "Comunidade")
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


DipOdo <- fulldata %>% filter(group == "1" | group == "2"| group == "3"| group == "4"| group == "5"| group == "6"| group == "7"| group == "8")

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
# Now we'll use that function on each of our groups. The x argument is the community number, and the y argument
# Crie elipses para todos os 8 grupos
# Resident community
Ellipse_Resident_1 <- pullEllipseStd(data=siber.dipodo, x=1, y=1)
Ellipse_Resident_2 <- pullEllipseStd(data=siber.dipodo, x=1, y=2)
Ellipse_Resident_3 <- pullEllipseStd(data=siber.dipodo, x=1, y=3)

# Migratory community
Ellipse_Migratory_1 <- pullEllipseStd(data=siber.dipodo, x=2, y=1)
Ellipse_Migratory_2 <- pullEllipseStd(data=siber.dipodo, x=2, y=2)
Ellipse_Migratory_3 <- pullEllipseStd(data=siber.dipodo, x=2, y=3)
Ellipse_Migratory_4 <- pullEllipseStd(data=siber.dipodo, x=2, y=4)
Ellipse_Migratory_5 <- pullEllipseStd(data=siber.dipodo, x=2, y=5)



# Atualize o gráfico ggplot
# Atualize o gráfico ggplot
# Atualize o gráfico ggplot
plot <- ggplot(DipOdo, aes(x=iso1, y=iso2)) +
  geom_point(aes(col=factor(community), shape=factor(group)), size=2) +
  scale_colour_manual(name = "Comunidade", labels = c("Residente", "Migratória"), values = c("#377EB8", "#E41A1C")) +
  scale_shape_manual(name ="Éspecies", labels = c("Ppa", "Pan", "Pad", "Calb", "Svi", "Mgi", "Sma", "Ldo"), values=c(1, 2, 3, 4, 5, 6, 7, 8)) +
  geom_path(data=Ellipse_Resident_1, aes(x=V1, y=V2), colour="#377EB8", size=1.25, linetype="dashed") +
  geom_path(data=Ellipse_Resident_2, aes(x=V1, y=V2), colour="#377EB8", size=1.25, linetype="dashed") +
  geom_path(data=Ellipse_Resident_3, aes(x=V1, y=V2), colour="#377EB8", size=1.25, linetype="dashed") +
  geom_path(data=Ellipse_Migratory_1, aes(x=V1, y=V2), colour="#E41A1C", size=1.25, linetype="dashed") +
  geom_path(data=Ellipse_Migratory_2, aes(x=V1, y=V2), colour="#E41A1C", size=1.25, linetype="dashed") +
  geom_path(data=Ellipse_Migratory_3, aes(x=V1, y=V2), colour="#E41A1C", size=1.25, linetype="dashed") +
  geom_path(data=Ellipse_Migratory_4, aes(x=V1, y=V2), colour="#E41A1C", size=1.25, linetype="dashed") +
  geom_path(data=Ellipse_Migratory_5, aes(x=V1, y=V2), colour="#E41A1C", size=1.25, linetype="dashed") +
  scale_x_continuous(name=expression({delta}^13*C~('\u2030'))) +
  scale_y_continuous(name=expression({delta}^15*N~('\u2030'))) +
  theme_bw()




# Exiba o gráfico
print(plot)

                                       

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
                 xlab = c("Comunidade | Grupo"),
                 ylab = expression(paste("Área da elipse padrão (", permil, "^2)")),
                 bty = "L",
                 las = 1,
                 main = "Estimativas do modelo de SEA",
                 ylims = c(0,9)
)

#And add red dots for the SEA computed earlier from the covariance matrix
points(1:8, group.ML[3,], col = "red", pch = 16)


#Now, let's compare SEA of D. brucei in the two communities.
Prob.diff.SEA.DB <- sum(SEA.B[,1]>SEA.B[,3])/ NROW(SEA.B[,1])

#Same thing but for O. validus. We can adapt this line to compare any pair of groups.
Prob.diff.SEA.OV <- sum(SEA.B[,2]>SEA.B[,4])/ NROW(SEA.B[,2])

#Another interesting thing to do would be to compare the overlaps between ellipses of both species in each community
#to estimate to which extent they share resources. Let's do that first for community 1:
sea.overlap.comm1 <- maxLikOverlap("1.1", "1.2", siber.dipodo, p.interval = NULL, n = 29)
sea.overlap.comm1
sea.overlap.comm1 <- maxLikOverlap("1.1", "1.3", siber.dipodo, p.interval = NULL, n = 29)
sea.overlap.comm1
sea.overlap.comm1 <- maxLikOverlap("1.2", "1.3", siber.dipodo, p.interval = NULL, n = 29)
sea.overlap.comm1
#We can also estimate this overlap in a relative way (i.e. as a percentage of the total niche area of the two species,
#combined)
sea.overlap.comm1.rel <- sea.overlap.comm1[3]/(sea.overlap.comm1[1]+sea.overlap.comm1[2]-sea.overlap.comm1[3])
sea.overlap.comm1.rel

#Now let's do the same thing for community 1 and 2.
sea.overlap.comm2 <- maxLikOverlap("2.5", "1.1", siber.dipodo, p.interval = NULL, n = 29)
sea.overlap.comm2
sea.overlap.comm2 <- maxLikOverlap("2.5", "1.2", siber.dipodo, p.interval = NULL, n = 29)
sea.overlap.comm2
sea.overlap.comm2 <- maxLikOverlap("2.5", "1.3", siber.dipodo, p.interval = NULL, n = 29)
sea.overlap.comm2

#Finally, we can use our model to generate Bayesian estimations of overlap, and compare the distributions to see how
#probable it is that niche overlap between the two sea stars is bigger in community 2 than in community 1.
#To do that, we first compute Bayesian overlaps for community 1.
bayes.overlap1 <- bayesianOverlap("1.1", "1.2", ellipses.post.dipodo, p.interval= NULL, draws=29, n = 29)

#Then for community 2.
bayes.overlap2 <- bayesianOverlap("2.5", "2.6", ellipses.post.dipodo, p.interval= NULL, draws=29, n = 29)

#And finally we can compare the two posterior distributions to compute our probability.
Prob.diff.overlap <- sum(bayes.overlap1[,3]>bayes.overlap2[,3])/ NROW(bayes.overlap1[,3])

