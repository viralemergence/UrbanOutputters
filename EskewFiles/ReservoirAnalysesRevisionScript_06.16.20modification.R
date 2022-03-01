

# Evan Eskew
# 02.22.2017

# This script was created to analyze reservoir data compiled by 
# Ben Plourde and the Foley lab. 

# The input data consist of two datasheets:
# One contains life history data for all mammals
# The other consists of data for mammalian reservoir species

library(dplyr)
library(pca3d)

#==============================================================================


# Import trait databases with imputed data (all species and only reservoirs)
d <- read.csv("Data/Fulldata_Imputed_FIXED.csv", na.string = c("NA", "NA "))
d2 <- read.csv("Data/res.imp.FIXED.csv", na.strings = c("NA", "NA "))
d2$MSW05_Order <- as.factor(d2$MSW05_Order)


# Modify the all species trait database to include log-transformed variables
# for key life history characteristics
d <- d %>%
  mutate(log_AdultBodyMass = log(AdultMass_g)) %>%
  mutate(log_GestationLen = log(Gest_d)) %>%
  mutate(log_LitterSize = log(LitterSize)) %>%
  mutate(log_NeonateBodyMass = log(BirthWt_g)) %>%
  mutate(log_InterbirthInterval = log(IBI_d)) %>%
  mutate(log_WeaningAge = log(Wean_d)) %>%
  mutate(log_SexualMaturityAge = log(SexMat_d))


# Add trait residuals to the d dataframe
d$res.GestationLen <- lm(d$log_GestationLen ~ d$log_AdultBodyMass)$residuals
d$res.LitterSize <- lm(d$log_LitterSize ~ d$log_AdultBodyMass)$residuals
d$res.NeonateBodyMass <- 
  lm(d$log_NeonateBodyMass ~ d$log_AdultBodyMass)$residuals
d$res.InterbirthInterval <- 
  lm(d$log_InterbirthInterval ~ d$log_AdultBodyMass)$residuals
d$res.WeaningAge <- 
  lm(d$log_WeaningAge ~ d$log_AdultBodyMass)$residuals
d$res.SexualMaturityAge <- 
  lm(d$log_SexualMaturityAge ~ d$log_AdultBodyMass)$residuals


# Modify the reservoir trait database to include log-transformed variables 
# and eliminate rows that have species as "spp." 
# Note that here you could also subset the reservoir data down. For example, 
# only pathogens with human targets
d2 <- d2 %>%
  mutate(log_AdultBodyMass = log(AdultMass_g)) %>%
  mutate(log_GestationLen = log(Gest_d)) %>%
  mutate(log_LitterSize = log(LitterSize)) %>%
  mutate(log_NeonateBodyMass = log(BirthWt_g)) %>%
  mutate(log_InterbirthInterval = log(IBI_d)) %>%
  mutate(log_WeaningAge = log(Wean_d)) %>%
  mutate(log_SexualMaturityAge = log(SexMat_d)) %>%
  filter(MSW05_Species != "spp." & MSW05_Species != "spp. ") %>%
  # filter(human.target == 1) %>%
  droplevels()


# Add residuals to the d2 dataframe by filling in the values from the 
# imputed d dataframe
for (row in 1:nrow(d2)){
  d2$res.GestationLen[row] = 
    d[d$binom == as.character(d2$species[row]),]$res.GestationLen
  d2$res.LitterSize[row] = 
    d[d$binom == as.character(d2$species[row]),]$res.LitterSize
  d2$res.NeonateBodyMass[row] =
    d[d$binom == as.character(d2$species[row]),]$res.NeonateBodyMass
  d2$res.InterbirthInterval[row] = 
    d[d$binom == as.character(d2$species[row]),]$res.InterbirthInterval
  d2$res.WeaningAge[row] = 
    d[d$binom == as.character(d2$species[row]),]$res.WeaningAge
  d2$res.SexualMaturityAge[row] = 
    d[d$binom == as.character(d2$species[row]),]$res.SexualMaturityAge
}

#==============================================================================


# Create a function to calculate a permutation test comparing random samples 
# of mean values from the trait database vs. the mean trait values for the
# reservoir species. The inputs are the trait of interest, the number of
# bootstrapped means to compute, and the x-axis limits for graphing the 
# density plots (the different traits sometimes need different values).
trait_test <- function(trait, n.means, x, y){

set.seed(8)

# Print trait (mainly for debugging purposes) and info about input dataframes
print(trait)
print(paste("Dimensions of full species database:", toString(dim(d))))
print(paste("Dimensions of reservoir database:", toString(dim(d2))))

# Find the columns representing the trait of interest in both the 
# d and d2 dataframes
d.col <- which(colnames(d) == trait)
d2.col <- which(colnames(d2) == trait)

# Select colors for plotting (col1 for reservoirs, col2 for other mammals)
col1 <- "black"
col2 <- "dimgrey"

# Draw mean trait values from the trait database (equal in number to "n.means").
# Each mean is computed from a randomly sampled number of species equal to the 
# number of unique species in the reservoir database. Note that this was our 
# initial idea but have since revised our analysis (see below)
samples <- sapply(1:n.means, function(z) 
  mean(d[,d.col][sample(1:nrow(d), length(unique(d2$species)))], na.rm = T))

# Or create a loop that will generate mean trait values while only drawing from
# Orders that are represented in the reservoir database. In this case, it's a
# bit more complicated. We still want "n.means" number of samples, but each
# sample should be the mean of the trait from a random subset of the trait
# database that has the same taxonomic structure as our reservoir database 
# (i.e., the same Orders are represented the same number of times).

# Set up a vector to hold the new samples
samples2 <- rep(NA, n.means)

# Create a loop to fill in these sampled mean trait values
for (i in 1:n.means) {
temp.df <- c()

# This inner loop works by grabbing a number of rows from the d dataframe to
# create a temporary dataframe that has the exact same taxonomic structure 
# as the d2 dataframe
for (order.name in levels(d2$MSW05_Order)){
  temp.d <- filter(d, Order == order.name)
  temp.d2 <- filter(d2, MSW05_Order == order.name)
  temp.d <- temp.d[sample(1:nrow(temp.d), length(unique(temp.d2$species))), ]
  temp.df <- rbind(temp.df, temp.d)
}

# We then summarize the trait value for this temporary dataframe, report 
# some stats about the progress of the loop, and repeat the process 
# "n.means" number of times
samples2[i] <- mean(temp.df[,d.col], na.rm = T)
print(paste("Sample:", i))
print(paste("Dimensions of temporary dataframe:", toString(dim(temp.df))))
print(paste("Number of missing values in temporary dataframe:", 
            sum(is.na(temp.df[,d.col]))))
}

# Plot the random samples along with vertical lines indicating the 2.5th 
# and 97.5th percentile

# First define x-axis labels based on the trait that is being plotted
xlabel <- trait
if(trait == "log_AdultBodyMass") {xlabel = "log(Adult Body Mass (g))"}
if(trait == "log_GestationLen") {xlabel = "log(Gestation Length (days))"}
if(trait == "log_LitterSize") {xlabel = "log(Litter Size)"}
if(trait == "log_NeonateBodyMass") {xlabel = "log(Neonate Body Mass (g))"}
if(trait == "log_InterbirthInterval") 
  {xlabel = "log(Interbirth Interval (days))"}
if(trait == "log_WeaningAge") {xlabel = "log(Weaning Age (days))"}
if(trait == "log_SexualMaturityAge") 
  {xlabel = "log(Sexual Maturity Age (days))"}
if(trait == "res.GestationLen") 
  {xlabel = "Residual of Gestation Length (days)"}
if(trait == "res.LitterSize") {xlabel = "Residual of Litter Size"}
if(trait == "res.NeonateBodyMass") 
  {xlabel = "Residual of Neonate Body Mass (g)"}
if(trait == "res.InterbirthInterval") 
  {xlabel = "Residual of Interbirth Interval (days)"}
if(trait == "res.WeaningAge") {xlabel = "Residual of Weaning Age (days)"}
if(trait == "res.SexualMaturityAge") 
  {xlabel = "Residual of Sexual Maturity Age (days)"}

# Initially I plotted results using a histogram, but I believe a density plot 
# is much preferable
# hist(samples2, col = col2, xlab = xlabel, main = "", xlim = c(x, y), 
# ylim = c(0, n.means/5), breaks = seq(from = x, to = y, length.out = 100))

# Plot a density kernel plot including a legend and vertical lines 
# representing the 2.5th and 97.5th percentiles
plot(density(samples2), xlab = xlabel, main = "", xlim = c(x, y), las = 1,
     ylim = c(0, ceiling(max(density(samples2)$y))), yaxs = "i", bty = "l")
polygon(density(samples2), col = adjustcolor(col2, alpha.f = 0.2), 
        border = col2, lwd = 2)
# legend("topleft", pch = 15, col = c(col1, col2), 
# legend = c("Reservoirs", "Representative Mammals"), bty = "n", cex = 1)
abline(v = sort(samples2)[0.025 * length(samples2)], 
       col = adjustcolor(col2, alpha.f = 0.4), lwd = 2, lty = 2)
abline(v = sort(samples2)[0.975 * length(samples2)], 
       col = adjustcolor(col2, alpha.f = 0.4), lwd = 2, lty = 2)

# Create a dataframe (d3) that represents only unique species found in d2
d3 <- c()
for (species in unique(d2$species)){
  row <- d2[which(d2$species == species)[1], ]
  d3 <- rbind(d3, row)
}

print(paste("Dimensions of reservoir database subset down to unique species:",
            toString(dim(d3))))
print(paste("Number of missing values:", sum(is.na(d3[,d2.col]))))

# Plot the mean value of the trait for reservoir species as a vertical line.
# Note that because we are calculating the mean on the d3 dataframe, we are
# excluding duplicated rows in d2 that represent the exact same reservoir
# species that many be invovled in multiple diseases
abline(v = mean(d3[ , d2.col], na.rm = T), col = col1, lwd = 3, lty = 1)

# Print a summary of the computed samples
print(summary(samples))
print(summary(samples2))
}

#==============================================================================


# Run trait tests on log-transformed life history variable residuals and 
# save the output to a pdf
sink("Trait_plots_ResidualsWImputation_log.txt")
pdf("Trait_plots_ResidualsWImputation.pdf", width = 8, height = 5)
n.iterations <- 1000

par(mfrow = c(2,3))

trait_test("res.GestationLen", n.iterations, -0.5, 0.5)
# legend("topright", legend = "A", bty = "n", cex = 2)
trait_test("res.LitterSize", n.iterations, -0.5, 0.5)
# legend("topright", legend = "B", bty = "n", cex = 2)
trait_test("res.NeonateBodyMass", n.iterations, -0.5, 0.5)
# legend("topright", legend = "C", bty = "n", cex = 2)
trait_test("res.InterbirthInterval", n.iterations, -0.5, 0.5)
# legend("topright", legend = "D", bty = "n", cex = 2)
trait_test("res.WeaningAge", n.iterations, -0.5, 0.5)
# legend("topright", legend = "E", bty = "n", cex = 2)
trait_test("res.SexualMaturityAge", n.iterations, -0.5, 0.5)
# legend("topright", legend = "F", bty = "n", cex = 2)

sink()
dev.off()

#==============================================================================


# PCA Analysis

# Create a dataframe representing only the columns of interest for 
# the PCA analysis
PCA.df <- select(d, log_AdultBodyMass, log_GestationLen, log_LitterSize,
                 log_NeonateBodyMass, log_InterbirthInterval, log_WeaningAge,
                 log_SexualMaturityAge)

# Note this dataframe needs to be subsetted down to only complete data for 
# the PCA analysis
nrow(PCA.df)
sum(complete.cases(PCA.df))
PCA.df <- PCA.df[complete.cases(PCA.df), ]


# Compute residuals on adult body mass because we want to know the effect of
# those other life history factors, controlling for the effect of adult body
# mass. Add those to the PCA.df dataframe
fit1 <- lm(PCA.df[,2] ~ PCA.df[,1])
fit2 <- lm(PCA.df[,3] ~ PCA.df[,1])
fit3 <- lm(PCA.df[,4] ~ PCA.df[,1])
fit4 <- lm(PCA.df[,5] ~ PCA.df[,1])
fit5 <- lm(PCA.df[,6] ~ PCA.df[,1])
fit6 <- lm(PCA.df[,7] ~ PCA.df[,1])

PCA.df$res.GestationLen <- fit1$residuals
PCA.df$res.LitterSize <- fit2$residuals
PCA.df$res.NeonateBodymass <- fit3$residuals
PCA.df$res.InterbirthInterval <- fit4$residuals
PCA.df$res.WeaningAge <- fit5$residuals
PCA.df$res.SexualMaturityAge <- fit6$residuals


# Perform the PCA, choosing PCA.df columns of interest (the residuals)
PCA.result <- prcomp(PCA.df[,8:13], center = T, scale. = T)


# Prepare the d and d2 dataframes so trait tests can be run on the PCA variables

# Subset d down to complete cases and add on the PCA scores
d <- d[complete.cases(select(d, log_AdultBodyMass, log_GestationLen, 
                            log_LitterSize, log_NeonateBodyMass, 
                            log_InterbirthInterval, log_WeaningAge, 
                            log_SexualMaturityAge)), ]
d$PC1 <- PCA.result$x[,1]
d$PC2 <- PCA.result$x[,2]
d$PC3 <- PCA.result$x[,3]
d$PC4 <- PCA.result$x[,4]
d$PC5 <- PCA.result$x[,5]
d$PC6 <- PCA.result$x[,6]
d <- droplevels(d)

# Subset d2 down to complete cases, and add on PCA scores taken from the 
# d dataframe
d2 <- d2[complete.cases(select(d2, log_AdultBodyMass, log_GestationLen, 
                              log_LitterSize, log_NeonateBodyMass, 
                              log_InterbirthInterval, log_WeaningAge, 
                              log_SexualMaturityAge)), ]
for (i in 1:nrow(d2)){
  d2$PC1[i] <- d$PC1[which(d$binom == toString(d2$species[i]))]
  d2$PC2[i] <- d$PC2[which(d$binom == toString(d2$species[i]))]
  d2$PC3[i] <- d$PC3[which(d$binom == toString(d2$species[i]))]
  d2$PC4[i] <- d$PC4[which(d$binom == toString(d2$species[i]))]
  d2$PC5[i] <- d$PC5[which(d$binom == toString(d2$species[i]))]
  d2$PC6[i] <- d$PC6[which(d$binom == toString(d2$species[i]))]}
d2 <- droplevels(d2)


# Output PCA results, plot PCAs, and run trait tests on PCA variables 

# Create a color vector for plotting PCA scores according to reservoir status
is.Reservoir <- sapply(1:nrow(d), function(z) 
  ifelse(sum(d$binom[z] == as.character(d2$species)) > 0, 1, 2))
plot.colors <- c("black", adjustcolor("dimgrey", alpha.f = 0.6))
color.vector <- plot.colors[is.Reservoir]

print(PCA.result)
plot(PCA.result, type = "l")
summary(PCA.result)

# Plot basic PCA biplot
biplot(PCA.result)

# Plot pretty PCA biplot
pca2d(PCA.result, col = color.vector, shape = 19, biplot = T)


# Plot a pretty PCA
pdf("PCA_plot_WImputation.pdf", width = 8, height = 5)

pca2d(PCA.result, col = color.vector, shape = 19)
points(PCA.result$x[is.Reservoir == 1, 1], 
       PCA.result$x[is.Reservoir == 1, 2], 
       col = plot.colors[1], pch = 19)
legend("topleft", pch = 19, col = plot.colors, 
       legend = c("Reservoirs", "Other Mammals"), bty = "n", cex = 1.1)

dev.off()


# Plot trait tests on PCA variables
sink("PCA_traitplots_WImputation_log.txt")
pdf("PCA_traitplots_WImputation.pdf", width = 8, height = 5)
n.iterations <- 1000

par(mfrow = c(1,2))

trait_test("PC1", n.iterations, -1.5, 0)
trait_test("PC2", n.iterations, -1, 0.5)
# trait_test("PC3", n.iterations, -1, 0.5)
# trait_test("PC4", n.iterations, -1, 0.5)
# trait_test("PC5", n.iterations, -0.5, 0.3)
# trait_test("PC6", n.iterations, -0.5, 0.3)

dev.off()
sink()

#==============================================================================


# Factor Analysis
# Can simply use the PCA.df dataframe from the previous analysis

# Perform the factor analysis with two extracted variables
factor.result <- factanal(PCA.df[,8:13], 2, rotation = "varimax")

# Show results
print(factor.result)
plot(factor.result$loadings[,1:2], type = "n", 
     xlim = c(-0.5,1), ylim = c(-0.5,1))
text(factor.result$loadings[,1:2], labels = names(PCA.df[,8:13]), cex = 0.8)
