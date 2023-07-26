library(ggplot2)
library(dplyr)
library(ggpubr)
library(BSDA)

#change working directory to data folder:
setwd("~/Box Sync/FLL/2023")
all_data <- read.csv(file.choose())
View(all_data)
all_data <- data.frame(all_data)
View(all_data)
random <- all_data$rand_val
location <- all_data$loc_val

#test for normality
hist(all_data$rand_val)
hist(all_data$loc_val)

random_sd <- sd(random, na.rm = TRUE)
location_sd <-sd(location, na.rm = TRUE)

View(random_sd)
View(location_sd)


par(mar = c(10, 5, 5, 5))

boxplot(random,location,
        #main = "Initial vs Interim, fasting food",
        names = c("random points", "occurence points"),
        las = 2,
        ylab = "Normalized Cumulative Current (Omniscape)",
        col = c("lightskyblue2", "lawngreen"),
        #        horizontal = TRUE,
        notch = FALSE
)

z.test(random, location, sigma.x=0.139, sigma.y=0.125, alternative="two.sided", mu=0)
