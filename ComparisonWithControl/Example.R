#######################################################################
## Generate the data frame 
## Read in data
day2 <- read.csv("Day2.csv", header = F)

## Add covariates
d2 <- as.data.frame(as.vector(as.matrix(day2)))
d2$dose <- rep(c(rep(0, 7), rep(10, 8), rep(50, 7), 100), 6)
d2$person <- c(rep("Mahi", nrow(day2)*2), rep("Nicole", nrow(day2)*2), 
               rep("Xochilt", nrow(day2)*2))
d2$side <- rep(c(rep("LV",nrow(day2)), rep("RV",nrow(day2))), 3)
## Change column names of the response variable.
colnames(d2)[1] <- "response"


## ANOVA test
a2 <- aov(response ~ dose + person + side, data = d2)
summary(a2)


## TukeyHSD test
d2$dose <- as.factor(d2$dose)
TukeyHSD(a2, which = "dose")



## Dunnett test
library(DescTools)
DunnettTest(response ~ dose, data = d2, control="0", conf.level=0.95)