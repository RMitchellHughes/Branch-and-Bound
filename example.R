# An example application of our algorithm (see section 4)

# Read in data and divide between x, y, and s
dat <- read.csv("NHANES_07_12.csv")
dat.x <- dat$SD.level
dat.y <- dat$bmi
dat.s <- dat[c('age','gender','fish.score','vitD.supplement','MET.score','white','black')]

# Create interaction terms
p <- ncol(dat.s)
for(var1 in 1:p) {
  for (var2 in var1:p) {
    if (var1 == var2) {
      next
    } else {
      var.name <- paste0(colnames(dat.s)[c(var1, var2)], collapse = '_')
    }
    dat.s[var.name] <- dat.s[var1] * dat.s[var2]
  }
}
# Eliminate interaction between racial indicator variables
dat.s <- subset(dat.s, select = !(colnames(dat.s) %in% 'white_black'))

# Run the BB algorithm
BB.confound(dat.x, dat.y, dat.s)

# Note: Because the data set is large, it may take a while for the algorithm to finish.
