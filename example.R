# An example application of our algorithm (see section 4)

# Read in data and divide between x, y, and w
dat <- read.csv("NHANES_07_12.csv")
dat.x <- dat$SD.level
dat.y <- dat$bmi
dat.w <- dat[c('age','gender','fish.score','vitD.supplement','MET.score','white','black')]

# Create interaction terms
p <- ncol(dat.w)
for(var1 in 1:p) {
  for (var2 in var1:p) {
    if (var1 == var2) {
      next
    } else {
      var.name <- paste0(colnames(dat.w)[c(var1, var2)], collapse = '_')
    }
    dat.w[var.name] <- dat.w[var1] * dat.w[var2]
  }
}
# Eliminate interaction between racial indicator variables
dat.w <- subset(dat.w, select = !(colnames(dat.w) %in% 'white_black'))

# Run the BBR algorithm
BB.confound.reorder(dat.x, dat.y, dat.w)

# Note: Because the data set is large, it may take a while for the algorithm to finish.
# In fact, the BB.confound algorithm never finished after 200 minutes of runtime.
