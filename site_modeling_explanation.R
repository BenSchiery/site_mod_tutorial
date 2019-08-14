rm(list = ls()) # remove everything from the environment
while(length(dev.list()) > 0){dev.off()} # turn off all graphics devices
library(raster) # loads the `raster` package, which needs to be installed first, i.e. run install.packages("raster")

###################
#### functions ####
###################

# generate random folds for a k-fold cross validation
k.fold <- function(n.obs, k){
  v <- rep(1:k, times = ceiling((n.obs / k)))[1:n.obs]
  sample(v)
}

# matthews correlation coefficient
mcc <- function(tp, fp, tn, fn){
  N <- tp + fp + fn + tn
  S <- (tp + fn) / N
  P <- (tp + fp) / N
  (tp / N - S * P) / sqrt(P * S * (1 - P) * (1 - S))
}

# pa = presence/absence data as 0s and 1s, prdxn = continuous data between 0 and 1
opt.threshold <- function(pa, prdxn){
  t.seq <- seq(0, 1, by = 0.001)
  cc <- sapply(X = t.seq,
               FUN = function(t){
                 pdx <- prdxn
                 pdx[pdx < t] <- 0
                 pdx[pdx >= t] <- 1
                 
                 a <- cbind(pa, 2 * pdx)
                 b <- rowSums(a)
                 
                 tp <- sum(b == 3) + 1 # adding 1 to each of these avoids division by zero that sometimes happens
                 fp <- sum(b == 2) + 1
                 fn <- sum(b == 1) + 1
                 tn <- sum(b == 0) + 1
                 
                 mcc(tp = tp,
                     fp = fp,
                     tn = tn,
                     fn = fn)
               })
  # plot(cc ~ t.seq, xlab = "threshold", ylab = "MCC", 
  #      main = "Threshold Optimization", cex = 0.7, pch = 16)
  t.seq[which(cc == max(cc))[1]]
}

# given a matrix and the position of an element in it, return the mean of the elements nearby
mean.adj.vals <- function(m, ro, co, rad = 1){
  ro. <- c((ro - rad):(ro + rad)); ro. <- ro.[ro. > 0 & ro. <= nrow(m)]
  co. <- c((co - rad):(co + rad)); co. <- co.[co. > 0 & co. <= ncol(m)]
  mean(as.numeric(m[ro., co.]), na.rm = T)
}

# returns a random spatial data surface. max.sd controls 
# the "smoothness," with smaller values giving smoother distributions. 
# iso controls how isolated high-probability areas tend to be
rand.surface <- function(n.row, n.col, max.sd = 0.02, iso = 1.4){
  rad <- round((n.row * n.col)^0.5 / 10)
  m <- matrix(rnorm(n.row * n.col), nrow = n.row, ncol = n.col)
  n.cell <- n.row * n.col
  v <- diag(m) # a sample path across the surface
  l <- length(v)
  d <- v[-1] - v[-l]
  # smooth until the sample path isn't too bumpy (small stdev) using large radius (`rad`)
  while(sd(d) > max.sd){
    ind.smp <- sample(1:n.cell)
    for(k in ind.smp){
      j <- ceiling(k / n.row)
      i <- k - (j - 1) * n.row
      m[i,j] <- mean.adj.vals(m, i, j, rad)
    }
    v <- diag(m)
    d <- v[-1] - v[-l]
  }
  # do one round of smoothing with small radius (1)
  for(k in ind.smp){
    j <- ceiling(k / n.row)
    i <- k - (j - 1) * n.row
    m[i,j] <- mean.adj.vals(m, i, j)
  }
  m <- (m - min(m))^iso
  m / max(m)
}

logistic <- function(x){
  1 / (1 + exp(-x))
}

pal <- colorRampPalette(colors = c("#5C33FF", "white", "#FF3333"))

###########################
#### generate surfaces ####
###########################

n.row <- n.col <- 60 # raster resolution; increasing this gets very computationally expensive very fast
n.cell <- n.row * n.col

# we generate some spatial environmental data
elevation <- rand.surface(n.row = n.row, n.col = n.col) * 300 # elevation in meters
ann.precip <- rand.surface(n.row = n.row, n.col = n.col) * 50 # average annual precipitation in cm
pop.dens <- rand.surface(n.row = n.row, n.col = n.col) * 100 # population density in individuals / sq. km
avg.temp <- rand.surface(n.row = n.row, n.col = n.col) * 40 # average daily temperature in degrees C
soil.sal <- rand.surface(n.row = n.row, n.col = n.col) * 5 # soil salinity in deciSiemens / meter

# we'll define site suitability as the following linear combination of our 
# environmental variables. normally this is the very thing we're trying to find when 
# we make a predictive model, so we wouldn't know it from the start. we'll pretend we don't 
# already know it and try to make our predictive models figure it out
# also, usually suitability won't be linearly dependent on our environmental variables, so 
# a nonlinear modeling algorithm might be better. We'll just be using a linear model because it's faster/simpler
suitability <- elevation / 300 + ann.precip / 50 - pop.dens / 100 + avg.temp / 40 - soil.sal / 5
suitability <- (suitability - min(suitability)) / (max(suitability) - min(suitability)) # scale to between 0 and 1

# the surfaces we just made are in the form of matrices. to make them "plot-able", we turn them into rasters
elev.rstr <- raster(elevation)
prec.rstr <- raster(ann.precip)
popd.rstr <- raster(pop.dens)
temp.rstr <- raster(avg.temp)
soil.rstr <- raster(soil.sal)
suit.rstr <- raster(suitability)

ext <- c("xmin" = 43.3917,
         "xmax" = 46.6766,
         "ymin" = 38.8474,
         "ymax" = 41.2357) # the bounding box for Armenia

# raster objects can contain geographical data, like this latitude/longitude bounding box, so we add it to our rasters
extent(elev.rstr) <- 
  extent(prec.rstr) <- 
  extent(popd.rstr) <- 
  extent(temp.rstr) <- 
  extent(soil.rstr) <- 
  extent(suit.rstr) <- ext 

p4s <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0") # a string specifying a map projection

proj4string(elev.rstr) <- 
  proj4string(prec.rstr) <- 
  proj4string(popd.rstr) <- 
  proj4string(temp.rstr) <- 
  proj4string(soil.rstr) <- 
  proj4string(suit.rstr) <- p4s # set the map projection of our rasters to which their coordinates are relative

############################
#### scatter some sites ####
############################

# uniformly distributed potential site locations
pot.sites <- cbind("long" = runif(n = round(n.cell^(0.9)), 
                                  min = ext["xmin"], 
                                  max = ext["xmax"]),
                   "lat" = runif(n = round(n.cell^(0.9)), 
                                 min = ext["ymin"], 
                                 max = ext["ymax"]))

# the suitability scores of those potential sites
pot.suits <- suit.rstr[cellFromXY(object = suit.rstr, 
                                  xy = pot.sites)]
# exclude some potential sites according to suitability, ending up with the following site locations
site.coords <- pot.sites[rbinom(n = length(pot.suits), 
                                size = 1, 
                                prob = pot.suits^4) == 1,]
plot(site.coords, xlim = ext[1:2], ylim = ext[3:4], 
     pch = 16, cex = 0.5, asp = 1, main = "Observed\n Site Locations")

# find what cells of our raster our sites fall into
site.rstr <- suit.rstr
values(site.rstr) <- 0
values(site.rstr)[cellFromXY(object = site.rstr, xy = site.coords)] <- 1
sites <- values(site.rstr)

#########################
#### plot everything ####
#########################

par(mfrow = c(2,3)) # tells the graphics device to make room to draw six plots in the same window arranged in 2 rows and 3 columns
plot(elev.rstr, col = pal(n.cell), asp = 1, main = "Elevation (m)")
plot(prec.rstr, col = pal(n.cell), asp = 1, main = "Annual\n Precipitation (cm)")
plot(popd.rstr, col = pal(n.cell), asp = 1, main = "Population Density")
plot(temp.rstr, col = pal(n.cell), asp = 1, main = "Average Daily Temperature\n (degrees C)")
plot(soil.rstr, col = pal(n.cell), asp = 1, main = "Soil Salinity (dS/m)")
plot(suit.rstr, col = pal(n.cell), asp = 1, main = "Site Suitability")

par(mfrow = c(1,1)) # go back to just plotting one thing in the plot window

##############################
#### GLM cross validation ####
##############################

# we now turn our spatial data into a single big matrix; each row represents a cell, each column, a variable
# make sure the rows of the `dat` object stay in the same order, otherwise we won't 
# be able to reassemble our predictions into a coherent heatmap of likely site locations
dat <- cbind("site" = as.numeric(sites),
             "elev" = as.numeric(elevation),
             "prec" = as.numeric(ann.precip),
             "popd" = as.numeric(pop.dens),
             "temp" = as.numeric(avg.temp),
             "soil" = as.numeric(soil.sal))

dat <- as.data.frame(dat) # convert dat from a matrix to a data.frame; the glm() function only works with data.frames

# we'll now perfolm a k-fold cross validation, k = 5 or 10 is standard.
# this is just to see how good the particular method we're using (General Linear 
# Modelling) is at making predictions from this particular dataset
k <- 10
folds <- k.fold(n.obs = n.cell, k = k)
thresholds <- vector("numeric", length = k)
mccs <- vector("numeric", length = k)
for(i in 1:k){
  training <- dat[folds != i,] # most of the data rows are used for training
  testing <- dat[folds == i,] # 1/k of them are used for testing
  
  mod.glm <- glm(site 
                 ~ elev 
                 + prec 
                 + popd 
                 + temp 
                 + soil, 
                 data = training, 
                 family = "binomial") # train a glm model using our training data
  
  prdxn <- predict(object = mod.glm, 
                   newdata = testing) # put the testing data into the model we just made and see what it predicts for site pres/abs
  prdxn <- logistic(prdxn) # the modelling method we used transforms its output, this is just the inverse transformation
  
  # predictions are between 0 and 1. we have to pick a number above which
  # predictions are counted as 1s and below which are counted as 0s. the 
  # opt.threshold() function finds the value of this threshold that makes 
  # our prediction get the closest to the observed site p/a data. these 
  # thresholds are stored and will be put through a weighted average to give 
  # us a single threshold value to use for all further predictions
  t <- opt.threshold(pa = testing[,"site"], prdxn = prdxn)
  prdxn[prdxn < t] <- 0
  prdxn[prdxn >= t] <- 2
  a <- cbind(testing[,"site"], prdxn)
  b <- rowSums(a)
  
  tp <- sum(b == 3) + 1 # true positives count # adding 1s helps avoid division by zero
  fp <- sum(b == 2) + 1 # false positives
  fn <- sum(b == 1) + 1 # false negatives
  tn <- sum(b == 0) + 1 # true negatives
  
  thresholds[i] <- t # store threshold that was used
  mccs[i] <- mcc(tp = tp, 
                 fp = fp, 
                 tn = tn, 
                 fn = fn) # store the matthews correlation coefficient (measures how good the prediction was)
}

######################################
#### build final prediction model ####
######################################

final.mod <- glm(site 
                 ~ elev 
                 + prec 
                 + popd 
                 + temp 
                 + soil, 
                 data = dat, 
                 family = "binomial") # train a new glm model using our whole dataset

# the "Estimate" column of this summary should *roughly* reflect the linear combination 
# of environmental variables we used to define our suitability score
summary(final.mod)

mean(mccs) # on average, this is how well the GLMs did at site prediction (0 = random chance, 1 = perfect)
# for any predictions made with this model, we'll use the following threshold to split predictions into 0s and 1s
final.threshold <- sum(thresholds * mccs / sum(mccs))

############################################
#### apply this model to a new location ####
############################################

# we assume that suitability is determined in this new location the same way it is in our old location.
# we get new values for the same variables we used before so that we can put them into the model we already trained

# we generate some NEW spatial environmental data
elevation.2 <- rand.surface(n.row = n.row, n.col = n.col) * 300 # elevation in meters
ann.precip.2 <- rand.surface(n.row = n.row, n.col = n.col) * 50 # average annual precipitation in cm
pop.dens.2 <- rand.surface(n.row = n.row, n.col = n.col) * 100 # population density in individuals / sq. km
avg.temp.2 <- rand.surface(n.row = n.row, n.col = n.col) * 40 # average daily temperature in degrees C
soil.sal.2 <- rand.surface(n.row = n.row, n.col = n.col) * 5 # soil salinity in deciSiemens / meter

# turn the above NEW objects from matrices into raster objects, making them plottable
elev.rstr.2 <- raster(elevation.2)
prec.rstr.2 <- raster(ann.precip.2)
popd.rstr.2 <- raster(pop.dens.2)
temp.rstr.2 <- raster(avg.temp.2)
soil.rstr.2 <- raster(soil.sal.2)

ext.2 <- c("xmin" = -88.0997,
           "xmax" = -84.7846,
           "ymin" = 37.7717,
           "ymax" = 41.7614) # the bounding box for Indiana

extent(elev.rstr.2) <- 
  extent(prec.rstr.2) <- 
  extent(popd.rstr.2) <- 
  extent(temp.rstr.2) <- 
  extent(soil.rstr.2) <- ext.2 # set the lat/lon region covered by our rasters

p4s <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0") # a string specifying a map projection

proj4string(elev.rstr.2) <- 
  proj4string(prec.rstr.2) <- 
  proj4string(popd.rstr.2) <- 
  proj4string(temp.rstr.2) <- 
  proj4string(soil.rstr.2) <- p4s # set the map projection to which coordinates are relative

##############################
#### plot the new rasters ####
##############################

par(mfrow = c(2,3))
plot(elev.rstr.2, col = pal(n.cell), asp = 1, main = "New Elevation (m)")
plot(prec.rstr.2, col = pal(n.cell), asp = 1, main = "New Annual\n Precipitation (cm)")
plot(popd.rstr.2, col = pal(n.cell), asp = 1, main = "New Population\n Density")
plot(temp.rstr.2, col = pal(n.cell), asp = 1, main = "New Average Daily\n Temperature (degrees C)")
plot(soil.rstr.2, col = pal(n.cell), asp = 1, main = "New Soil Salinity\n (dS/m)")

############################
#### run the prediction ####
############################

dat.2 <- cbind("elev" = as.numeric(elevation.2),
               "prec" = as.numeric(ann.precip.2),
               "popd" = as.numeric(pop.dens.2),
               "temp" = as.numeric(avg.temp.2),
               "soil" = as.numeric(soil.sal.2))

dat.2 <- as.data.frame(dat.2)

# we run our whole new dataset through the glm trained on the old data to 
# get our suitability predictions for the new, unexplored location
pred.suitability <- predict(object = final.mod, newdata = dat.2)
pred.suitability <- logistic(pred.suitability)
pred.suitability <- matrix(pred.suitability, nrow = n.row)

pdct.rstr <- raster(pred.suitability) # turn the prediction into a raster to match the others
extent(pdct.rstr) <- ext.2
proj4string(pdct.rstr) <- p4s

# plot the predicted suitability surface
par(mfrow = c(1,1))
plot(pdct.rstr, col = pal(n.cell), main = "Predicted Site\n Suitability")

# which cells from the raster we just plotted meet or exceed our model's suitability threshold?
# regions where black dots get plotted are where you should look for sites
suitable.cells <- which(values(pdct.rstr) > final.threshold) 
suitable.coords <- xyFromCell(object = pdct.rstr, cell = suitable.cells)
points(suitable.coords, pch = 16, cex = 0.4)
