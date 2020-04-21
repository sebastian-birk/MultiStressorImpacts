source("functions.R")

# Transformation of continuous variable (y) and continuous stressors variables (x1 and x2)----

library(car)

P.y = estimateBC(y)
yT = applyBC(y, P.y)

P.x1 = estimateBC(x1)
x1T = applyBC(x1, P.x1)

P.x2 = estimateBC(x2)
x2T = applyBC(x2, P.x2)

# Hint: categorical response variables have to be an ordered factor

################### Select type of model #######################################
# depends on the type of response variable
# decide if you need a random effect, decision help in table 3 in the
#  MARS WP6 Synthesis guidance (DOI: 10.13140/RG.2.2.10494.95040)

## Legend
# M = fitted model
# y = response variable
# yT = transformed values of the response variable
# x1T, x2T = transformed values of the stressor variables
# RE1, RE2 = random effects (grouping factors, e.g. block, site, year)


## continous response variable

### linear model (no random effects)
### generalised linear model (GLM)
M = glm(yT ~ x1T*x2T)
### mixed effects model (with random effects)
### generalised linear mixed model (GLMM)
library(lmerTest)
M = lmer(yT ~ x1T*x2T + (1|RE1) + (1|RE2), REML = FALSE)

## binary response variable

### linear model (no random effects)
### generalised linear model (GLM)
M = glm(yT ~ x1T*x2T, family = "binomial")

### mixed effects model (with random effects)
### generalised linear mixed model (GLMM)
library(lme4)
M = glmer(yT ~ x1T*x2T + (1|RE1) + (1|RE2), family = "binomial")

## count (response variable)

### linear model (no random effects)
### generalised linear model (GLM)
M = glm(yT ~ x1T*x2T, family = "poisson")

### mixed effects model (with random effects)
### generalised linear mixed model (GLMM)
library(lme4)
M = glmer(yT ~ x1T*x2T + (1|RE1) + (1|RE2), family = "poisson")

## ordered categories (response variable)

### linear model (no random effects)
### cumulative link model (CLM)
library(ordinal)
M = clm(y ~ x1T*x2T)

### mixed effects model (with random effects)
### cumulative link mixed model (CLMM)
library(ordinal)
M = clmm(y ~ x1T*x2T + (1|RE1) + (1|RE2))

################### Test and correct for residual autocorrelation ##############

# Test only necessary for linear models without random effects
# see MARS WP6 guidance to check, if you need to test for residual autocorrelation

# 1st step: extract residuals from the fitted model
# for consistency across models, we will use response residuals
# response residual = observed response minus the fitted response

# calculation of response residuals r depends on the type of response variable

## continous response variable
r = residuals(M, type = "response")

## binary response variable
r = residuals(M, type = "response")

## count (response variable)
r = residuals(M, type = "response")

## ordered categories (response variable)
### matrix of fitted probabilities for each category
classProbs = predict(M, newdata=data.frame(x1T,x2T))$fit 
### binary dummy matrix of observed categories
classObs = sapply(levels(y), function(x) { as.numeric(x==y) }) 
### response residuals for each category
r = classObs - classProbs 

# 2nd step: test for temporal and/or spatial autocorrelation

## temporal autocorrelation

### create inverse weights from the differences in years between each pair of observations
w = 1/as.matrix(dist(year))
diag(w) = 0
### perform Moran’s test on the residuals (r)
library(ape)
Moran.I(x=r, weight=w)

## spatial autocrrelation 

### calculate the distances between each pair of sites, using their decimal 
###  longitudes and latitudes
library(sp)
d = spDists(cbind(lon,lat), longlat=T)
### convert the distance matrix into inverse weights
w = 1/d
diag(w) = 0
### perform Moran’s test on the residuals (r)
library(ape)
Moran.I(x=r, weight=w)

# We suggest: substantive problems will be indicated by Moran´s I > 0.1

# 3rd step (a): implement the temporal trend surface linear models

## continous response variable
### estimate the temporal trend of the linear model
library(mgcv)
G = gam(yT ~ x1T*x2T + s(year))
trendSurface = predict(G, 
                       newdata = data.frame(x1T = 0, x2T = 0, year))-coef(G)[1]
### update the linear model with the temporal trend surface
M = update(M, offset = trendSurface)

## binary response variable
### estimate the temporal trend of the linear model
library(mgcv)
G = gam(yT ~ x1T*x2T + s(year),family = "binomial")
trendSurface = predict(G, 
                       newdata = data.frame(x1T = 0, x2T = 0, year))-coef(G)[1]
### update the linear model with the temporal trend surface
M = update(M, offset = trendSurface)

## count (response variable)
### estimate the temporal trend of the linear model
library(mgcv)
G = gam(yT ~ x1T*x2T + s(year),family = "poisson")
trendSurface = predict(G, 
                       newdata = data.frame(x1T = 0, x2T = 0, year))-coef(G)[1]
### update the linear model with the temporal trend surface
M = update(M, offset = trendSurface)

## ordered categories (response variable)
### estimate the temporal trend of the linear model
#### NOT NEEDED
### update the linear model with the temporal trend surface
M = clm(y ~ x1T*x2T + poly(year,3))


# 3rd step (b): implement the spatial trend surface linear models

## continous response variable
### estimate the spatial trend of the linear model
library(mgcv)
G = gam(yT ~ x1T*x2T + s(lon,lat))
trendSurface = predict(G, 
                       newdata = data.frame(x1T = 0, x2T = 0, long,lat))-coef(G)[1]
### update the linear model with the spatial trend surface
M = update(M, offset = trendSurface)

## binary response variable
### estimate the spatial trend of the linear model
library(mgcv)
G = gam(yT ~ x1T*x2T + s(lon, lat),family = "binomial")
trendSurface = predict(G, 
                       newdata = data.frame(x1T = 0, x2T = 0, lon, lat))-coef(G)[1]
### update the linear model with the spatial trend surface
M = update(M, offset = trendSurface)

## count (response variable)
### estimate the spatial trend of the linear model
library(mgcv)
G = gam(yT ~ x1T*x2T + s(lon, lat),family = "poisson")
trendSurface = predict(G, 
                       newdata = data.frame(x1T = 0, x2T = 0, lon, lat))-coef(G)[1]
### update the linear model with the spatial trend surface
M = update(M, offset = trendSurface)

## ordered categories (response variable)
### estimate the spatial trend of the linear model
#### NOT NEEDED
### update the linear model with the spatial trend surface
M = clm(y ~ x1T*x2T + poly(lon, 3)*poly(lat,3))

################### Model evaluation ###########################################

# 1. The marginal and conditional model R².
## marginal R² uses only the fixed effects
## conditional R² also includes the random effects (not possible for categorical data)
library(MuMIn)
r.squaredGLMM(M)

# 2. Pearson´s product-moment correlation between the model response residuals (r) 
# and the fitted values
## code for extracting residuals from most model types in autocorrelation sections

## for GLMM: use
r = residuals(M, type = "response")

## for the linear mixed model of ordered categorical responses (CLMM):
### residuals can not be extracted so no test can be performded.
### to test the correlation:
cor.test(r, fitted(M))

## for the linear model of ordered categorical responses (CLM) 
## where there are residuals for each category calculated as in the tables above, use:
for(i in 1:ncol(r)) print(cor.test(r[,i], classProbs[,i]))

# 3. The Shapiro-Wilk test for normality of residuals (r)
## Remember that models for categorical responses have residuals for each
## different category and require separate testing. 
## The test should be non-significant, unless the dataset is very large
## in which case very small deviations from normality are statistically significant.
shapiro.test(r)

## For the linear model of ordered categorical responses (CLM)
## where there are residuals for each category, use:
apply(r, 1, shapiro.test)

## For the linear mixed model of ordered categorical responses (CLMM),
## residuals cannot be extracted so no test can be performed.


################### Determine importance of the interaction term ###############

# Three simple methods to estimate the importance of the interaction term will be used:

## 1. The Z-score for the interaction term
### large absolute values indicate an important interaction
(summary(M)$coef[,1]/summary(M)$coef[,2])["x1T:x2T"]

## 2. The change in model AIC if the interaction is removed
### large positive values indicate an important interaction
AIC(update(M, ~.-x1T:x2T)) - AIC(M)

## 3. The drop in marginal and conditional R² if the interaction is removed
### large negative values indicate an important interaction
r.squaredGLMM(update(M, ~.-x1T:x2T)) - r.squaredGLMM(M)

## Note that R² cannot be calculated for models of categorical responses.


################### Visualise the interaction ##################################

# The fitted model equation can be used to plot the response surface
# for the two main stressors, x1T and x2T.
## First, you will need to extract the model fixed effect coefficients:

### model function used: glm
### fitted fixed effect for model object M
B = coef(M)

### model function used: lmer
### fitted fixed effect for model object M
B = fixef(M)

### model function used: glmer
### fitted fixed effect for model object M
B = fixef(M)

### model function used: clm
### fitted fixed effect for model object M
#### This code plots the transition between the first two ordered categories
#### of the response variable. Change [1] to [2],[3], etc.
#### to plot subsequent transitions
B = c(coef(M)[1], coef(M)[c("x1T", "x2T", "x1T:x2T")])

### model function used: clmm
### fitted fixed effect for model object M
#### This code plots the transition between the first two ordered categories
#### of the response variable. Change [1] to [2],[3], etc.
#### to plot subsequent transitions
B = c(coef(M)[1], coef(M)[c("x1T", "x2T", "x1T:x2T")])

# plot the response from a fitted mit with our function interactionPlot()
## description below the function

interactionPlot = function(B, X1, X2, Y,
                           TP=list(P.x1=if(exists("P.x1")) P.x1 else NA,
                                   P.x2=if(exists("P.x2")) P.x2 else NA,
                                   P.y=if(exists("P.y")) P.y else NA,
                                   family="gaussian"),
                           responseLab="z", x1Lab="x1", x2Lab="x2") {
  # Function to plot interactions from fitted models.
  # B = vector of model fixed effect coeffcients
  # X1, X2 = vectors with values of the stressors (not transformed)
  # Y = vector with values of the response (not transformed)
  # TP = list of transformation parameters containing the following elements;
  # P.x1 = output from estimateBC() for x1, or NA if no transformation applied
  # P.x2 = output from estimateBC() for x2, or NA if no transformation applied
  # P.y = output from estimateBC() for y, or NA if no transformation applied
  # family = family of the generalised model. It should be one of
  # "gaussian" (for continuous response),
  # "poisson" (for count response),
  # "binomial" (for binary or ordered categorical response), or
  # NA (if you want to plot the response on the linear model scale)
  # responseLab = label for the response variables
  # x1Lab, x2Lab = labels for the stressor variables
  
  require(emdbook)
  if(is.numeric(X1) & is.numeric(X2)){ # X1 and X2 are both continuous
    myF <<- function(X1=X1, X2=X2, transPar=TP) {
      if(sum(is.na(transPar$P.x1))==0) X1 = applyBC(X1, transPar$P.x1)
      if(sum(is.na(transPar$P.x2))==0) X2 = applyBC(X2, transPar$P.x2)
      z = B[1] + B[2]*X1 + B[3]*X2 + B[4]*X1*X2
      if(!is.na(transPar$family)){
        if(transPar$family=="gaussian") z = backBC(z, transPar$P.y)
        if(transPar$family=="poisson") z = exp(z)
        if(transPar$family=="binomial") z = 1/(1+exp(-z))
      }
      return(z)
    }
    curve3d(myF, from=c(min(X1),min(X2)), to=c(max(X1),max(X2)),
            sys3d="image", col=colorRampPalette(c("blue","white","red"))(64),
            xlab=x1Lab, ylab=x2Lab, main=responseLab, varnames=c("X1","X2"))
    points(X1, X2, pch=21, bg="grey50")
    curve3d(myF, from=c(min(X1),min(X2)), to=c(max(X1),max(X2)), sys3d="contour",
            add=T, labcex=1, varnames=c("X1","X2"))
    box()
  }
  
  if(is.numeric(X1) & is.factor(X2)){ # X1 continuous, X2 is a factor
    myF <<- function(X1, X2=0, transPar=TP) {
      if(sum(is.na(TP$P.x1))==0) X1 = applyBC(X1, TP$P.x1)
      z = B[1] + B[2]*X1 + B[3]*X2 + B[4]*X1*X2
      if(!is.na(TP$family)){
        if(TP$family=="gaussian") z = backBC(z, TP$P.y)
        if(TP$family=="poisson") z = exp(z)
        if(TP$family=="binomial") z = 1/(1+exp(-z))
      }
      return(z)
    }
    curve(myF(X1=x, X2=0), from=min(X1), to=max(X1), col="blue", lwd=2, xlab=x1Lab, ylab=responseLab,
          ylim=if(is.numeric(y) & !is.na(TP$family)) range(y) else NULL)
    curve(myF(X1=x, X2=1), add=T, col="red", lwd=2, lty=2)
    rug(x1)
    legend("topright", legend=levels(X2), lty=1:2, lwd=2, col=c("blue","red"), title=x2Lab)
  }
  
  if(is.factor(X1) & is.numeric(X2)){ # X1 is a factor, X2 continuous
    myF <<- function(X1=0, X2, transPar=TP) {
      if(sum(is.na(TP$P.x2))==0) X2 = applyBC(X2, TP$P.x2)
      z = B[1] + B[2]*X1 + B[3]*X2 + B[4]*X1*X2
      if(!is.na(TP$family)){
        if(TP$family=="gaussian") z = backBC(z, TP$P.y)
        if(TP$family=="poisson") z = exp(z)
        if(TP$family=="binomial") z = 1/(1+exp(-z))
      }
      return(z)
    }
    curve(myF(X1=0, X2=x), from=min(X2), to=max(X2), col="blue", lwd=2, xlab=x2Lab,
          ylab=responseLab, ylim=if(is.numeric(y) & !is.na(TP$family)) range(y) else NULL)
    curve(myF(X1=1, X2=x), add=T, col="red", lwd=2, lty=2)
    rug(x2)
    legend("topright", legend=levels(X1), lty=1:2, lwd=2, col=c("blue","red"),
           title=x1Lab)
  }
  
  if(is.factor(X1) & is.factor(X2)){ # X1 is a factor, X2 is a factor
    z = c(B[1], B[1]+B[2], B[1]+B[3], B[1]+B[2]+B[3], sum(B))
    names(z) = c(paste(levels(X1)[1], levels(X2)[1], sep=" / "),
                 paste(levels(X1)[2], levels(X2)[1], sep=" / "),
                 paste(levels(X1)[1], levels(X2)[2], sep=" / "),
                 paste("E(", paste(levels(X1)[2], levels(X2)[2], sep=" + "), ")", sep=""),
                 paste(levels(X1)[2], levels(X2)[2], sep=" x "))
    #names(z) = paste(rep(levels(X1),2), rep(levels(X2),each=2), sep=" / ")
    if(!is.na(TP$family)){
      if(TP$family=="gaussian") z = backBC(z, TP$P.y)
      if(TP$family=="poisson") z = exp(z)
      if(TP$family=="binomial") z = 1/(1+exp(-z))
    }
    z2 = z - z[1] # z values minus the control
    par(mfrow=c(2,1))
    barplot(z, ylab=responseLab, xlab=paste(x1Lab,x2Lab,sep=" / "))
    abline(h=0)
    barplot(z2, ylab=paste(responseLab,"- control"), xlab=paste(x1Lab,x2Lab,sep=" / "))
    abline(h=0)
    par(mfrow=c(1,1))
  }
  
  return(NULL)
}

# To use interactionPlot()

## you first need to create a list that holds the variable transformation
## parameters: outputs of estimateBC() if continous variables or NA for
## factors where no transformation is applied

## and the family model: this determines how to transform the predicted
## response variables ('gaussian', 'poisson', 'binomial') to plot
## responses on the scale of the raw data or NA to plot responses on the
## scale of the model linear function.

# For example, if both stressors and response are continuous,
# so a Gaussian model was used:
myTP = list(P.x1=P.x1, P.x2=P.x2, P.y=P.y, family="gaussian")

print(myTP)

library(emdbook)

interactionPlot(B=B, X1=x1, X2=x2, Y=y, TP=myTP, 
                responseLab="z (actual values)", 
                x1Lab="x1", x2Lab="x2")

# To make the equivalent plot on the linear scale of the model:
myTP = list(P.x1=P.x1, P.x2=P.x2, P.y=P.y, family=NA)

interactionPlot(B=B, X1=x1, X2=x2, Y=y, TP=myTP, 
                responseLab="z (model scale)", 
                x1Lab="x1", x2Lab="x2")


# For more information check the guidance document:
# DOI: 10.13140/RG.2.2.10494.95040
























