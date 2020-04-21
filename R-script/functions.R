# R functions to estimate and apply the variable transformation
# Remove missing values before running the transformation
# (Section 2 in Guidance document; DOI: 10.13140/RG.2.2.10494.95040)

estimateBC = function(x){ 
  # function to estimate transformation parameters for continuous variable x 
  require(car) 
  gamma = min(x, na.rm=T) - 0.001 # offset (min value minus a small number) 
  x = x - gamma # subtract gamma from x, so that it is strictly positive 
  lambda = powerTransform(x~1, family="bcPower")$lambda # estimate lambda of Box-Cox transformation... 
  xT = bcPower(x, lambda=lambda) # apply box-cox transform 
  xT.mean = mean(xT) # mean of transformed values, for centring 
  xT.sd = sd(xT) # sd of transformed values, for scaling 
  # return the transformation parameters 
  return(c(gamma=gamma, lambda=lambda, xT.mean=xT.mean, xT.sd=xT.sd)) 
}

applyBC = function(x, P=estimateBC(x)){ 
  # function to transform continuous variable x using transformation parameters P 
  require(car) 
  gamma = P[1] 
  lambda = P[2] 
  xT.mean = P[3] 
  xT.sd = P[4] 
  xT = bcPower(x-gamma, lambda) # apply box-cox transform 
  xT = (xT-xT.mean)/xT.sd # centre and scale 
  return(xT) 
}

backBC = function(xT, P){ # function to back transform transformed variable xT using transformation parameters P 
  gamma=P[1] 
  lambda=P[2] 
  xT.mean=P[3] 
  xT.sd=P[4]
  xT.unscaled = xT*xT.sd + xT.mean 
  x.original = exp(log(lambda*xT.unscaled + 1)/lambda) + gamma 
  return(x.original) 
}
  
  
  
  
  
  