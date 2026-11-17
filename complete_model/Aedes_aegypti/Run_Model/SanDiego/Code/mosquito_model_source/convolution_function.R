## Convolution For Age Progession in Adult Mosquito Life Stages

### Description: ###
# creates function to convolve a randomly generated numerical pdf of age shift and current count vector of individuals in respective life stage

### Arguments: ###
# 1: x1 - vector of adults for each development stage
# 2: y1 - density of exponentially distributed adult development rate
# 3: padsize - extends vector size to account for individuals surpassing 100% development

### Final_Output: ###
# conv: new population of individuals at each development bin
# emigrants: individuals leaving life stage

# This function performs a convolution of two variables of the same length
# padding on the is included on the right hand side
ConvolveFunc = function(x1, y1, padsize){

  x2 = c(x1, rep(0.0, padsize))
  y2 = c(y1, rep(0.0, padsize))
  
  # fast Fourier transforming
  fx2 = fft(x2)
  fy2 = fft(y2)
  
  # Perform convolution and inverse transforming
  Convolution1 = Re(fft(fx2*fy2, inverse = T)/length(x2))
  
  # Now clipping
  Convolution2 = Convolution1[1:length(x1)]
  
  # compute how many have moved into the padded region
  Emigrants = sum(Convolution1[(length(x1)+1):length(x2)])
  
  return(list(conv = Convolution2, emigrants = Emigrants))
}
