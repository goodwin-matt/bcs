# Import the "signal" library to plot spectogram. 
library(signal)
# Import the "pracma" library to use the function "findpeaks".
library(pracma)
# This function is a modified version of the function 'specgram' from the library 'signal'
FindSpecgram <- function (x, n = min(256, length(x)), Nfft=0, Fs = 2, window = hanning(n), 
                          overlap = ceiling(length(window)/2), real=TRUE)
{
  # Finds the spectrogram of a given vector x. This is a modification of the function 
  # 'specgram' from the library 'signal'.
  #
  # Args:
  #   x: the vector of samples.
  #   n: the size of the Fourier transform window.
  #   Nfft: this is an added argument from the original library. 
  #         Represents how many DFT coefficients to return.
  #   Fs: the sample rate, Hz.
  #   window: shape of the fourier transform window, defaults to hanning(n).
  #   overlap: overlap with previous window, defaults to half the window length.
  #   real: whether x is real or not. Defaults to TRUE
  #
  # Returns:
  #   The matrix S, representing the spectrogram of x. Also the vectors f and t
  #   representing the frequency and time vectors of the spectrogram.
  if (!is.numeric(x)) 
    stop("'x' has to be a numeric.")
  if (length(window) == 1) 
    window <- hanning(window)
  if (length(n) > 1) 
    stop("specgram does not handle frequency vectors yet")
  win_size <- length(window)
  if (win_size > n) {
    n <- win_size
    warning("specgram fft size adjusted to", n)
  }
  step <- win_size - overlap
  if (length(x) > win_size) 
    offset <- seq(1, length(x) - win_size, by = step)
  else offset <- 1
  S <- matrix(0, n, length(offset))
  for (i in seq_along(offset)) {
    S[1:win_size, i] <- x[offset[i]:(offset[i] + win_size - 
                                       1)] * window
  }
  
  # This part has been added in from the original 'signal' library. The idea is 
  # to 'pad' the matrix S with zeros before calculating the fft of S.
  if (Nfft>dim(S)[1]) {
    added.zeros <- (Nfft-n)*length(offset)
    padded.zeros <- matrix(rep(0,added.zeros),ncol=length(offset),nrow=(Nfft-n))
    S <- rbind(S,padded.zeros)
  }
  
  S <- mvfft(S)
  
  # Readjust the dimensions. Also added in from the original 'signal' library.
  if (Nfft > n & real) {
    if (Nfft%%2 == 1)
      ret_n <- (Nfft+1)/2
    else
      ret_n <- Nfft/2 + 1
  }
  else
    ret_n <- n/2
  S <- S[1:ret_n, ]
  f <- (0:(ret_n - 1)) * Fs/n
  t <- offset/Fs
  res <- list(S = S, f = f, t = t)
  class(res) <- "specgram"
  return(res)
}

# Make a chirp signal
Fs <- 8000  # sampling frequency    						
Ts <- 1/Fs  # sample period							
Tfinal <- 10  # sumber of seconds of data
tlist <- seq(0,Tfinal,by=Ts)  # time list
Nt <- length(tlist)
N_eta <- 15  # number of harmonics
alist <- rep(1,N_eta)  # signal amplitudes
sigma <- .001  # noise variance

# Parameters of the fundamental frequency
toff <- -1.5  # shift of tanh
tfineff <- 1  # final tanh time
Fmin <- 480  # limit inital frequency
Fmax <- 2100  # limit final frequency

# Create fundamental frequencey
tau <- (tfineff - toff)/Tfinal * tlist + toff  # time list
flist <- Fmin + (Fmax-Fmin) * (tanh(tau)+1)/2  # frequency list
phit <- 2*pi*cumsum(flist)*Ts  # phase (integral of freq.)
plot(tlist, flist, type='l') # plot the fundamental frequency 

# Generate signal and all harmonics
j <- complex(real=0, imaginary=1.0)
s <- rep(0,Nt)

for (i in 1:N_eta){
  s <- s + alist[i]*cos(i*phit)
}
if (sigma) {  # add noise
  s <- s #+ sigma*rnorm(Nt)
}

# # Compute and plot the spectrogram
# windowlen <- 500
# noverlap <- windowlen/2
# sampnew <- windowlen - noverlap  # number of new samples each time
# Tchunk <- sampnew * Ts  # time increment for each chunk
# NFFT <- 512  # this argument determines how many DFT coefficients to find
# spec <- FindSpecgram(s, n=windowlen, Nfft=NFFT, Fs=Fs, overlap=noverlap, 
#                      window=hamming(windowlen))
# signal:::plot.specgram(spec,col=heat.colors(8000))
# 
# # Find data for analysis
# S <- spec$S
# f <- spec$f
# time <- spec$t
# 
# # step through the rows of S
# ns <- dim(S)[1]  # number of rows in spectogram, ns = NFFT/2 + 1 for real data
# ms <- dim(S)[2]  # number of columns in spectogram, ms = number of time slots
# 
# flist <- seq(0,Fs/2,length=ns)
# sizethresh <- 1
# MINPEAKHEIGHT <- 10
# t0 <- windowlen/2 * Ts  # time of middle of data
# nu <- 0  # initial frequency change rate
# 
# for (k in 1:ms) {  # 1:ms for each chunk of data
#   t <- t0 + (k-1)*Tchunk
#   cat(sprintf("k=%s t=%g\n", k, t))
#   # plot(c(t, t),c(0, Fs/2),type='l')  # plot the current time step
#   
#   plot(flist, abs(S[,k]), type='l', col='blue')
#   findpeaks.output <- findpeaks(abs(S[,k]))
#   pks <- findpeaks.output[,1]
#   locs <- findpeaks.output[,2]
#   points(flist[locs], pks, col='red', pch=18)
#   
#   npeaks <- length(pks)
#   px <- c()  # the peak locations we keep
#   p <- c()  # the peak values we keep
#   keeplist <- c()
#   i1 <- 0
#   for (i2 in 1:npeaks) {
#     if (pks[i2] > MINPEAKHEIGHT) {
#       keeplist <- c(keeplist, i2)
#       i1 <- i1+1
#       px[i1] <- locs[i2]
#       p[i1] <- pks[i2]
#     }
#   }
#   npeaksfound <- i1
#   points(flist[px], p, cex=1.5)
#   cat(sprintf("npeaks=%d npeaksfound=%d\n", npeaks, npeaksfound))
# }  # for k (for each frequency chunk)
