library(devtools)
document()
load_all()

path <- "./inst/fastdm_30_2"
asc <- Vectorize(function(x) strtoi(charToRaw(x),16L))
chr <- Vectorize(function(n) rawToChar(as.raw(n)))

#generate example data (rwiener from RWiener package)
IDs <- chr(65:74)
dat1 <- list()
for(i in IDs){
  rts <- RWiener::rwiener(1000,2,.25,.5,2) 
  rts$q <- round(rts$q, 3)
  rts$correct <- ifelse(rts$resp =="upper",1,0)
  names(rts) <- c("TIME","resp","RESPONSE")
  dat1[[i]] <- rts
  rm(rts)
}

#analyse data to recover parameters
out1 <- fastdm(~ zr*.5 + szr*0 + sv*0 + st0*0 + d*0, dat1, wd_fastdm = path, method = "ks", precision = 3, res_prefix = "ID")
#run() #use compiled library

print(out1)

#generate example data (construct-samples from fast-dm)
dat2 <- rdm(100, a=1, zr=.5, v=1.25, t0=0, N = 10) 

#analyse data to recover parameters
out2 <- fastdm(~ zr*.5 + szr*0 + sv*0 + st0*0 + d*0, dat2, method = "ks", precision = 3, res_prefix = "ID")
print(out1)

#visualize distribution function
predict.dm(type="cdf", output="plot", a=mean(out$a), v=mean(out$v), t0=mean(out$t0))