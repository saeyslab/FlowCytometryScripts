# Use same files as generated for FLOCK

# Read results generated with python
res_hdpgmm_10000 <- read.csv("hdpgmm/028_10k.csv",header = FALSE)
res_hdpgmm <- read.csv("hdpgmm/028_all.csv",header = FALSE)

# Plot on Rtsne mapping
load("rtsne.Rdata")
plot(rtsne_res$Y,col=rainbow(max(res_hdpgmm_10000))[res_hdpgmm_10000[,1]])
plot(rtsne_res$Y,col=rainbow(max(res_hdpgmm))[res_hdpgmm[1:10000,]])
