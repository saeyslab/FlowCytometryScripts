load("FR-FCM-ZZQY/21-10-15_Tube_028.Rdata")
write.table(exprs(ff_t[selected,colsToCluster])*100,sep = "\t",
            file="flock/21-10-15_Tube_028.txt",
            row.names = FALSE, quote = FALSE)

write.table(exprs(ff_t[selected,colsToCluster][1:10000,])*100,sep = "\t",
            file="flock/21-10-15_Tube_028_10000.txt",
            row.names = FALSE, quote = FALSE)

# Run on command line
# cd "C:\Users\svgassen\Documents\Flow Cytometry\data\NatureReview\flock"
# Measure-Command{.\flock2.exe .\21-10-15_Tube_028_10000.txt}
# Measure-Command{.\flock2.exe .\21-10-15_Tube_028.txt}

# File is always called flock_results.txt, manually renamed
res_flock_10000 <- read.table("flock/flock_results_10000.txt",header = TRUE,sep = "\t")[,"Population"]
res_flock <- read.table("flock/flock_results_all.txt",header = TRUE,sep = "\t")[,"Population"]

load("rtsne.Rdata")
plot(rtsne_res$Y,col=rainbow(max(res_flock_10000))[res_flock_10000])
plot(rtsne_res$Y,col=rainbow(max(res_flock))[res_flock[1:10000]])
