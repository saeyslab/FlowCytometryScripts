library(flowCore)

removeMargins <- function(flowFrame,dimensions){
    # Look op the accepted ranges for the dimensions
    meta <- pData(flowFrame@parameters)
    rownames(meta) <- meta[,"name"]
    
    # Initialize variables
    selection <- rep(TRUE,times=nrow(flowFrame))
    e <- exprs(flowFrame)
    
    # Make selection
    for(d in dimensions){
        selection <- selection & 
            e[,d] > max(meta[d,"minRange"],min(e[,d])) &
            e[,d] < min(meta[d,"maxRange"],max(e[,d]))
    }
    return(selection)
}

selectGoodQuality <- function(dir,file,comp,toTransform,
                              checkMargins,
                              markersToPlot,
                              binTime = 0.1,
                              useCountLower = TRUE, useCountUpper= TRUE, countMax=600, countThreshold = 100,
                              timeMin = NULL, timeMax=NULL,
                              removeMargins = TRUE,
                              writeOutput = TRUE, outputDir = "QC",
                              showPlot = TRUE){
    # Read the file
    ff <- read.FCS(file.path(dir,file)) 
    ff <- compensate(ff,comp)
    colnames(ff)[which(colnames(ff) %in% colnames(comp))] <- paste0("Comp-",colnames(comp))
    
    # Change time to seconds
    toSec <- function(str){
        str <- as.numeric(strsplit(str,':')[[1]])
        60*60*str[1]+60*str[2]+str[3]
    }
    totalTime <- 
        toSec(ff@description$`$ETIM`)-toSec(ff@description$`$BTIM`)
    time <- 
        matrix(exprs(ff)[,"Time"]*totalTime/max(exprs(ff)[,"Time"]),
               ncol=1,dimnames = list(NULL,"Time"))
    exprs(ff)[,"Time"] <- time
    
    ff_t <- transform(ff,transformList(colnames(ff[,toTransform]),logicleTransform()))
    
    # Split the data in bins
    nIntervals <- round(totalTime/binTime)
    cuts <- cut(exprs(ff_t)[,"Time"],breaks=nIntervals)
    splitted <- split(as.data.frame(exprs(ff_t)),cuts)
    timePoints <- as.numeric(gsub("\\(","",gsub(",.*","",names(splitted))))
    
    # Logical vector to indicate good cells
    selected <- rep(TRUE,nrow(ff))
    
    counts <- sapply(splitted,function(flowDataFrame){dim(flowDataFrame)[1]})
    if(showPlot){
        plot(cbind(timePoints, counts),ylim=c(0,countMax))
    }
    if(max(counts)>countMax){warning(paste0("Counts in ",file," higher than ",countMax,"!"))}
    
    # Counts over time
    if(useCountLower | useCountUpper){
        sampleMedian <- median(counts)
        if(showPlot){
            if(useCountLower) abline(h=sampleMedian-countThreshold,col="#FF0000",lwd=2)
            if(useCountUpper) abline(h=sampleMedian+countThreshold,col="#FF0000",lwd=2)
        }
        selectedBins <- names(which((counts > sampleMedian-countThreshold) & 
                                        (counts < sampleMedian+countThreshold)))
        selected <- selected & (cuts %in% selectedBins)
        message(paste0(table(selected)["FALSE"]," cells removed due to counts"))
    } 
    
    if(!is.null(timeMin)){
        if(showPlot) abline(v=timeMin,col="#FF0000",lwd=2)
        selected <- selected & (ff[,"Time"] > timeMin)
        message(paste0(table((ff[,"Time"] > timeMin))["FALSE"]," cells removed due to time min"))
    }
    
    if(!is.null(timeMax)){
        if(showPlot) abline(v=timeMax,col="#FF0000",lwd=2)
        selected <- selected & (ff[,"Time"] < timeMax)
        message(paste0(table((ff[,"Time"] < timeMax))["FALSE"]," cells removed due to time max"))
    }
    
    toRemove <- removeMargins(ff,checkMargins)
    selected <- selected & toRemove
    message(paste0(table(toRemove)["FALSE"]," cells removed due to removeMargins"))
    
    ff_new <- cbind2(ff,matrix(as.numeric(selected),dimnames = list(NULL,"good")))
    if(writeOutput){
        suppressWarnings(dir.create(file.path(dir,outputDir)))
        #write.table(selected,file=file.path(dir,"QC",paste0("QC_",file,".csv")),row.names = FALSE,col.names = c("Selected"))
        
        suppressWarnings(write.FCS(ff_new, file.path(dir, outputDir, gsub(".fcs","_good.fcs",file))))
    }
    
    if(showPlot){
        
        markers <- ff@parameters@data[,"desc"][markersToPlot]
        markers[is.na(markers)] <- ff@parameters@data[,"name"][markersToPlot[is.na(markers)]]
        names(markers) <- colnames(ff)[markersToPlot]
        
        # Plot marker intensities over Time
        for(marker in c(names(markers))){
            medians <- sapply(splitted,function(flowDataFrame){
                median(flowDataFrame[,marker])
            })
            qtl_25 <- sapply(splitted,function(flowDataFrame){
                quantile(flowDataFrame[,marker],c(.25))
            })
            qtl_75 <- sapply(splitted,function(flowDataFrame){
                quantile(flowDataFrame[,marker],c(.75))
            })
            if(marker=="FSC-A" || marker=="SSC-A"){
                plot(exprs(ff_t)[,c("Time",marker)],
                    pch=".",col=c("#FF0000","#00000055")[selected+1],
                    #ylim=c(min,max),#axes=FALSE,
                    main=paste0(markers[marker]," (",gsub("Comp-","",gsub("-A$","",marker)),")"))
            } else {
                plot(exprs(ff_t)[,c("Time",marker)],
                     pch=".",col=c("#FF0000","#00000055")[selected+1],
                     ylim=c(-2,4),#axes=FALSE,
                     main=paste0(markers[marker]," (",gsub("Comp-","",gsub("-A$","",marker)),")"))
            }
            points(cbind(timePoints,medians),
                   pch="o",col="#3498db")
            points(cbind(timePoints,qtl_25),
                   pch="o",col="#a6bddb")
            points(cbind(timePoints,qtl_75),
                   pch="o",col="#a6bddb")
        }
    }
    
    return(ff_new)
}

dir  <- "FR-FCM-ZZQY"
files <- list.files(dir,pattern="Tube_[0-9]*.fcs$")
dir.create(file.path(dir,"QC"))
    

compensationFile = "FR-FCM-ZZQY/attachments/CompensationFlowJo.csv"
comp <- read.csv(compensationFile,row.names=1,check.names = FALSE)
colnames(comp) <- rownames(comp) <- gsub(" ::.*","",colnames(comp))

for(file in files){
    print(file)
    png(file.path(dir,"QC",paste0("QC_1s_",file,".png")),width=2000,height=1000)
        
    layout(matrix(c(1,1,1,3,4,5,
                    1,1,1,3,4,5,
                    1,1,1,6,7,8,
                    2,2,2,6,7,8,
                    2,2,2,9,10,11,
                    2,2,2,9,10,11), 6, 6, byrow = TRUE))
    ff <- selectGoodQuality(dir=dir,file=file,comp=comp,toTransform=7:19,
                            countMax = 10000, countThreshold = 2000,
                            checkMargins = c(1:6),binTime = 1,
                            markersToPlot = c(4,1,10,11,12,14,15,17,18,19))
    dev.off()
}