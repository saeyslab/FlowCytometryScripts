# Load the preprocessed data:
# See script_preprocessing.R
#    ff:            Compensated flowFrame
#    ff_t:          Compensated and logicle transformed flowFrame
#    manual:        Array with label for each cell
#    selected:      Array with TRUE/FALSE whether cell falls in single live 
#                   cells
#    gatingMatrix:  Matrix with rows corresponding to cells and a column for 
#                   each manual gate. Each column contains TRUE/FALSE values
#                   indicating whether the cells fall in the specific gate
#    colsToCluster: Columns to use for clustering
load("FR-FCM-ZZQY/21-10-15_Tube_028.Rdata")

# Load the plot settings
# See script_plotSettings.R
#    cellTypeColors,markerColors,markerNames,
#    circular_markerOrder,grid_markerOrder
load("plotSettings.Rdata")

# Load the FlowSOM library
library(Rtsne)

# Set some parameters
tSNE_subsample <- 10000

# Set seed for reproducable results
set.seed(42)
# Record start time
start <- Sys.time()
    # Run the rtsne algorithm
    rtsne_res <- Rtsne(exprs(ff_t[selected,colsToCluster])[seq_len(tSNE_subsample),])
    # Plot the results, marker plots in separate png because pdf is too big
    
    # Plot manual
    pdf("Rtsne_manual.pdf",useDingbats = FALSE)
    plot(rtsne_res$Y,col=c("#888888",cellTypeColors)[manual[selected][seq_len(tSNE_subsample)]],pch=19,
         bty="n",axes=F,xlab="",ylab="")
    dev.off()
    
    # Plot individual markers
    png("Rtsne_markers.png",width = 1200,height=1200)
    par(mfrow=c(3,3))
    for(m in grid_markerOrder){
        channel <- colnames(ff)[m]
        plot(rtsne_res$Y,
             col=colorRampPalette(c("#dddddd",markerColors[channel]))(100)[
                 as.numeric(cut(exprs(ff_t[selected,channel])[seq_len(tSNE_subsample),],
                                breaks = 100))],
             main=markerNames[channel],bty="n",axes=F,xlab="",ylab="",pch=19)
    }
    par(mfrow=c(1,1))
    dev.off()
    
# Record end time
t_Rtsne_10000<- Sys.time() - start
# Save results
save(t_Rtsne_10000, rtsne_res, file="rtsne_10000.Rdata")