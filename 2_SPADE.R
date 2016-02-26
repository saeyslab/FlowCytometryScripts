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

# Load the SPADE library
library(spade)

# Set some parameters
SPADE_nClusters = 50

# Set seed for reproducable results
set.seed(42)
# Record start time
start <- Sys.time()
# Run the SPADE algorithm
out_dir <- "SPADE"
SPADE.driver("fcs/21-10-15_Tube_028_selected.fcs",
             out_dir=out_dir,
             cluster_cols=colnames(ff)[colsToCluster], 
             comp=FALSE, transforms=logicleTransform(),
             k=SPADE_nClusters)

# Plot the results
mst_graph <- igraph:::read.graph(paste(out_dir,"mst.gml",
                                       sep=.Platform$file.sep),
                                 format="gml")
layout <- read.table(paste(out_dir,"layout.table",sep=.Platform$file.sep))
SPADE.plot.trees(mst_graph,out_dir,layout=as.matrix(layout),
                 out_dir=paste(out_dir,"pdf",sep=.Platform$file.sep))
# Record end time
t_SPADE<- Sys.time() - start


# Fit the SPADE results in a FlowSOM object for the visualization 
# used in the paper.
library(FlowSOM)
spade_fcs <- read.FCS(paste0(out_dir,"/21-10-15_Tube_028_selected.fcs.density.fcs.cluster.fcs"))
spade_res <- list("MST"=list(),"map"=list(),"data"=exprs(ff_t[selected,]))
spade_res$MST$graph <- mst_graph
spade_res$MST$l <- as.matrix(layout)
spade_res$map$mapping <- exprs(spade_fcs[,"cluster"])
vsize <- table(spade_res$map$mapping)/nrow(spade_res$map$mapping)
spade_res$MST$size <- (vsize/(max(vsize, na.rm = TRUE)) * 3 + 2) * 4 # same resizing as in SPADE library
spade_res$map$codes <- spade_res$map$medianValues <- t(sapply(seq_len(SPADE_nClusters), function(i) {
    apply(subset(spade_res$data, spade_res$map$mapping[, 1] == i), 
          2, median)
}))
spade_res$map$medianValues[is.nan(spade_res$map$medianValues)] <- 0
spade_res$map$grid <- matrix(1:SPADE_nClusters,nrow=SPADE_nClusters,ncol=1)


# Save results
save(t_SPADE, spade_res, file="SPADE.Rdata")

# Plot results
pdf("SPADE.pdf",useDingbats = FALSE)
PlotPies(UpdateNodeSize(spade_res,reset=T),manual[selected],
         colorPalette = colorRampPalette(c("#FFFFFF",cellTypeColors)))

# Plot individual markers
par(mfrow=c(3,3))
for(m in grid_markerOrder){
    channel <- colnames(ff)[m]
    PlotMarker(spade_res,channel,
               colorPalette = colorRampPalette(c("#dddddd",markerColors[channel])),
               main=markerNames[channel])
}
par(mfrow=c(1,1))
dev.off()

