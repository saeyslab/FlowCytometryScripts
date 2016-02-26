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
library(FlowSOM)

# Set some parameters
flowSOM_metaClusters = 10
flowSOM_xdim = 7
flowSOM_ydim = 7



# Set seed for reproducable results
set.seed(42)
# Record start time
start <- Sys.time()
    # Run the flowSOM algorithm
    fsom <- FlowSOM(ff_t[selected,],
                    compensate = FALSE, transform = FALSE, scale = TRUE,
                    colsToUse = colsToCluster, xdim=flowSOM_xdim, ydim=flowSOM_xdim,
                    nClus = flowSOM_metaClusters)
    # Plot the results
    pdf("FlowSOM.pdf",useDingbats = FALSE)
        PlotPies(UpdateNodeSize(fsom[[1]],reset=T),manual[selected],
                 view="MST",
                 colorPalette = colorRampPalette(c('#FFFFFF',cellTypeColors)))
        PlotStars(UpdateNodeSize(fsom[[1]],reset=T),
                  view = "MST",
                  backgroundValues = as.factor(fsom[[2]]),
                  markers = circular_markerOrder,
                  colorPalette = colorRampPalette(markerColors[
                      colnames(ff)[circular_markerOrder]]),
                  backgroundColor = c("#ff7f0055","#FECC0844",
                                      "#bdc9e133","#9b59b644",
                                      '#00000022',"#de2d2655",
                                      "#e74c3c22","#3498db44",
                                      "#addd8e44","#FF000066") 
                  # Colors adapted to manual annotation
        )
        PlotStars(fsom[[1]],
                  view = "MST",
                  backgroundValues = as.factor(fsom[[2]]),
                  markers = circular_markerOrder,
                  colorPalette = colorRampPalette(markerColors[
                      colnames(ff)[circular_markerOrder]]),
                  backgroundColor = c("#ff7f0055","#FECC0844",
                                      "#bdc9e133","#9b59b644",
                                      '#00000022',"#de2d2655",
                                      "#e74c3c22","#3498db44",
                                      "#addd8e44","#FF000066") 
        )
    dev.off()
# Record end time
t_flowSOM<- Sys.time() - start
# Save results
save(t_flowSOM, fsom, file="flowSOM.Rdata")