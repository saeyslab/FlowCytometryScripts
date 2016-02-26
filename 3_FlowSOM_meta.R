# Load the preprocessed data:
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

# Load the FlowSOM library
library(FlowSOM)

# Set seed for reproducable results
set.seed(42)
# Record start time
start <- Sys.time()
    # Run the FlowSOM algorithm on the selected cells from the flowFrame
    # Use only the specified columns
    fsom <- FlowSOM(ff_t[selected,],colsToUse = colsToCluster,nClus=10)
    res_FlowSOM <- fsom[[2]][fsom[[1]]$map$mapping[,1]]
# Record end time
t_FlowSOM<- Sys.time() - start
# Save results
save(t_FlowSOM, res_FlowSOM, file="FlowSOM_meta.Rdata")


# Repeat analysis with only 10.000 cells

# Set seed for reproducable results
set.seed(42)
# Record start time
start <- Sys.time()
    # Run the FlowSOM algorithm on the selected cells from the flowFrame
    # Use only the specified columns
    fsom <- FlowSOM(ff_t[selected,][1:10000,],colsToUse = colsToCluster,nClus=10)
    res_FlowSOM_10000 <- fsom[[2]][fsom[[1]]$map$mapping[,1]]
# Record end time
t_FlowSOM_10000<- Sys.time() - start
# Save results
save(t_FlowSOM_10000, res_FlowSOM_10000, file="FlowSOM_meta_10000.Rdata")
