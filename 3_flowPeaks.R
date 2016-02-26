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

# Load the flowPeaks library
library(flowPeaks)

# Set seed for reproducable results
set.seed(42)
# Record start time
start <- Sys.time()
    # Run the flowPeaks algorithm on the selected cells from the flowFrame
    # Use only the specified columns
    fP <- flowPeaks(exprs(ff_t[selected,colsToCluster]))
    res_flowPeaks <- fP$peaks.cluster
# Record end time
t_flowPeaks<- Sys.time() - start
# Save results
save(t_flowPeaks, res_flowPeaks, file="flowPeaks.Rdata")


# Repeat analysis with only 10.000 cells
# Set seed for reproducable results
set.seed(42)
# Record start time
start <- Sys.time()
    # Run the flowPeaks algorithm on the selected cells from the flowFrame
    # Use only the specified columns
    fP <- flowPeaks(exprs(ff_t[selected,colsToCluster][1:10000,]))
    res_flowPeaks_10000 <- fP$peaks.cluster
# Record end time
t_flowPeaks_10000 <- Sys.time() - start
# Save results
save(t_flowPeaks_10000, res_flowPeaks_10000, file="flowPeaks_10000.Rdata")
