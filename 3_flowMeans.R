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

# Load the flowMeans library
library(flowMeans)

# Set seed for reproducable results
set.seed(42)
# Record start time
start <- Sys.time()
    # Run the flowMeans algorithm on the selected cells from the flowFrame
    # Use only the specified columns
    # Look for 10 clusters
    fM <- flowMeans(x = ff_t[selected,],
                    varNames = colnames(ff)[colsToCluster],
                    NumC = 10)
    res_flowMeans <- fM@Label
# Record end time
t_flowMeans <- Sys.time() - start
# Save results
save(fM, t_flowMeans, res_flowMeans, file="flowMeans.Rdata")


# Repeat analysis with only 10.000 cells
# Set seed for reproducable results
set.seed(42)
# Record start time
start <- Sys.time()
    # Run the flowMeans algorithm on the selected cells from the flowFrame
    # Use only the specified columns
    # Look for 10 clusters
    fM <- flowMeans(x = ff_t[selected,][1:10000,],
                    varNames = colnames(ff)[colsToCluster],
                    NumC = 10)
    res_flowMeans_10000 <- fM@Label
# Record end time
t_flowMeans_10000 <- Sys.time() - start
# Save results
save(fM, t_flowMeans_10000, res_flowMeans_10000, file="flowMeans_10000.Rdata")
