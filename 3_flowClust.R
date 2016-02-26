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

# Load the flowClust and flowMerge libraries
library(flowClust)
library(flowMerge)

# Set seed for reproducable results
set.seed(42)
# Record start time
start <- Sys.time()
    # Run the flowClust algorithm on the selected cells from the flowFrame
    # Use only the specified columns
    # Look for 10 clusters
    fC <- flowClust(x = ff_t[selected,],
                    varNames = colnames(ff)[colsToCluster],
                    K = 10)
    # Run the flowMerge algorithm on the flowClust result 
    mergeRes<-merge(flowObj(fC,ff[selected,]))
    # Choose the optimal flowMerge result
    fM <- mergeRes[[fitPiecewiseLinreg(mergeRes)]]
    # Group outliers in one extra cluster
    fM@label[is.na(fM@label)] <- max(fM@label,na.rm = TRUE)+1
    # Final result
    res_flowClust <- fM@label
# Record end time
t_flowClust <- Sys.time() - start
# Save results
save(fC,mergeRes, fM,t_flowClust, res_flowClust, file="flowClust.Rdata")


# Repeat analysis with only 10.000 cells

# Set seed for reproducable results
set.seed(42)
# Record start time
start <- Sys.time()
    # Run the flowClust algorithm on the selected cells from the flowFrame
    # Use only the specified columns
    # Look for 10 clusters
    fC <- flowClust(x = ff_t[selected,][1:10000,],
                    varNames = colnames(ff)[colsToCluster],
                    K = 10)
    # Run the flowMerge algorithm on the flowClust result 
    mergeRes<-merge(flowObj(fC,ff[selected,][1:10000,]))
    # Choose the optimal flowMerge result
    fM <- mergeRes[[fitPiecewiseLinreg(mergeRes)]]
    # Group outliers in one extra cluster
    fM@label[is.na(fM@label)] <- max(fM@label,na.rm = TRUE)+1
    # Final result
    res_flowClust_10000 <- fM@label
# Record end time
t_flowClust_10000 <- Sys.time() - start
# Save results
save(fC,mergeRes, fM,t_flowClust_10000,res_flowClust_10000, file="flowClust_10000.Rdata")