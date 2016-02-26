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

# Load the SamSPECTRAL library
library(SamSPECTRAL)

# Set seed for reproducable results
set.seed(42)
# Record start time
start <- Sys.time()
    # Run the SamSPECTRAL algorithm on the selected cells from the flowFrame
    # Use only the specified columns
    # Look for 10 clusters
    ss <- SamSPECTRAL(exprs(ff_t[selected,]),
                      dimensions = colsToCluster,
                      normal.sigma = 200, separation.factor = 0.39,
                      number.of.clusters = 10)
    # Add outliers as another cluster
    ss[is.na(ss)] <- max(ss,na.rm=TRUE)+1
    res_SamSPECTRAL <- ss
# Record end time
t_SamSPECTRAL <- Sys.time() - start
# Save results
save(t_SamSPECTRAL, res_SamSPECTRAL, file="SamSPECTRAL.Rdata")


# Repeat analysis with only 10.000 cells
# Set seed for reproducable results
set.seed(42)
# Record start time
start <- Sys.time()
    # Run the SamSPECTRAL algorithm on the selected cells from the flowFrame
    # Use only the specified columns
    # Look for 10 clusters
    ss <- SamSPECTRAL(exprs(ff_t[selected,][1:10000,]),
                      dimensions = colsToCluster,
                      normal.sigma = 200, separation.factor = 0.39,
                      number.of.clusters = 10)
    # Add outliers as another cluster
    ss[is.na(ss)] <- max(ss,na.rm=TRUE)+1
    res_SamSPECTRAL_10000 <- ss
# Record end time
t_SamSPECTRAL_10000 <- Sys.time() - start
# Save results
save(t_SamSPECTRAL_10000, res_SamSPECTRAL_10000, file="SamSPECTRAL_10000.Rdata")
