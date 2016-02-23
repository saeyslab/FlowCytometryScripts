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
load("fcs/21-10-15_Tube_028.Rdata")

# Load the immunoClust library
library(immunoClust)

# Set seed for reproducable results
set.seed(42)
# Record start time
start <- Sys.time()
    # Run the immunoClust algorithm on the selected cells from the flowFrame
    # Use only the specified columns
    iC <- cell.process(ff[selected,],parameters = colnames(ff)[colsToCluster],classify.all = TRUE)
    res_immunoClust <- iC@label
# Record end time
t_immunoClust <- Sys.time() - start
# Save results
save(t_immunoClust, res_immunoClust, file="immunoClust.Rdata")


# Repeat analysis with only 10.000 cells
# Set seed for reproducable results
set.seed(42)
# Record start time
start <- Sys.time()
    # Run the immunoClust algorithm on the selected cells from the flowFrame
    # Use only the specified columns
    iC <- cell.process(ff[selected,][1:10000,],
                       parameters = colnames(ff)[colsToCluster],
                       classify.all = TRUE,
                       N=10000)
    res_immunoClust_10000 <- iC@label
# Record end time
t_immunoClust_10000 <- Sys.time() - start
# Save results
save(t_immunoClust_10000, res_immunoClust_10000, file="immunoClust_10000.Rdata")
