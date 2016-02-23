
cellTypeColors <- c("B cells" = "#e74c3c",
                    "T cells" = "#a4cc2e",
                    "NK cells" = "#3498db",
                    "NK T cells" = "#a65628",
                    "Macrophages" = "#ff7f00",
                    "DCs" = "#FECC08",
                    "Neutrophils" = "#9b59b6")

markerNames <- c("AmCyan-A" = "Autofluorescence",
                 "Comp-BV711-A" = "CD64",
                 "Comp-PE-Cy5-A" = "CD19",
                 "Comp-PerCP-Cy5-5-A" = "MHCII",
                 "Comp-PE-Cy7-A" = "CD11c",
                 "Comp-BV605-A" = "CD11b",
                 "Comp-Alexa Fluor 700-A" = "Ly-6G",
                 "Comp-APC-A"="NK1.1",
                 "Comp-Pacific Blue-A"="CD49b",
                 "Comp-PE-A"="CD3",
                 "Comp-BV786-A"="FcERI")
markerColors <- c("AmCyan-A" = "#cc4c02",
                  "Comp-BV711-A" = "#fdae6b",
                  "Comp-PE-Cy5-A" = "#e41a1c",
                  "Comp-PerCP-Cy5-5-A" = "#ff7f00",
                  "Comp-PE-Cy7-A" = "#FECC08",
                  "Comp-BV605-A" = "#ae017e",
                  "Comp-Alexa Fluor 700-A" = "#9b59b6",
                  "Comp-APC-A"="#3498db",
                  "Comp-Pacific Blue-A"="#1f78b4",
                  "Comp-PE-A"="#a4cc2e",
                  "Comp-BV786-A"="#fb9a99")

circular_markerOrder = c(10,12,18,8,19,11,15,14,17)
grid_markerOrder = c(8,18,12,19,17,10,11,15,14)

save(cellTypeColors,markerColors,markerNames,circular_markerOrder,grid_markerOrder,
     file="plotSettings.Rdata")
