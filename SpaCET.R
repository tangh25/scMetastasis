# Install required packages if needed
# BiocManager::install("genefilter")
# devtools::install_github("JEFworks/MUDAN")
# devtools::install_github("data2intelligence/SpaCET")

# Specify Python environment for reticulate (required by SpaCET)
Sys.setenv(RETICULATE_PYTHON="/home/htang/.conda/envs/TumorBoundary/bin/python.exe")

library(SpaCET)

############################################################
# Load Visium spatial transcriptomics data
############################################################

# Example Visium dataset from SpaCET package
visiumPath = file.path(system.file(package = "SpaCET"), "extdata/Visium_BC")

# Create SpaCET object from 10X Visium data
SpaCET_obj = create.SpaCET.object.10X(visiumPath = visiumPath)

SpaCET_obj@input$counts[1:8,1:6]

############################################################
# Quality control
############################################################

SpaCET_obj = SpaCET.quality.control(SpaCET_obj)

# Visualize spatial distribution of UMI counts and gene numbers
SpaCET.visualize.spatialFeature(
  SpaCET_obj,
  spatialType = "QualityControl",
  spatialFeatures = c("UMI", "Gene"),
  imageBg = TRUE
)

############################################################
# Cell type deconvolution
############################################################

SpaCET_obj = SpaCET.deconvolution(
  SpaCET_obj,
  cancerType = "BRCA",   # Breast cancer reference
  coreNo = 8             # Number of CPU cores
)

SpaCET_obj@results$deconvolution$propMat[1:13, 1:6]

############################################################
# Visualization of cell fractions
############################################################
# Visualize all primary and secondary cell types
SpaCET.visualize.spatialFeature(
  SpaCET_obj,
  spatialType = "CellFraction",
  spatialFeatures = "All",
  sameScaleForFraction = TRUE,
  pointSize = 0.1,
  nrow = 5
)

# Interactive visualization
SpaCET.visualize.spatialFeature(
  SpaCET_obj,
  interactive = TRUE
)

############################################################
# Cell–cell interaction (CCI) analysis
############################################################

# Colocalization analysis
SpaCET_obj = SpaCET.CCI.colocalization(SpaCET_obj)

SpaCET.visualize.colocalization(SpaCET_obj)

# Ligand–receptor network score
SpaCET_obj = SpaCET.CCI.LRNetworkScore(
  SpaCET_obj,
  coreNo = 8
)

SpaCET.visualize.spatialFeature(
  SpaCET_obj,
  spatialType = "LRNetworkScore",
  spatialFeatures = c("Network_Score", "Network_Score_pv")
)

############################################################
# Tumor boundary analysis
############################################################

# Identify tumor–stroma interface
SpaCET_obj = SpaCET.identify.interface(SpaCET_obj)

# Visualize tumor–stroma interface
SpaCET.visualize.spatialFeature(
  SpaCET_obj,
  spatialType = "Interface",
  spatialFeatures = "Interface"
)
