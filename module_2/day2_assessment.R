options(java.parameters = "-Xmx8g")
library(VoltRon)

####
# Visium Analysis ####
####

####
## Import Data ####
####

# Dependencies
if(!requireNamespace("rhdf5"))
  BiocManager::install("rhdf5")
library(rhdf5)

# import Visium data
BC_sample <- importVisium("./module_2/workshop/data/BreastCancer/Visium",
                         sample_name = "Breast_Cancer")
# features
head(vrFeatures(BC_sample))
length(vrFeatures(BC_sample))

# normalize and select features
BC_sample <- normalizeData(BC_sample)
BC_sample <- getFeatures(BC_sample, n = 3000)

# selected features
head(vrFeatureData(BC_sample))
selected_features <- getVariableFeatures(BC_sample)
head(selected_features, 20)

####
### Dimensional Reduction ####
####

# embedding
BC_sample <- getPCA(BC_sample, features = selected_features, dims = 30)
BC_sample <- getUMAP(BC_sample, dims = 1:30)
vrEmbeddingNames(BC_sample)

# embedding visualization
vrEmbeddingPlot(BC_sample, embedding = "umap")

####
### Clustering ####
####

# graph for neighbors
BC_sample <- getProfileNeighbors(BC_sample, dims = 1:30, k = 10, method = "SNN")
vrGraphNames(BC_sample)

# clustering
BC_sample <- getClusters(BC_sample, resolution = 0.5, label = "Clusters", graph = "SNN")

####
### Visualization ####
####

# embedding
pdf("./module_2/KH_quiz_1.pdf")
vrEmbeddingPlot(BC_sample, embedding = "umap", group.by = "Clusters")
vrSpatialPlot(BC_sample, group.by = "Clusters")
dev.off()

