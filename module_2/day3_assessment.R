options(java.parameters = "-Xmx32g")
library(VoltRon)
library(dplyr)


# Import analyzed data
BC_vis_sample <- importVisium("./module_2/workshop/data/BreastCancer/Visium",
                          sample_name = "Breast_Cancer")


# load H&E data from disk
Xen_img <- loadVoltRon("./module_2/workshop/data/ondisk/Xen_R1/")


# automated alignment
xen_reg <- registerSpatialData(object_list = list(BC_vis_sample, Xen_img))
