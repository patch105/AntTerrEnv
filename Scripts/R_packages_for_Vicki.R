
library(dplyr)
library(purrr)
library(here)
library(terra) # Note - I have an HPC terra earlier version not on lib_loc 
library(sf) # Also updated sf so don't need lib.loc=lib_loc
library(ecospat, lib.loc = lib_loc)
library(usdm, lib.loc = lib_loc)
library(ggplot2)
library(tidyterra, lib.loc = lib_loc)
library(dismo, lib.loc = lib_loc)
library(predicts, lib.loc = lib_loc)
library(blockCV, lib.loc = lib_loc)
library(scales)
library(mgcv)
library(randomForest, lib.loc = lib_loc)
library(precrec, lib.loc = lib_loc)
library(glmnet, lib.loc = lib_loc)
library(flexsdm, lib.loc = lib_loc)
library(tidyr)
library(patchwork)
library(viridis)
library(glue)
library(gdalraster, lib.loc = lib_loc)
library(vapour, lib.loc = lib_loc)
library(corrplot)
library(ncdf4, lib.loc = lib_loc)
library(RColorBrewer)
library(ggpubr)
library(forcats)
library(ggtext)
library(arrow)
library(stringr)
library(readr)
library(tibble)
library(lubridate)
library(data.table)
library(DescTools)
library(inlabru, lib.loc=lib_loc)  
library(sf) # NOTE SF MUST BE LOADED BEFORE RISDM
library(RISDM,lib.loc=lib_loc)
library(fmesher,lib.loc=lib_loc)
library(DescTools, lib.loc = lib_loc)
library(precrec, lib.loc = lib_loc)
library(kuenm, lib.loc = lib_loc)
library(prg, lib.loc = lib_loc)


# NOTE * - Also required to have rJava  installed

# # Installing and loading packages
# if(!require(devtools)){
#   install.packages("devtools")
# }
# 
# if(!require(kuenm)){
#   devtools::install_github("marlonecobos/kuenm")
# }

library(kuenm, lib.loc = lib_loc)

# library(devtools)
# install_github('meeliskull/prg/R_package/prg')
library(prg, lib.loc = lib_loc)

# remotes::install_github("rvalavi/myspatial")
library(myspatial, lib.loc = lib_loc)

