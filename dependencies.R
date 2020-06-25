# R/Cran packages 
install.packages(c("readxl", "enrichR","tidyverse", "stringr", "reshape2"), 
                dependencies = TRUE, repos='http://cran.us.r-project.org')
# R/Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", dependencies = TRUE, repos='http://cran.us.r-project.org')
  BiocManager::install(c("GEOquery", "maanova", "AnnotationDbi", "mgug4122a.db",  "CTDquerier"))
