[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7064124.svg)](https://doi.org/10.5281/zenodo.7064124)

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/WormTensor)](https://cran.r-project.org/package=WormTensor)
[![Downloads](https://cranlogs.r-pkg.org/badges/WormTensor)](https://CRAN.R-project.org/package=WormTensor)
[![Total Downloads](https://cranlogs.r-pkg.org/badges/grand-total/WormTensor?color=orange)](https://CRAN.R-project.org/package=WormTensor)
--------------------------------------------------
# WormTensor
R package for a Clustering Method for Time-Series Whole-Brain Activity Data of C. elegans

Installation of Dependent Packages
======
~~~~
# GHCR
docker pull ghcr.io/yamaken37/wormtensor
~~~~
or 
~~~~
# Docker Hub
docker pull yamaken37/wormtensor
~~~~

Installation
======
~~~~
git clone https://github.com/rikenbit/WormTensor/
R CMD INSTALL WormTensor
~~~~
or type the code below in the R console window
~~~~
library(devtools)
devtools::install_github("rikenbit/WormTensor")
~~~~
or type the code below in the R console window
~~~~
# CRAN
install.packages(c("WormTensor"), repos="http://cran.r-project.org")
~~~~

## License
Copyright (c) 2022 Kentaro Yamamoto and Laboratory for Bioinformatics Research, RIKEN Center for Biosystems Dynamics Reseach
Released under the [MIT License](https://choosealicense.com/licenses/mit/).

## Authors
- Kentaro Yamamoto
- Koki Tsuyuzaki
- Itoshi Nikaido
