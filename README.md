[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7064124.svg)](https://doi.org/10.5281/zenodo.7064124)

[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/WormTensor)](https://cran.r-project.org/package=WormTensor)
[![Downloads](https://cranlogs.r-pkg.org/badges/WormTensor)](https://CRAN.R-project.org/package=WormTensor)
[![Total Downloads](https://cranlogs.r-pkg.org/badges/grand-total/WormTensor?color=orange)](https://CRAN.R-project.org/package=WormTensor)

[![WormTensor status badge](https://rikenbit.r-universe.dev/badges/WormTensor)](https://rikenbit.r-universe.dev)
[![:name status badge](https://rikenbit.r-universe.dev/badges/:name)](https://rikenbit.r-universe.dev)
[![:registry status badge](https://rikenbit.r-universe.dev/badges/:registry)](https://rikenbit.r-universe.dev)
[![:total status badge](https://rikenbit.r-universe.dev/badges/:total)](https://rikenbit.r-universe.dev)
![GitHub Actions](https://github.com/rikenbit/WormTensor/actions/workflows/build_push_test.yml/badge.svg)
--------------------------------------------------
# WormTensor
R package for a Clustering Method for Time-Series Whole-Brain Activity Data of C. elegans

<img src="WormTensor.png" width="120px" height="139px" />

Installation of Dependent Packages
======
~~~~
# GHCR
docker pull ghcr.io/rikenbit/wormtensor
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
