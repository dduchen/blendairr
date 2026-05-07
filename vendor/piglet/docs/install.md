# Download & Installation

PIgLET is available for installation from the development version.

To build from the source code, first install the build dependencies:
  
  ```{r, eval=FALSE}
install.packages(c("devtools", "roxygen2", "testthat", "knitr", "rmarkdown", "plotly"))
```

To install the latest version via devtools:
  
  ```{r, eval=FALSE}
library(devtools)
install_bitbucket("yaarilab/piglet")
```

Note, installing from bitbucket does not generate the documentation. To
generate them, first clone the repository and then build:
  
  ```{r, eval=FALSE}
library(devtools)
install_deps()
document()
build()
install()
```
