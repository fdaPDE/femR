# femR
The R wrapper to the fdaPDE finite element solver for Partial Differential Equations 

> This is a first pre-release of the package, as such, bugs might occur. Feel free to open an issue in case of problems.

## Installation

Make sure to have the following dependencies installed on your system:

* a C++17 compliant compiler
* the `Rcpp` and `RcppEigen` packages 

then, to install the latest stable version of `femR`, you can either:

1. use the `devtools` package. From the R console, execute

      ```
	  devtools::install_github("fdaPDE/femR", ref="stable") 
	  ```

2. clone this repository and install. From a terminal, execute

      ``` 
	  git clone --recurse-submodules -b stable git@github.com:fdaPDE/femR.git 
      cd path/to/femR 
	  ```

	and install the package from the R console

	``` 
	  install.packages(".", type="source", repos=NULL) 
	  ```

Both procedures will automatically pull the [fdaPDE-core](https://github.com/fdaPDE/fdaPDE-core) submodule dependence required by `femR`. It is not recommended to download the source code directly from Github, as this won't include any submodule dependence, making the installation procedure to fail.
