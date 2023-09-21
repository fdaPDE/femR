# fem-R
The R wrapper to the fdaPDE finite element solver for Partial Differential Equations 

## Installation:

1. Quick installation relying on `devtools` package, from R console: 

      ```> devtools::install_github("fdaPDE/fem-R", ref="develop") ```

2. Cloning the development version of the package.
  
      from terminal:
  
      ``` $ git clone --recurse-submodules -b develop git@github.com:fdaPDE/fem-R.git ```
  
      ``` $ cd path/to/fem-R ```
  
      from R console:
        
      ``` > install.packages(".", type="source", repos=NULL) ```

**Remark** 
femR makes use of git submodules, hence do not download the .zip file from the repository, unzip it and try to install it because the installation procedure will fail. 

