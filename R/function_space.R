.BasisFunctionCtr <- R6Class("BasisFunction",
                              private = list(
                                cpp_handler = "ANY",  ## pointer to cpp backend
                                mesh = "Mesh"
                              ),
                              public = list(
                                initialize = function(cpp_handler, mesh){
                                  private$cpp_handler = cpp_handler
                                  private$mesh = mesh
                                },
                                size = function(){
                                  private$cpp_handler$size()
                                },
                                get_dofs_coordinates = function(){
                                  private$cpp_handler$get_dofs_coordinates()
                                },
                                ## evaluates the basis system on the given set of locations
                                eval = function(locations) {
                                  return(private$cpp_handler$eval(0L, locations))
                                },
                                ## integrates a function expressed as basis expansion with respect to this basis system over the
                                ## whole domain of definition
                                integrate = function(func) {
                                  if (!(is.function(func) || is.vector(func))) stop("invalid argument, should be either a function or a vector")
                                  ## required dof coordinates...
                                  c <- func
                                  if (is.function(func)) { ## recover fem expansion
                                    c <- func(private$cpp_handler$get_dofs_coordinates())
                                  }
                                  return(private$cpp_handler$integrate(c))
                                }
  )
)

#' Finite Element Basis 
#'
#' @name Basis
#'
#' @exportClass BasisFunction
setOldClass(c("BasisFunction", "R6"))

setGeneric("BasisFunction", function(mesh, fe_order) standardGeneric("BasisFunction"))

setMethod("BasisFunction", signature = signature("Mesh","integer"),
          function(mesh, fe_order){
            basis_function <- .BasisFunctionCtr$new(cpp_handler = 
                                                  new(eval(parse(text = paste("cpp_lagrange_basis_2d_fe", 
                                                                              as.character(fe_order), sep = ""))), extract_private(mesh)$cpp_handler, 0),
                                                  mesh=mesh)
})

#' @export 
.FunctionSpaceCtr <- R6Class("FunctionSpace",
                              private= list(
                                mesh = "Mesh",
                                basis_function = "BasisFunction",     ## cpp backend
                                fe_order = "integer"                  ## this is specific for fem, must be generalized
                              ),
                              public = list(
                                initialize = function(mesh, basis_function, fe_order){
                                  private$mesh <- mesh 
                                  private$basis_function <- basis_function
                                  private$fe_order <- fe_order
                                },
                                get_basis = function(){ 
                                  return(private$basis_function) 
                                }
                              )
)

#' Finite Element Functional Space
#'
#' @name FunctionSpace
#'
#' @exportClass FunctionSpace
setOldClass(c("FunctionSpace", "R6"))

#' Create FunctionSpace object
#'
#' @param mesh A mesh object created by \code{Mesh}:
#' @param fe_order Either '1' or '2'. It specifies the finite element order.
#' @return An S4 object representing a Function Space.
#' @export 
#' @rdname FunctionSpace
setGeneric("FunctionSpace", function(mesh, fe_order) standardGeneric("FunctionSpace"))

#' @rdname FunctionSpace
#' @examples
#' \dontrun{
#' library(femR)
#' data("unit_square")
#' mesh <- Mesh(unit_square)
#' Vh <- FunctionSpace(mesh = mesh, fe_order = 1)
#' }
setMethod("FunctionSpace",
          signature = c(mesh = "Mesh", fe_order = "numeric"),
          function(mesh, fe_order) {
            .FunctionSpaceCtr$new(
              mesh = mesh,
              basis_function = BasisFunction(mesh, as.integer(fe_order)),
              fe_order = as.integer(fe_order)
            )
          }
)

#' @rdname FunctionSpace
#' @examples
#' \dontrun{
#' library(femR)
#' data("unit_square")
#' mesh <- Mesh(unit_square)
#' Vh <- FunctionSpace(mesh)
#' }
setMethod("FunctionSpace", signature = c(mesh="Mesh", fe_order="missing"),
          function(mesh){
              return(.FunctionSpaceCtr$new(mesh=mesh, 
                                           basis_function = BasisFunction(mesh, 1L), 
                                           fe_order=1L))
})

## finite element function
.FunctionCtr <- R6Class("Function",
                        private = list(
                          FunctionSpace  = "FunctionSpace" 
                        ),
                        public = list(
                          coeff = "matrix",
                          initialize = function(FunctionSpace, coeff){
                            private$FunctionSpace <- FunctionSpace
                            self$coeff <- coeff
                          },
                          eval = function(X) {
                            M = dim(extract_private(private$FunctionSpace)$mesh$get_nodes())[2]
                            if(is.vector(X)) X <- t(as.matrix(X)) 
                            
                            if(dim(X)[2] != M) {
                              stop(paste("matrix of evaluation points should be a 2 columns matrix"))
                            }
                            evals <- NULL  
                            if(length(extract_private(private$FunctionSpace)$mesh$get_times()) == 0){
                              evals <- apply(private$FunctionSpace$get_basis()$eval(as.matrix(X)) %*% self$coeff, 
                                             MARGIN=1, FUN=sum)
                            }else{
                              evals <- matrix(nrow=nrow(X),ncol=length(extract_private(private$FunctionSpace)$mesh$get_times()))
                              for(t in 1:length(extract_private(private$FunctionSpace)$mesh$get_times())){
                                evals[,t] <- apply(private$FunctionSpace$get_basis()$eval(as.matrix(X)) %*% self$coeff[,t], 
                                                   MARGIN=1, FUN=sum) 
                              }  
                            }
                            return(evals)
                          },
                          set_coefficients = function(coefficients){
                            if( nrow(coefficients) != private$FunctionSpace$get_basis$size())
                              stop("Input parameter dimensions are different from the size of the function space.")
                            self$coeff <- coefficients
                          }
                        )
)

#' Finite Element Function
#'
#' @name Function
#'
#' @exportClass Function
setOldClass(c("Function", "R6"))

## constructor

#' Create Function object
#'
#' @param FunctionSpace a FunctionSpace object created by \code{FunctionSpace}:
#' @return An S4 object representing a Function belonging to the FunctionSpace passed as parameter.
#' @export 
#' @examples
#' \dontrun{
#' library(femR)
#' data("unit_square")
#' mesh <- Mesh(unit_square)
#' Vh <- FunctionSpace(mesh)
#' f <- Function(Vh)
#' }
Function <- function(FunctionSpace) {
  coeff = matrix(ncol = 1, nrow = 0)
  .FunctionCtr$new(coeff = coeff, FunctionSpace = FunctionSpace)
}

## Function plot overload
# setMethod("plot", signature=c(x="Function"),

#' Plot a Function object
#'
#' @param x A \code{Function} generated by \code{Function}
#' @param ... Arguments representing graphical options to be passed to \code{\link[plotly]{plot_ly}}.
#' @return A plotly object
#' @importFrom graphics plot
#' @export
plot.Function <- function(x, ...){
  mesh <- extract_private(extract_private(x)$FunctionSpace)$mesh
  times <- mesh$get_times()
  is_parabolic <- FALSE
  if(length(times)!=0) is_parabolic <- TRUE
  if(!is_parabolic){
  plot_data <- data.frame(X=mesh$get_nodes()[,1], 
                          Y=mesh$get_nodes()[,2],
                          coeff=x$coeff[1:nrow(mesh$get_nodes())])
  I=mesh$get_elements()[,1]-1
  J=mesh$get_elements()[,2]-1
  K=mesh$get_elements()[,3]-1
  fig<- plot_ly(plot_data, x=~X, y=~Y, z=~coeff,
          i = I, j = J, k = K,
          intensity=~coeff, color=~coeff, type="mesh3d", 
          colorbar=list(title=""), ...) %>%
      layout(scene = list(
        aspectmode = "data", 
        xaxis = list(
          title = '', showgrid = F, zeroline = F, showticklabels = F),
        yaxis = list(
          title = '', showgrid = F, zeroline = F, showticklabels = F),
        zaxis = list(
          title = '', showgrid = F, zeroline = F, showticklabels = F)),
        camera = list(
          eye = list(x = 1.25, y = -1.25, z = 1.25))) %>%
      colorbar(len = 1, title="")
  }else{
    plot_data <- data.frame(X=rep(mesh$get_nodes()[,1], times=length(times)), 
                            Y=rep(mesh$get_nodes()[,2], times=length(times)),
                            Z=rep(0, times=length(times)),
                            coeff=as.vector(x$coeff[1:nrow(mesh$get_nodes()),]),
                            times = round(rep(times, each=nrow(mesh$get_nodes())),3))
    coeff <- matrix(nrow=nrow(mesh$get_elements()), 
                    ncol=length(times))
    for(t in 1:length(times)){
      coeff[,t] <- apply(mesh$get_elements(), MARGIN=1, FUN=
                           function(edge){
                             mean(x$coeff[edge,t])
                           })
    }
    limits = c(min(coeff), max(coeff))
    cmin = limits[1]; cmax=limits[2]
    I=mesh$get_elements()[,1]-1
    J=mesh$get_elements()[,2]-1
    K=mesh$get_elements()[,3]-1
    fig<- plot_ly(plot_data, x=~X, y=~Y, z=~Z, frame=~times,
                  i = I, j = J, k = K, cmin = limits[1], cmax=limits[2],
                  intensity=~coeff, color=~coeff, type="mesh3d",
                  colorbar=list(title=""), ...) %>%
        layout(scene = list(
          aspectmode = "data", 
          xaxis = list(
            title = '', showgrid = F, zeroline = F, showticklabels = F),
          yaxis = list(
            title = '', showgrid = F, zeroline = F, showticklabels = F),
          zaxis = list(
            title = '', showgrid = F, zeroline = F, showticklabels = F),
          camera = list(
            eye = list(x = 0, y = -0.01,  z = 1.25))), dragmode="zoom") %>%
        colorbar(len = 1, title="") %>%
        animation_slider(currentvalue = list(prefix ="t = "))
        #animation_opts(frame=5) %>% 
  }
  fig 
}

## Function contour overload

#' Create a contour plot of a Function
#'
#' @param x A \code{Function} generated by \code{Function}
#' @param ... Arguments representing graphical options to be passed to \code{\link[plotly]{plot_ly}}.
#' @return A plotly object
#' @export
setMethod("contour", signature=c(x="Function"), function(x, ...){
  mesh <- extract_private(extract_private(x)$FunctionSpace)$mesh
  times <- mesh$get_times()
  is_parabolic <- FALSE
  if(length(times)!=0) is_parabolic <- TRUE
  
  if(!is_parabolic){
  xrange <- range(mesh$get_nodes()[,1])
  yrange <- range(mesh$get_nodes()[,2])
  Nx <- 40
  Ny <- 40
  eval_x <- seq(from=xrange[1], to=xrange[2], length.out=Nx)
  eval_y <- seq(from=yrange[1], to=yrange[2], length.out=Ny)
  eval_points <- expand.grid(eval_x, eval_y)
  Z <- matrix(x$eval(eval_points), nrow=Nx,ncol=Ny)
  fig <- plot_ly(type="contour", x=eval_x, y=eval_y, z=Z, 
                 intensity=Z, color = Z,
                 contours=list(showlabels = TRUE),
                 colorbar=list(title=""), ...) %>%
    layout(xaxis = list(title = "", showgrid=F, zeroline=F, ticks="", showticklabels=F),
           yaxis = list(title = "", showgrid=F, zeroline=F, ticks="", showticklabels=F))
  }else{
    
    xrange <- range(mesh$get_nodes()[,1])
    yrange <- range(mesh$get_nodes()[,2])
    Nx <- 40
    Ny <- 40
    eval_x <- seq(from=xrange[1], to=xrange[2], length.out=Nx)
    eval_y <- seq(from=yrange[1], to=yrange[2], length.out=Ny)
    eval_points <- expand.grid(eval_x, eval_y)
    
    eval_x <- seq(from=xrange[1], to=xrange[2], length.out=Nx)
    eval_y <- seq(from=yrange[1], to=yrange[2], length.out=Ny)
    eval_points <- expand.grid(eval_x, eval_y)
    X <- rep(eval_points[,1], times=length(times))
    Y <- rep(eval_points[,2], times=length(times))
    Z <- x$eval(eval_points)

    plot_data <- data.frame(X=X, Y=Y,Z=as.vector(Z),
                            times = rep(times, each=nrow(eval_points)))
    limits = c(min(plot_data$Z), max(plot_data$Z))
    fig <- plot_ly(plot_data, type="contour", x=~X, y=~Y, z=~Z, frame=~times, 
                   zmin=limits[1], zmax=limits[2],
                   intensity=~Z, color = ~Z,
                   contours=list(showlabels = TRUE),
                   colorbar=list(title=""), ...) %>%
      layout(xaxis = list(title = "", showgrid=F, zeroline=F, ticks="", showticklabels=F),
             yaxis = list(title = "", showgrid=F, zeroline=F, ticks="", showticklabels=F)) %>%
      animation_slider(currentvalue = list(prefix ="t = ")) #animation_opts(frame=5) %>% 
    
  }
  fig
}
)