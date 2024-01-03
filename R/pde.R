
pde_type_list <- list("laplacian" = 1, "elliptic" = 2, "parabolic" = 3)

# Pde Class wraps C++ R_PDE class
.PdeCtr <- R6Class("Pde",
  public = list(
    is_dirichletBC_set = "logical",      
    is_initialCondition_set = "logical",
    is_parabolic ="logical",          
    cpp_handler = "ANY",
    DifferentialOp = "DiffOpObject",
    initialize = function(cpp_handler, DifferentialOp, is_parabolic, is_dirichletBC_set, is_initialCondition_set){
      self$cpp_handler <- cpp_handler
      self$DifferentialOp <- DifferentialOp 
      self$is_parabolic <- is_parabolic
      self$is_dirichletBC_set <- is_dirichletBC_set
      self$is_initialCondition_set <- is_initialCondition_set
    },
    solve = function(){
      if(!self$is_dirichletBC_set){
        if(!self$is_parabolic){
          dirichletBC_ = as.matrix(rep(0,times=nrow(self$cpp_handler$get_dofs_coordinates())))
        }else 
          dirichletBC_ = matrix(0, nrow=nrow(self$cpp_handler$get_dofs_coordinates()), 
                                   ncol=length(self$DifferentialOp$f$FunctionSpace$mesh$times))
      self$cpp_handler$set_dirichlet_bc(dirichletBC_)
      self$is_dirichletBC_set <- TRUE
      }
      if(self$is_parabolic & (!self$is_initialCondition_set)){
        stop("initialCondition must be provided.")
      }
      self$cpp_handler$solve()
    },
    get_dofs_coordinates = function(){
      self$cpp_handler$get_dofs_coordinates()
    },
    set_boundary_condition = function(fun, type="dirichlet", on=NULL){
      if(!any(type == c("dirichlet", "Dirichlet", "d"))) 
        stop("Only Dirichlet boundary condtions allowed.")
      dirichletBC_ <- NULL
      if(any(typeof(fun) == c("function", "closure"))){
        if(!self$is_parabolic){
          dirichletBC_ <- as.matrix(fun(self$cpp_handler$get_dofs_coordinates()))
        }else{
          dirichletBC_ <- fun(self$cpp_handler$get_dofs_coordinates(), 
                              self$DifferentialOp$f$FunctionSpace$mesh$times)
        }
      }else if(any(typeof(fun) == c("matrix","vector", "numeric" ,"double"))){
        if(nrow(as.matrix(fun)) == nrow(self$cpp_handler$get_dofs_coordinates())){
          dirichletBC_ <- fun  
        }else if (nrow(as.matrix(fun)) == 1L){
          if(!self$is_parabolic)
            dirichletBC_ <- matrix(fun, nrow=nrow(self$cpp_handler$get_dofs_coordinates()), ncol=1)
          else 
            dirichletBC_ <- matrix(fun, nrow=nrow(self$cpp_handler$get_dofs_coordinates()), 
                                        ncol=length(self$DifferentialOp$f$FunctionSpace$mesh$times))
          }   
        }
      self$is_dirichletBC_set <- TRUE
      self$cpp_handler$set_dirichlet_bc(dirichletBC_)
    },
    set_initial_condition = function(fun){
      if(!self$is_parabolic)
        stop("Cannot set initial condition for elliptic problem.")
      self$is_initialCondition_set <- TRUE
      if(typeof(fun) == "closure"){ 
        self$cpp_handler$set_initial_condition(fun(self$cpp_handler$get_dofs_coordinates()))
      }else if(any(typeof(fun) == c("matrix","vector", "numeric" ,"double"))){
        if(nrow(as.matrix(fun)) == nrow(self$cpp_handler$get_dofs_coordinates())){
          self$cpp_handler$set_initial_condition(fun)
        }else if(nrow(as.matrix(fun)) == 1L){
          self$cpp_handler$set_initial_condition(matrix(fun, nrow=nrow(self$cpp_handler$get_dofs_coordinates()),
                                                 ncol=1))   
        }
      }
    },
    get_mass = function(){
      self$cpp_handler$mass()
    },
    get_stiff = function(){
      self$cpp_handler$stiff()
    }
  )
)

#' @name pde
#'
#' @exportClass Pde
setOldClass(c("Pde", "R6"))

## infers the type of a pde
extract_pde_type <- function(DifferentialOp) {
  if ("time" %in% names(DifferentialOp$params)) {
    return(pde_type_list$parabolic)
  }
  if ("diffusion" %in% names(DifferentialOp$params) && !is.matrix(DifferentialOp$params$diffusion)) {
    return(pde_type_list$laplacian)
  }
  return(pde_type_list$elliptic)
}

## parse pde parameters
parse_pde_parameters <- function(DifferentialOp){
  pde_type <- extract_pde_type(DifferentialOp)
  
  pde_parameters <- NULL
  pde_parameters$diffusion <- 1.0
  pde_parameters$transport <- matrix(0, nrow = 2, ncol = 1)
  pde_parameters$reaction  <- 0.0
  
  for (i in 1:length(DifferentialOp$params)) {
    pde_parameters[[DifferentialOp$tokens[i]]] <- DifferentialOp$params[[DifferentialOp$tokens[i]]]
  }
  
  if(pde_type == pde_type_list$parabolic)
    pde_parameters[["time_mesh"]] <- as.vector(DifferentialOp$f$FunctionSpace$mesh$times)
  
  if(pde_type == pde_type_list$parabolic & is(pde_parameters[["diffusion"]], "numeric"))
    pde_parameters[["diffusion"]] <- pde_parameters[["diffusion"]]*matrix(c(1,0,0,1), nrow=2, ncol=2, byrow=T)
  
  return(pde_parameters)
}

## build cpp object
make_pde <- function(DifferentialOp) {
  pde_type <- extract_pde_type(DifferentialOp)
  
  ## pde parameters
  pde_parameters <- parse_pde_parameters(DifferentialOp)
  
  ## define Rcpp module
  D <- DifferentialOp$f$FunctionSpace$mesh$cpp_handler ## domain
  fe_order <- DifferentialOp$f$FunctionSpace$fe_order  ## finite element order
  cpp_handler <- NULL
  if (fe_order == 1) { ## linear finite elements
    cpp_handler <- new(cpp_pde_2d_fe1, D, pde_type - 1, pde_parameters, DifferentialOp$f)
  }
  if (fe_order == 2) { ## quadratic finite elements
    cpp_handler <- new(cpp_pde_2d_fe2, D, pde_type - 1, pde_parameters, DifferentialOp$f)
  }
  
  return(cpp_handler)
}

#' A PDEs object
#'
#' @param DifferentialOp a differential operator.
#' @param f a standard R function representing the forcing term of the PDE or a numeric value, in case of constant forcing term.
#' @return A R6 class representing a PDE.
#' @rdname pde
#' @export 
setGeneric("Pde", function(DifferentialOp,f) standardGeneric("Pde"))

#' @rdname pde
setMethod("Pde", signature=c(DifferentialOp="DiffOpObject", f="function"),
          function(DifferentialOp,f){
            
            pde_type = extract_pde_type(DifferentialOp)
            cpp_handler <- make_pde(DifferentialOp)
            
            quad_nodes <- cpp_handler$get_quadrature_nodes()
            ## evaluate forcing term on quadrature nodes
            if(pde_type == pde_type_list$parabolic) {
              times <- DifferentialOp$f$FunctionSpace$mesh$times
              cpp_handler$set_forcing(as.matrix(f(quad_nodes, times)))
            }else{
              cpp_handler$set_forcing(as.matrix(f(quad_nodes)))
            }
            
            ## initialize solver 
            cpp_handler$init()
            
            is_parabolic = pde_type == pde_type_list$parabolic
            
            ## return
            .PdeCtr$new(is_dirichletBC_set = FALSE,
                        is_initialCondition_set = FALSE,
                        is_parabolic = is_parabolic,
                        cpp_handler = cpp_handler,
                        DifferentialOp = DifferentialOp)
})

#' @rdname pde
setMethod("Pde", signature=c(DifferentialOp="DiffOpObject", f="numeric"),
          function(DifferentialOp,f){
            pde_type = extract_pde_type(DifferentialOp)
            
            fun <- NULL
            if(pde_type!=pde_type_list$parabolic){
              fun <- function(points){
                return( matrix(f, nrow=nrow(points), ncol=1))
              }
            }else{
              fun <- function(points, times){
                return( matrix(f, nrow=nrow(points), ncol=length(times)))
              }  
            }
            
            return(Pde(DifferentialOp,fun))
})
