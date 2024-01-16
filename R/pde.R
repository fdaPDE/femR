
pde_type_list <- list("laplacian" = 1, "elliptic" = 2, "parabolic" = 3)

# Pde Class wraps C++ R_PDE class
.PdeCtr <- R6Class("Pde",
  private = list(
    is_dirichletBC_set = "logical",      
    is_initialCondition_set = "logical",
    is_parabolic ="logical",          
    cpp_handler_ = "ANY",
    DifferentialOp_ = "DifferentialOp"
  ),
  public = list(
    initialize = function(cpp_handler, DifferentialOp, is_parabolic, is_dirichletBC_set, is_initialCondition_set){
      private$cpp_handler_ <- cpp_handler
      private$DifferentialOp_ <- DifferentialOp 
      private$is_parabolic <- is_parabolic
      private$is_dirichletBC_set <- is_dirichletBC_set
      private$is_initialCondition_set <- is_initialCondition_set
    },
    solve = function(){
      if(!private$is_dirichletBC_set){
        if(!private$is_parabolic){
          dirichletBC_ = as.matrix(rep(0,times=nrow(private$cpp_handler_$dofs_coordinates())))
        }else{ 
          times <- private$DifferentialOp_$Function()$FunctionSpace()$mesh()$time_nodes()
          dirichletBC_ = matrix(0, nrow=nrow(private$cpp_handler_$dofs_coordinates()), 
                                   ncol=length(times))
        }
      private$cpp_handler_$set_dirichlet_bc(dirichletBC_)
      private$is_dirichletBC_set <- TRUE
      }
      if(private$is_parabolic & (!private$is_initialCondition_set)){
        stop("initialCondition must be provided.")
      }
      private$cpp_handler_$solve()
      private$DifferentialOp_$Function()$set_coefficients(private$cpp_handler_$solution())
      invisible(self)
    },
    dofs_coordinates = function(){
      private$cpp_handler_$dofs_coordinates()
    },
    set_boundary_condition = function(boundary_condition, type="dirichlet", on=NULL){
      if(!any(type == c("dirichlet", "Dirichlet", "d"))) 
        stop("Only Dirichlet boundary condtions allowed.")
      dirichletBC_ <- NULL
      if(any(typeof(boundary_condition) == c("function", "closure"))){
        if(!private$is_parabolic){
          dirichletBC_ <- as.matrix(boundary_condition(private$cpp_handler_$dofs_coordinates()))
        }else{
          times <- private$DifferentialOp_$Function()$FunctionSpace()$mesh()$time_nodes()
          dirichletBC_ <- boundary_condition(private$cpp_handler_$dofs_coordinates(), times)
        }
      }else if(any(typeof(boundary_condition) == c("matrix","vector", "numeric" ,"double"))){
        if(nrow(as.matrix(boundary_condition)) == nrow(private$cpp_handler_$dofs_coordinates())){
          dirichletBC_ <- boundary_condition  
        }else if (nrow(as.matrix(boundary_condition)) == 1L){
          if(!private$is_parabolic)
            dirichletBC_ <- matrix(boundary_condition, nrow=nrow(private$cpp_handler_$dofs_coordinates()), ncol=1)
          else{
            times <- private$DifferentialOp_$Function()$FunctionSpace()$mesh()$time_nodes()
            dirichletBC_ <- matrix(boundary_condition, nrow=nrow(private$cpp_handler_$dofs_coordinates()), 
                                        ncol=length(times))
          }
        }   
        }
      private$is_dirichletBC_set <- TRUE
      private$cpp_handler_$set_dirichlet_bc(dirichletBC_)
      invisible(self)
    },
    set_initial_condition = function(initial_condition){
      if(!private$is_parabolic)
        stop("Cannot set initial condition for elliptic problem.")
      private$is_initialCondition_set <- TRUE
      if(typeof(initial_condition) == "closure"){ 
        private$cpp_handler_$set_initial_condition(initial_condition(private$cpp_handler_$dofs_coordinates()))
      }else if(any(typeof(initial_condition) == c("matrix","vector", "numeric" ,"double"))){
        if(nrow(as.matrix(initial_condition)) == nrow(private$cpp_handler_$dofs_coordinates())){
          private$cpp_handler_$set_initial_condition(initial_condition)
        }else if(nrow(as.matrix(initial_condition)) == 1L){
          private$cpp_handler_$set_initial_condition(matrix(initial_condition, nrow=nrow(private$cpp_handler_$dofs_coordinates()),
                                                 ncol=1))   
        }
      }
      invisible(self)
    },
    mass = function(){
      private$cpp_handler_$mass()
    },
    stiff = function(){
      private$cpp_handler_$stiff()
    }
  )
)

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
    pde_parameters[["time_nodes"]] <- as.vector(DifferentialOp$Function()$FunctionSpace()$mesh()$time_nodes())
  
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
  D <- extract_private(DifferentialOp$Function()$FunctionSpace()$mesh())$cpp_handler_
  fe_order <- DifferentialOp$Function()$FunctionSpace()$fe_order()
  cpp_handler <- NULL
  if (fe_order == 1) { ## linear finite elements
    cpp_handler <- new(cpp_pde_2d_fe1, D, pde_type - 1, pde_parameters)
  }
  if (fe_order == 2) { ## quadratic finite elements
    cpp_handler <- new(cpp_pde_2d_fe2, D, pde_type - 1, pde_parameters)
  }
  return(cpp_handler)
}

#' A PDEs object
#'
#' @param DifferentialOp a differential operator.
#' @param forcing a standard R function representing the forcing term of the PDE or a numeric value, in case of constant forcing term.
#' @return A R6 class representing a PDE.
#' @rdname pde
#' @export 
setGeneric("Pde", function(DifferentialOp, forcing) standardGeneric("Pde"))

#' @rdname pde
setMethod("Pde", signature=c(DifferentialOp="DifferentialOp", forcing="function"),
          function(DifferentialOp,forcing){
            
            pde_type = extract_pde_type(DifferentialOp)
            cpp_handler <- make_pde(DifferentialOp)
            
            quad_nodes <- cpp_handler$quadrature_nodes()
            ## evaluate forcing term on quadrature nodes
            if(pde_type == pde_type_list$parabolic) {
              times <- DifferentialOp$Function()$FunctionSpace()$mesh()$time_nodes()
              cpp_handler$set_forcing(as.matrix(forcing(quad_nodes, times)))
            }else{
              cpp_handler$set_forcing(as.matrix(forcing(quad_nodes)))
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
setMethod("Pde", signature=c(DifferentialOp="DifferentialOp", forcing="numeric"),
          function(DifferentialOp,forcing){
            pde_type = extract_pde_type(DifferentialOp)
            
            forcing_ <- NULL
            if(pde_type!=pde_type_list$parabolic){
              forcing_ <- function(points){
                return( matrix(forcing, nrow=nrow(points), ncol=1))
              }
            }else{
              forcing_ <- function(points, times){
                return( matrix(forcing, nrow=nrow(points), ncol=length(times)))
              }  
            }
            
            return(Pde(DifferentialOp,forcing_))
})
