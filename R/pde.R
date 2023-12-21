
pde_type_list <- list("laplacian" = 1, "elliptic" = 2, "parabolic" = 3)

# Pde Class wraps C++ R_PDE class
.PdeCtr <- setRefClass(
  Class = "PdeObject",
  fields = c(
    is_dirichletBC_set = "logical",      
    is_initialCondition_set = "logical",
    is_parabolic ="logical",          
    cpp_handler = "ANY"),
  methods = c(
    solve = function(){
      if(!is_dirichletBC_set){
        if(!is_parabolic){
          dirichletBC_ = as.matrix(rep(0,times=nrow(cpp_handler$get_dofs_coordinates())))
        }else 
          dirichletBC_ = matrix(0, nrow=nrow(cpp_handler$get_dofs_coordinates()), ncol=length(times))
      cpp_handler$set_dirichlet_bc(dirichletBC_)
      is_dirichletBC_set <<- TRUE
      }
      if(is_parabolic & (!is_initialCondition_set)){
        stop("initialCondition must be provided.")
      }
      cpp_handler$solve()
    },
    get_solution = function(){
      cpp_handler$solution()
    },
    get_dofs_coordinates = function(){
      cpp_handler$get_dofs_coordinates()
    },
    set_boundary_condition = function(fun, type="dirichlet", on=NULL){
      if(!any(type == c("dirichlet", "Dirichlet", "d"))) 
        stop("Only Dirichlet boundary condtions allowed.")
      dirichletBC_ <- NULL
      if(typeof(fun) == "closure"){
        if(!is_parabolic){
          dirichletBC_ <- as.matrix(fun(cpp_handler$get_dofs_coordinates()))
        }else{
          dirichletBC_ <- fun(cpp_handler$get_dofs_coordinates(),times)
        }
      }else if(any(typeof(fun) == c("matrix","vector", "numeric" ,"double"))){
        if(nrow(as.matrix(fun)) == nrow(cpp_handler$get_dofs_coordinates())){
          dirichletBC_ <- fun  
        }else if (nrow(as.matrix(fun)) == 1L){
          if(!is_parabolic)
            dirichletBC_ <- matrix(fun, nrow=nrow(cpp_handler$get_dofs_coordinates()), ncol=1)
          else 
            dirichletBC_ <- matrix(fun, nrow=nrow(cpp_handler$get_dofs_coordinates()), ncol=length(times))
          }   
        }
      is_dirichletBC_set <<- TRUE
      cpp_handler$set_dirichlet_bc(dirichletBC_)
    },
    set_initial_condition = function(fun){
      if(!is_parabolic)
        stop("Cannot set initial condition for elliptic problem.")
      is_initialCondition_set <<- TRUE
      if(typeof(fun) == "closure"){ 
        cpp_handler$set_initial_condition(fun(cpp_handler$get_dofs_coordinates()))
      }else if(any(typeof(fun) == c("matrix","vector", "numeric" ,"double"))){
        if(nrow(as.matrix(fun)) == nrow(cpp_handler$get_dofs_coordinates())){
          cpp_handler$set_initial_condition(fun)
        }else if(nrow(as.matrix(fun)) == 1L){
          cpp_handler$set_initial_condition(matrix(fun, nrow=nrow(cpp_handler$get_dofs_coordinates()),
                                                 ncol=1))   
        }
      }
    },
    get_mass = function(){
      cpp_handler$mass()
    },
    get_stiff = function(){
      cpp_handler$stiff()
    }
  )
)

## infers the type of a pde
extract_pde_type <- function(L) {
  if ("time" %in% names(L$params)) {
    return(pde_type_list$parabolic)
  }
  if ("diffusion" %in% names(L$params) && !is.matrix(L$params$diffusion)) {
    return(pde_type_list$laplacian)
  }
  return(pde_type_list$elliptic)
}

## parse pde parameters
parse_pde_parameters <- function(L){
  pde_parameters <- NULL
  pde_parameters$diffusion <- 1.0
  pde_parameters$transport <- matrix(0, nrow = 2, ncol = 1)
  pde_parameters$reaction  <- 0.0
  
  for (i in length(L$params)) {
    pde_parameters[[paste(L$tokens[i])]] <- L$params[[paste(L$tokens[i])]]
  }
  return(pde_parameters)
}

## build cpp object
make_pde <- function(L) {
  pde_type = extract_pde_type(L)
  
  ## pde parameters
  pde_parameters <- parse_pde_parameters(L)
  
  if(pde_type == pde_type_list$parabolic) pde_parameters["time_mesh"] <- L$f$FunctionSpace$mesh$times
  
  ## define Rcpp module
  D <- L$f$FunctionSpace$mesh$cpp_handler ## domain
  fe_order <- L$f$FunctionSpace$fe_order  ## finite element order
  cpp_handler <- NULL
  if (fe_order == 1) { ## linear finite elements
    cpp_handler <- new(cpp_pde_2d_fe1, D, pde_type - 1, pde_parameters, L$f)
  }
  if (fe_order == 2) { ## quadratic finite elements
    cpp_handler <- new(cpp_pde_2d_fe2, D, pde_type - 1, pde_parameters, L$f)
  }
  
  return(cpp_handler)
}

#' A PDEs object
#'
#' @param L a differential operator.
#' @param f a standard R function representing the forcing term of the PDE or a numeric value, in case of constant forcing term.
#' @return A S4 object representing a PDE.
#' @rdname pde
#' @export 
setGeneric("Pde", function(L,f) standardGeneric("Pde"))

#' @rdname pde
setMethod("Pde", signature=c(L="DiffOpObject", f="ANY"),
          function(L,f){
            
            pde_type = extract_pde_type(L)
            cpp_handler <- make_pde(L)
            
            quad_nodes <- cpp_handler$get_quadrature_nodes()
            ## evaluate forcing term on quadrature nodes
            if(pde_type == pde_type_list$parabolic) {
              times <- L$f$FunctionSpace$mesh$times
              cpp_handler$set_forcing(as.matrix(f(quad_nodes, times)))
            }else{
              cpp_handler$set_forcing(as.matrix(f(quad_nodes)))
            }
            
            ## initialize solver 
            cpp_handler$init()
            
            is_parabolic = pde_type == pde_type_list$parabolic
            
            ## return
            .PdeCtr(is_dirichletBC_set = FALSE,
                    is_initialCondition_set = FALSE,
                    is_parabolic = is_parabolic,
                    cpp_handler = cpp_handler)
})

#' @rdname pde
setMethod("Pde", signature=c(L="DiffOpObject", f="numeric"),
          function(L,f){
            pde_type = extract_pde_type(L)
            
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
            
            return(Pde(L,fun))
})
