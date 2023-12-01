# Pde Class wraps C++ R_PDE class
.PdeCtr <- setRefClass(
  Class = "PdeObject",
  fields = c(
   # L  = "ANY",                         # DiffOpObject
    is_dirichletBC_set = "logical",      
    times =       "vector",              # vector of time instants (parabolic problem only)
    is_initialCondition_set = "logical", # initial condition       (parabolic problem only)
    pde_ = "ANY" ,                       # Wraps of C++ class
    is_parabolic = "logical",
    is_init = "logical"                  
  ),
  methods = c(
    solve = function(){
      if(!is_dirichletBC_set){
        if(!is_parabolic){
          dirichletBC_ = as.matrix(rep(0,times=nrow(pde_$get_dofs_coordinates())))
        }else 
          dirichletBC_ = matrix(0, nrow=nrow(pde_$get_dofs_coordinates()), ncol=length(times))
      pde_$set_dirichlet_bc(dirichletBC_)
      is_dirichletBC_set <<- TRUE
      }
      if(is_parabolic & (!is_initialCondition_set)){
        stop("initialCondition must be provided.")
      }
      if(!is_init) pde_$init()
      pde_$solve()
    },
    get_solution = function(){
      pde_$solution()
    },
    get_dofs_coordinates = function(){
      pde_$get_dofs_coordinates()
    },
    set_boundary_condition = function(fun, type="dirichlet", on=NULL){
      if(!any(type == c("dirichlet", "Dirichlet", "d"))) 
        stop("Only Dirichlet boundary condtions allowed.")
      dirichletBC_ <- NULL
      if(typeof(fun) == "closure"){
        if(!is_parabolic){
          dirichletBC_ <- as.matrix(fun(pde_$get_dofs_coordinates()))
        }else{
          dirichletBC_ <- fun(pde_$get_dofs_coordinates(),times)
        }
      }else if(any(typeof(fun) == c("matrix","vector", "numeric" ,"double")) & 
               nrow(as.matrix(fun)) == nrow(pde_$get_dofs_coordinates())){ # si controllerà "on"
        dirichletBC_ <- fun
      }
      is_dirichletBC_set <<- TRUE
      pde_$set_dirichlet_bc(dirichletBC_)
    },
    set_initial_condition = function(fun){
      if(!is_parabolic)
        stop("Cannot set initial condition for elliptic problem.")
      is_initialCondition_set <<- TRUE
      if(typeof(fun) == "closure"){ 
        pde_$set_initial_condition(fun(pde_$get_dofs_coordinates()))
      }else if(any(typeof(fun) == c("matrix","vector", "numeric" ,"double")) & 
              nrow(as.matrix(fun)) == nrow(pde_$get_dofs_coordinates())){
        pde_$set_initial_condition(fun)
      }
    },
    get_mass = function(){
      pde_$get_mass()
    },
    get_stiff = function(){
      pde_$get_stiff()
    }
  )
)

#' A PDEs object
#'
#' @param L a differential operator.
#' @param u a standard R function representing the forcing term of the PDE.
#' @return A S4 object representing a PDE.
#' @rdname pde
#' @export 
setGeneric("Pde", function(L,u) standardGeneric("Pde"))

#' @rdname pde
setMethod("Pde", signature=c(L="DiffOpObject", u="ANY"),
          function(L,u){
            D = L$f$FunctionSpace$mesh$data ## C++ R_Mesh class
            
            is_parabolic = FALSE
            times <- vector(mode="numeric", length=0L)
            if( length(L$f$FunctionSpace$mesh$times) != 0 ){ 
              is_parabolic = TRUE
              times <- L$f$FunctionSpace$mesh$times
            }
            
            ## set pde type
            pde_type <- 0
            pde_parameters <- NULL
            if ("diffusion" %in% names(L$params) & !is.matrix(L$params$diffusion)){
              ## specialized implementation for laplace operator
              pde_type <- 1L 
              
              pde_parameters$diffusion <- 0.0
              
            } else {
              ## general diffusion-transport-reaction problem, constant coefficients
              pde_type <- 2L
              
              pde_parameters$diffusion <- matrix(0, nrow = 2, ncol = 2)
              
            }
            
            if(is_parabolic){
              if(pde_type == 1L){
                pde_type = 3L
              }else{
                pde_type = 4L
              }
            }
            
            ## prepare pde_parameters list
            pde_parameters$transport <- matrix(0, nrow = 2, ncol = 1)
            pde_parameters$reaction  <- 0.0
            pde_parameters$time <- 0L
            for(i in 1:length(L$tokens)) {
              pde_parameters[[paste(L$tokens[i])]] <- L$params[[paste(L$tokens[i])]]
            }
            
            if(is_parabolic) pde_parameters$time <- times
            
            fe_order <- L$f$FunctionSpace$fe_order
            ## define Rcpp module
            pde_ <- NULL
            if (fe_order == 1) { ## linear finite elements
              pde_ <- new(PDE_2D_ORDER_1, D, pde_type, pde_parameters, L$f)
            }
            if (fe_order == 2) { ## quadratic finite elements
              pde_ <- new(PDE_2D_ORDER_2, D, pde_type, pde_parameters, L$f)
            }
            
            L$f$pde = pde_
            
            quad_nodes <- as.matrix(pde_$get_quadrature_nodes())
            ## evaluate forcing term on quadrature nodes
            if(!is_parabolic){
              pde_$set_forcing(as.matrix(u(quad_nodes)))
            }else{
              pde_$set_forcing(u(quad_nodes, times))
            }
            ## initialize solver 
            pde_$init()
            is_init = TRUE
            
            ## return
            .PdeCtr(times = times,
                    is_dirichletBC_set = FALSE,
                    is_initialCondition_set = FALSE,
                    pde_ = pde_,                  # Wraps of C++ class
                    is_parabolic = is_parabolic,
                    is_init=is_init)
})

