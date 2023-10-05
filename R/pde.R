## create pde object backed by Cpp_pde_module
pde <- function(L, u, dirichletBC) {
    D = L$f$FunctionalSpace$mesh ## domain object

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
    
    ## prepare pde_parameters list
    pde_parameters$transport <- matrix(0, nrow = 2, ncol = 1)
    pde_parameters$reaction  <- 0.0
    for(i in 1:length(L$tokens)) {
      pde_parameters[[paste(L$tokens[i])]] <- L$params[[paste(L$tokens[i])]]
    }
    
    fe_order <- L$f$FunctionalSpace$fe_order
    ## define Rcpp module
    pde_ <- NULL
    if (fe_order == 1) { ## linear finite elements
        pde_ <- new(PDE_2D_ORDER_1, D, pde_type, pde_parameters, L$f)
    }
    if (fe_order == 2) { ## quadratic finite elements
        pde_ <- new(PDE_2D_ORDER_2, D, pde_type, pde_parameters, L$f)
    }

    L$f$pde = pde_
    
    ## evaluate forcing term on quadrature nodes
    quad_nodes <- as.matrix(pde_$get_quadrature_nodes())
    pde_$set_forcing(as.matrix(u(quad_nodes)))
    
    ## initialize and return
    pde_$init()
    
    ## set Dirichlet boundary conditions
    pde_$set_dirichlet_bc(as.matrix(dirichletBC(pde_$get_dofs_coordinates())))
    
    return(pde_)
}
