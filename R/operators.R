## finite element function
.FunctionCtr <- setRefClass(
    Class = "FunctionObject",
    fields = c(
        mesh  = "ANY", 
        coeff = "matrix",
        pde = "ANY"
    ),
    methods = list(
        eval_at = function(X) {
            M = dim(mesh$nodes)[2]
            if(is.vector(X)) {
                pde$eval(coeff, as.matrix(t(X)))
            } else {           
            if(dim(X)[2] != M) {
                stop(paste("matrix of evaluation points should be an N x", M, "matrix"))
            }
            pde$eval(coeff, as.matrix(X))
            }
        }
    )
    
)
## constructor
Function <- function(domain) {
    coeff = matrix(ncol = 1, nrow = 0)
    .FunctionCtr(coeff = coeff, mesh = domain)
}

## gradient of Function
.FunctionGradCtr <- setRefClass(
    Class = "FunctionGradObject",
    fields = c(
        f = "FunctionObject", ## Function of which the gradient is taken
        K = "ANY" ## set by product operator
    )
)

## take gradient of Function
setGeneric("grad", function(f) standardGeneric("grad"))
setMethod("grad", signature(f = "FunctionObject"), function(f) {
    .FunctionGradCtr(f = f)
})
## diffusion tensor - FunctionGrad product overload
`*.FunctionGradObject` <- function(op1, op2) {
    if (!is.matrix(op1) && !is.function(op1))
        stop("bad diffusion tensor type")
    .FunctionGradCtr("FunctionGrad", f = op2$f, K = op1)
}

## base class for differential operators
.DiffOpCtr <- setRefClass(
    Class = "DiffOpObject",
    fields = c(
        ## single blocks composing the overall operator
        tokens = "vector",
        params = "list",
        ## Function to which the operator is applied
        f = "FunctionObject"
    )
)
## sum of differential operators
`+.DiffOpObject` <- function(op1, op2) {
    if (tracemem(op1$f) != tracemem(op2$f)) {
        stop("operator arguments must be the same")
    }
    ## check for duplicated operator tokens
    if (any(duplicated(c(op1$tokens, op2$tokens)))) {
        stop("duplicated operator")
    }
    .DiffOpCtr(
        tokens = c(op1$tokens, op2$tokens),
        params = c(op1$params, op2$params),
        f = op1$f
    )
}
## differential operator product by scalar
`*.DiffOpObject` <- function(op1, op2){
    if (!is.numeric(op1)) stop("bad product")
    op2$params <- op1*op2$params
    op2
}

## diffusion term
.DiffusionCtr <- setRefClass(
    Class = "DiffusionOperator",
    contains = "DiffOpObject"
)
## laplace() returns a special operator for the case of
## isotropic  and stationary diffusion
laplace <- function(f) {
    if (!is(f, "FunctionObject")) {
        stop("wrong argument type")
    }
    .DiffusionCtr(
        tokens = "diffusion",
        params = list(diffusion = 1),
        f = f
    )
}
## the general non-isotrpic, non-stationary diffusion operator
div <- function(f) {
    if (is(f, "FunctionGradObject")) {
        if (!is.null(f$K)) {
            return(.DiffusionCtr(
                       tokens = "diffusion",
                       params = list(diffusion = f$K),
                       f = f$f)
                   )
        }
    }
    stop("wrong argument type")
}

## transport term
.TransportCtr <- setRefClass(
    Class = "TransportOperator",
    contains = "DiffOpObject"
)
setGeneric("dot", function(op1, op2) standardGeneric("dot"))
setMethod("dot", signature(op1 = "vector", op2 = "FunctionGradObject"),
          function(op1, op2) {
              .TransportCtr(
                  tokens = "transport",
                  params = list(transport = as.matrix(op1)),
                  f = op2$f
              )
})

## reaction term
.ReactionCtr <- setRefClass(
    Class = "ReactionOperator",
    contains = "DiffOpObject"
)
`*.FunctionObject` <- function(c, f) {
    if (!is.function(c) && !is.numeric(c)) {
        stop("wrong argument type")
    }
    .ReactionCtr(
        tokens = "reaction",
        params = list(reaction = c),
        f = f
    )
}
