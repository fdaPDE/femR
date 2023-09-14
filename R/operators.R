## finite element function
setClass(
    "Function",
    slots = c(
        mesh  = "list", 
        coeff = "matrix"
    )
)
## constructor
Function <- function(domain) {
    coeff = matrix(ncol = 1, nrow = 0)
    new("Function", coeff = coeff, mesh = domain)
}

## gradient of Function
setClass(
    "FunctionGrad",
    slots = c(
        ## Function of which the gradient is taken
        f = "Function",
        K = "ANY" ## set by product operator
    )
)

## take gradient of Function
setGeneric("grad", function(f) standardGeneric("grad"))
setMethod("grad", signature(f = "Function"), function(f) {
    new("FunctionGrad", f = f)
})
## diffusion tensor - FunctionGrad product overload
`*.FunctionGrad` <- function(op1, op2) {
    if (!is.matrix(op1) && !is.function(op1))
        stop("bad diffusion tensor tyoe")
    new("FunctionGrad", f = op2@f, K = op1)
}

## base class for differential operators
setClass(
    "diff_op",
    slots = c(
        ## single blocks composing the overall operator
        tokens = "vector",
        params = "list",
        ## Function to which the operator is applied
        f = "Function"
    )
)
## sum of differential operators
`+.diff_op` <- function(op1, op2) {
    if (tracemem(op1@f) != tracemem(op2@f)) {
        stop("operator arguments must be the same")
    }
    ## check for duplicated operator tokens
    if (any(duplicated(c(op1@tokens, op2@tokens)))) {
        stop("duplicated operator")
    }
    new("diff_op",
        tokens = c(op1@tokens, op2@tokens),
        params = c(op1@params, op2@params),
        f = op1@f)
}
## differential operator product by scalar
`*.diff_op` <- function(op1, op2){
    if (!is.numeric(op1)) stop("bad product")
    op2@params <- op1*op2@params
    op2
}

## diffusion term
setClass(
    "diffusion",
    contains = "diff_op"
)
## laplace() returns a special operator for the case of
## isotropic  and stationary diffusion
laplace <- function(f) {
    if (!is(f, "Function")) {
        stop("wrong argument type")
    }
    new("diffusion",
        tokens = "diffusion",
        params = list(diffusion = 1),
        f = f)
}
## the general non-isotrpic, non-stationary diffusion operator
div <- function(f) {
    if (is(f, "FunctionGrad")) {
        if (!is.null(f@K)) {
            return(new("diffusion",
                       tokens = "diffusion",
                       params = list(diffusion = f@K),
                       f = f@f)
                   )
        }
    }
    stop("wrong argument type")
}

## transport term
setClass(
    "transport",
    contains = "diff_op"
)
setGeneric("dot", function(op1, op2) standardGeneric("dot"))
setMethod("dot", signature(op1 = "vector", op2 = "FunctionGrad"), function(op1, op2) {
    new("transport",
        tokens = "transport",
        params = list(transport = as.matrix(op1)),
        f = op2@f)
})

## reaction term
setClass(
    "reaction",
    contains = "diff_op"
)
`*.Function` <- function(c, f) {
    if (!is.function(c) && !is.numeric(c)) {
        stop("wrong argument type")
    }
    new("reaction",
        tokens = "reaction",
        params = list(reaction = c),
        f = f)
}
