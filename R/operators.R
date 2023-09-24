## take gradient of Function
setGeneric("grad", function(f) standardGeneric("grad"))
setMethod("grad", signature(f = "feFunctionObject"), function(f) {
  .feFunctionGradCtr(f = f)
})

## diffusion tensor - FunctionGrad product overload
`*.feFunctionGradObject` <- function(op1, op2) {
  if (!is.matrix(op1) && !is.function(op1))
    stop("bad diffusion tensor type")
  .feFunctionGradCtr(f = op2$f, K = op1)
}

## base class for differential operators
.DiffOpCtr <- setRefClass(
    Class = "DiffOpObject",
    fields = c(
        ## single blocks composing the overall operator
        tokens = "vector",
        params = "list",
        ## Function to which the operator is applied
        f = "feFunctionObject"
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

## differential operators minus (unary) operator
`-.DiffOpObject` <- function(op) {
  op$params[[1]] <- -op$params[[1]]
  op
}

## differential operator product by scalar
`*.DiffOpObject` <- function(op1, op2){
    if (!is.numeric(op1)) stop("bad product")
    op2$params[[1]] <- op1*op2$params[[1]]
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
    if (!is(f, "feFunctionObject")) {
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
    if (is(f, "feFunctionGradObject")) {
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
setMethod("dot", signature(op1 = "vector", op2 = "feFunctionGradObject"),
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
`*.feFunctionObject` <- function(c, f) {
    if (!is.function(c) && !is.numeric(c)) {
        stop("wrong argument type")
    }
    .ReactionCtr(
        tokens = "reaction",
        params = list(reaction = c),
        f = f
    )
}
