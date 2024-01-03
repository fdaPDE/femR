# gradient of Function
.FunctionGradCtr <- R6Class("FunctionGrad",
  public = list(
    f = "Function", ## Function of which the gradient is taken
    K = "ANY",      ## set by product operator
    initialize = function(f, K=matrix(c(1,0,0,1), nrow=2,ncol=2,byrow=T)){
      self$f <- f
      self$K <- K
    }
  )
)

#' Function gradient
#'
#' @name grad
#'
#' @exportClass FunctionGrad
setOldClass(c("FunctionGrad", "R6"))

## take gradient of Function

#' compute gradient of Function
#'
#' @param f a Function created by \code{Function}:
#' @return An S4 object representing the gradient of the Function passed as parameter.
#' @export 
#' @rdname grad
#' @examples
#' \dontrun{
#' library(femR)
#' data("unit_square")
#' mesh <- Mesh(unit_square)
#' Vh <- FunctionSpace(mesh)
#' f <- Function(Vh)
#' grad_f <- grad(f)
#' }
setGeneric("grad", function(f) standardGeneric("grad"))

#' @rdname grad
setMethod("grad", signature(f = "Function"), function(f) {
  .FunctionGradCtr$new(f = f)
})

# setMethod("*", signature=c(e1="matrix", e2="FunctionGrad")

#' product overload for FunctionGradObejct
#'
#' @param e1 a numeric matrix.
#' @param e2 a FunctionGrad class created by \code{grad} function.
#' @return A FunctionGrad class.
#' @export 
`*.FunctionGrad` <- function(e1, e2){
  if (!is.matrix(e1)) stop("First argument must be a matrix!")
        .FunctionGradCtr$new(f = e2$f, K = e1)    
}

## base class for differential operators
.DiffOpCtr <- R6Class("DiffOpObject",
    public = list(
        initialize = function(tokens, params, f){
          self$tokens <- tokens
          self$params <- params
          self$f <- f
        },
        ## single blocks composing the overall operator
        tokens = "vector",
        params = "list",
        ## Function to which the operator is applied
        f = "Function"
    )
)

#' R6 class representing a differential operator
#' 
#' @name Differential Operator
#'
#' @exportClass DiffOpObject
setOldClass(c("DiffOpObject", "R6"))

##setMethod("+", signature = c(e1="DiffOpObject", e2="DiffOpObject"),

## sum of differential operators

#' plus operator overload for DiffOpObject
#'
#' @param e1 a DiffOpObject.
#' @param e2 a DiffOpObject.
#' @return A S4 object representing the sum of two differential operators.
#' @name DifferentialOps
#' @export 
`+.DiffOpObject` <- function(e1, e2){
            if (!identical(e1$f, e2$f)) {
              stop("operator arguments must be the same")
            }
            ## check for duplicated operator tokens
            if (any(duplicated(c(e1$tokens, e2$tokens)))) {
              stop("duplicated operator")
            }
            .DiffOpCtr$new(
              tokens = c(e1$tokens, e2$tokens),
              params = c(e1$params, e2$params),
              f = e1$f
            )
}

# setMethod("-", signature = c(e1="DiffOpObject", e2="DiffOpObject"),

#' difference of differential operators
#'
#' @param e1 a DiffOpObject.
#' @param e2 a DiffOpObject.
#' @return A S4 object representing the difference of two differential operators.
#' @rdname minus_DiffOb_op
#' @name DifferentialOps
#' @export 
`-.DiffOpObject` <- function(e1, e2){
            if(missing(e2)){
              e1$params[[1]] <- -e1$params[[1]]
              return(e1)
            }else if (!identical(e1$f, e2$f)) {
              stop("operator arguments must be the same")
            }
            ## check for duplicated operator tokens
            if (any(duplicated(c(e1$tokens, e2$tokens)))) {
              stop("duplicated operator")
            }
            for(i in 1:length(e2$params)){
              e2$params[[i]] = -e2$params[[i]] 
            }
            .DiffOpCtr$new(
              tokens = c(e1$tokens, e2$tokens),
              params = c(e1$params, e2$params),
              f = e1$f
            )
  }

#setMethod("-", signature =c(e1 = "DiffOpObject", e2 = "missing"),

# minus (unary) operator for differential operators
# 
# @param e1 a DiffOpObject.
# @return A S4 object differential operator whose paramter has changed sign.
# @export

# @rdname minus_DiffOb_op
# `-.DiffOpObject` <- function(e1){
#             e1$params[[1]] <- -e1$params[[1]]
#             e1
# }

# setMethod("*", signature=c(e1="numeric", e2="DiffOpObject"),

#' product by scalar for differential operators
#'
#' @param e1 a numeric.
#' @param e2 a DiffOpObject.
#' @return A S4 object representing the product by scalar for a differential operator.
#' @name DifferentialOps
#' @export 
`*.DiffOpObject` <- function(e1,e2){  
            if (!is.numeric(e1)) stop("bad product")
            e2$params[[1]] <- e1*e2$params[[1]]
            e2
          }
  
## diffusion term
.DiffusionCtr <- R6Class("DiffusionOperator",
    inherit = .DiffOpCtr
)

#' Diffusion Operator
#'
#' @name DifferentialOperator
#'
#' @exportClass DiffusionOperator
setOldClass(c("DiffusionOperator", "DiffOpObject"))

## laplace() returns a special operator for the case of
## isotropic  and stationary diffusion

#' laplace operator for Function class
#'
#' @param f a Function.
#' @return A S4 object representing the laplace operator applied to the function passed as parameter.
#' @rdname DiffusionOperator
#' @export 
#' @examples
#' \dontrun{
#' library(femR)
#' data("unit_square")
#' mesh <- Mesh(unit_square)
#' Vh <- FunctionSpace(mesh)
#' f <- Function(Vh)
#' laplace_f <- laplace(f)
#' }
laplace <- function(f) {
    if (!is(f, "Function")) {
        stop("wrong argument type")
    }
    .DiffusionCtr$new(
        tokens = "diffusion",
        params = list(diffusion = 1),
        f = f
    )
}

## the general non-isotrpic, non-stationary diffusion operator

#' divergence operator FunctionGrad
#'
#' @param f a Function.
#' @return A S4 object representing the diffusion term of a second order linear differential operator.
#' @rdname DiffusionOperator
#' @export
#' @examples
#' \dontrun{
#' library(femR)
#' data("unit_square")
#' mesh <- Mesh(unit_square)
#' Vh <- FunctionSpace(mesh)
#' f <- Function(Vh)
#' K <- matrix(c(1,2,1,0),nrow=2,ncol=2)
#' diffusion <- div(K*grad(f))
#' } 
div <- function(f) {
    if (is(f, "FunctionGrad")) {
        if (!is.null(f$K)) {
            return(.DiffusionCtr$new(
                       tokens = "diffusion",
                       params = list(diffusion = f$K),
                       f = f$f)
                   )
        }
    }
    stop("wrong argument type")
}

## transport term
.TransportCtr <- R6Class("TransportOperator",
    inherit = .DiffOpCtr
)

#' Transport Operator
#'
#' @name TransportOperator
#'
#' @exportClass TransportOperator
setOldClass(c("TransportOperator", "DiffOpObject"))

#' dot product between vector and FunctionGrad
#'
#' @param op1 a numeric vector.
#' @param op2 a FunctionGrad.
#' @return A S4 object representing the advection term of a second order linear differential operator.
#' @rdname TransportOperator
#' @export
#' @examples
#' \dontrun{
#' library(femR)
#' data("unit_square")
#' mesh <- Mesh(unit_square)
#' Vh <- FunctionSpace(mesh)
#' f <- Function(Vh)
#' b <- c(1,1)
#' advection <- dot(b,grad(f))
#' } 
setGeneric("dot", function(op1, op2) standardGeneric("dot"))

#' @rdname TransportOperator
setMethod("dot", signature(op1 = "vector", op2 = "FunctionGrad"),
          function(op1, op2) {
              .TransportCtr$new(
                  tokens = "transport",
                  params = list(transport = as.matrix(op1)),
                  f = op2$f
              )
})

## reaction term
.ReactionCtr <- R6Class("ReactionOperator",
    inherit = .DiffOpCtr
)

#' Reaction Operator
#'
#' @name ReactionOperator
#'
#' @exportClass ReactionOperator
setOldClass(c("ReactionOperator", "DiffOpObject"))

# setMethod("*", signature = c(e1="numeric", e2="Function"),
          
#' product by scalar for FunctionObejct
#'
#' @param e1 a numeric.
#' @param e2 a FunctioObject created by \code{Function}.
#' @return A S4 object representing the reaction term of a second order linear differential operator.
#' @rdname ReactionOperator
#' @export 
#' @examples
#' \dontrun{
#' library(femR)
#' data("unit_square")
#' mesh <- Mesh(unit_square)
#' Vh <- FunctionSpace(mesh)
#' f <- Function(Vh)
#' reaction <- 2*f
#' }
`*.Function` <- function(e1,e2){
            if(!is(e1,"numeric")) stop("First argument must be a scalar!")
            .ReactionCtr$new(
              tokens = "reaction",
              params = list(reaction = e1),
              f = e2
            )
}

.TimeDerivativeCtr <- R6Class("TimeDerivative",
  inherit = .DiffOpCtr
)

#' @name TimeDerivative
#'
#' @exportClass TimeDerivative
setOldClass(c("TimeDerivative", "DiffOpObject"))

# overloading stats::dt function

#' time derivate of FunctionObejct
#'
#' @param x a Function.
#' @param df missing.
#' @param ncp missing.
#' @return A S4 object representing the time derivative of a Function.
#' @rdname TimeDerivative
#' @export
#' @examples
#' \dontrun{
#' library(femR)
#' data("unit_square")
#' mesh <- Mesh(unit_square)
#' Vh <- FunctionSpace(mesh)
#' f <- Function(Vh)
#' dt(f)
#' }
setMethod("dt", signature = c(x="Function", df="missing", ncp="missing"),
           function(x,df,ncp){
             .TimeDerivativeCtr$new(
               tokens="time",
               params = list(time=1L),
              f=x
            )
})
