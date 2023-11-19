# gradient of Function
.FunctionGradCtr <- setRefClass(
  Class = "FunctionGradObject",
  fields = c(
    f = "FunctionObject", ## Function of which the gradient is taken
    K = "ANY" ## set by product operator
  )
)

## take gradient of Function

#' compute gradient of FunctionObject
#'
#' @param f a FunctionObject created by \code{Function}:
#' @return An S4 object representing the gradient of the FunctionObject passed as parameter.
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
setMethod("grad", signature(f = "FunctionObject"), function(f) {
  .FunctionGradCtr(f = f)
})

#' product overload for FunctionGradObejct
#'
#' @param e1 a numeric matrix.
#' @param e2 a FunctionGradObject created by \code{grad} function.
#' @return A FunctionGradObject.
#' @export 
setMethod("*", signature=c(e1="matrix", e2="FunctionGradObject"),
          function(e1,e2){
        .FunctionGradCtr(f = e2$f, K = e1)    
})

# ## diffusion tensor - FunctionGrad product overload
# `*.FunctionGradObject` <- function(op1, op2) {
#   if (!is.matrix(op1) && !is.function(op1))
#     stop("bad diffusion tensor type")
#   .FunctionGradCtr(f = op2$f, K = op1)
# }

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

#' plus operator overload for DiffOpObject
#'
#' @param e1 a DiffOpObject.
#' @param e2 a DiffOpObject.
#' @return A S4 object representing the sum of two differential operators.
#' @export 
setMethod("+", signature = c(e1="DiffOpObject", e2="DiffOpObject"),
          function(e1, e2) {
            if (tracemem(e1$f) != tracemem(e2$f)) {
              stop("operator arguments must be the same")
            }
            ## check for duplicated operator tokens
            if (any(duplicated(c(e1$tokens, e2$tokens)))) {
              stop("duplicated operator")
            }
            .DiffOpCtr(
              tokens = c(e1$tokens, e2$tokens),
              params = c(e1$params, e2$params),
              f = e1$f
            )
          }
)

#' difference of differential operators
#'
#' @param e1 a DiffOpObject.
#' @param e2 a DiffOpObject.
#' @return A S4 object representing the difference of two differential operators.
#' @rdname minus_DiffOb_op
#' @export 
setMethod("-", signature = c(e1="DiffOpObject", e2="DiffOpObject"),
          function(e1, e2) {
            if (tracemem(e1$f) != tracemem(e2$f)) {
              stop("operator arguments must be the same")
            }
            ## check for duplicated operator tokens
            if (any(duplicated(c(e1$tokens, e2$tokens)))) {
              stop("duplicated operator")
            }
            for(i in 1:length(e2$params)){
              e2$params[[i]] = -e2$params[[i]] 
            }
            .DiffOpCtr(
              tokens = c(e1$tokens, e2$tokens),
              params = c(e1$params, e2$params),
              f = e1$f
            )
})

 
# minus (unary) operator for DiffOpObject
# 
# @param e1 a DiffOpObject.
# @return A S4 object differential operator whose paramter has changed sign.
# @export

#' @rdname minus_DiffOb_op
setMethod("-", signature =c(e1 = "DiffOpObject", e2 = "missing"),
          function(e1, e2){
            e1$params[[1]] <- -e1$params[[1]]
            e1
          }
)

#' product by scalar for differential operators
#'
#' @param e1 a numeric.
#' @param e2 a DiffOpObject.
#' @return A S4 object representing the product by scalar for a differential operator.
#' @export 
setMethod("*", signature=c(e1="numeric", e2="DiffOpObject"),
          function(e1,e2){
            #if (!is.numeric(e1)) stop("bad product")
            e2$params[[1]] <- e1*e2$params[[1]]
            e2
            
})

## diffusion term
.DiffusionCtr <- setRefClass(
    Class = "DiffusionOperator",
    contains = "DiffOpObject"
)

## laplace() returns a special operator for the case of
## isotropic  and stationary diffusion

#' laplace operator for FunctionObject
#'
#' @param f a FunctionObject.
#' @return A S4 object representing the laplace operator applied to the function passed as parameter.
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

#' divergence operator FunctionGradObject
#'
#' @param f a FunctionObject.
#' @return A S4 object representing the diffusion term of a second order linear differential operator.
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

#' dot product between vector and FunctionGradObject
#'
#' @param op1 a numeric vector.
#' @param op2 a FunctionGradObject.
#' @return A S4 object representing the advection term of a second order linear differential operator.
#' @rdname dot_product
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

#' @rdname dot_product
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

#' product by scalar for FunctionObejct
#'
#' @param e1 a numeric.
#' @param e2 a FunctioObject created by \code{Function}.
#' @return A S4 object representing the reaction term of a second order linear differential operator.
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
setMethod("*", signature = c(e1="numeric", e2="FunctionObject"),
          function(e1,e2){
            .ReactionCtr(
              tokens = "reaction",
              params = list(reaction = e1),
              f = e2
            )
})

.TimeDerivativeCtr <- setRefClass(
  Class = "TimeDerivative",
  contains ="DiffOpObject"
)

# overloading stats::dt function

# time derivate of FunctionObejct
#
# @param x a FunctionObject.
# @return A S4 object representing the time derivative of a FunctionObject.
# @export
# @examples
# \dontrun{
# library(femR)
# data("unit_square")
# mesh <- Mesh(unit_square)
# Vh <- FunctionSpace(mesh)
# f <- Function(Vh)
# dt(f)
# }
# dt.FunctionObject <- function(x){
#   .TimeDerivativeCtr(
#     tokens="time",
#     params = list(time=1L),
#     f=x
#   )
# }

#' time derivate of FunctionObejct
#'
#' @param x a FunctionObject.
#' @param df missing.
#' @param ncp missing.
#' @return A S4 object representing the time derivative of a FunctionObject.
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
setMethod("dt", signature = c(x="FunctionObject", df="missing", ncp="missing"),
           function(x,df,ncp){
             .TimeDerivativeCtr(
               tokens="time",
               params = list(time=1L),
              f=x
            )
})
