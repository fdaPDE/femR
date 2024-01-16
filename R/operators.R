# gradient of Function
.FunctionGradCtr <- R6Class("FunctionGrad",
  private = list(
    Function_ = "Function"
  ),
  public = list(
    K = "ANY",      ## set by product operator
    initialize = function(Function, K=matrix(c(1,0,0,1), nrow=2,ncol=2,byrow=T)){
      private$Function_ <- Function
      self$K <- K
    },
    Function = function(){
      private$Function_
    }
  )
)

setOldClass(c("FunctionGrad", "R6"))

## take gradient of Function

#' compute gradient of Function
#'
#' @param Function a Function created by \code{Function}:
#' @return A R6 object representing the gradient of the Function passed as parameter.
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
setGeneric("grad", function(Function) standardGeneric("grad"))

#' @rdname grad
setMethod("grad", signature(Function = "Function"), function(Function) {
  .FunctionGradCtr$new(Function = Function)
})

# #`*.FunctionGrad` <- function(e1, e2){ setMethod("*", signature=c(e1="matrix", e2="FunctionGrad"),

#' product overload for FunctionGradObejct
#'
#' @param e1 a numeric matrix.
#' @param e2 a FunctionGrad class created by \code{grad} function.
#' @return A FunctionGrad class.
#' @export 
`*.FunctionGrad` <- function(e1, e2){
  if (!is.matrix(e1)) stop("First argument must be a matrix!")
        .FunctionGradCtr$new(Function = e2$Function(), K = e1)    
}

## base class for differential operators
.DifferentialOpCtr <- R6Class("DifferentialOp",
    private = list(
      ## Function to which the operator is applied
      Function_ = "Function"
    ),
    public = list(
        initialize = function(tokens, params, Function){
          self$tokens <- tokens
          self$params <- params
          private$Function_ <- Function
        },
        ## single blocks composing the overall operator
        tokens = "vector",
        params = "list",
        Function = function(){
          private$Function_
        }
    )
)

setOldClass(c("DifferentialOp", "R6"))

##setMethod("+", signature = c(e1="DifferentialOp", e2="DifferentialOp"),
## `+.DifferentialOp` <- function(e1, e2){

## sum of differential operators

#' Sum of DifferentialOp
#'
#' @param e1 a \code{DifferentialOp}.
#' @param e2 a \code{DifferentialOp}.
#' @return A R6 class representing the sum of two differential operators.
#' @export 
`+.DifferentialOp` <- function(e1, e2){
            if (!identical(e1$Function(), e2$Function())) {
              stop("operator arguments must be the same")
            }
            ## check for duplicated operator tokens
            if (any(duplicated(c(e1$tokens, e2$tokens)))) {
              stop("duplicated operator")
            }
            .DifferentialOpCtr$new(
              tokens = c(e1$tokens, e2$tokens),
              params = c(e1$params, e2$params),
              Function = e1$Function()
            )
}

# setMethod("-", signature = c(e1="DifferentialOp", e2="DifferentialOp"),
# `-.DifferentialOp` <- function(e1, e2){

#' Difference of differential operators 
#'
#' @param e1 a \code{DifferentialOp}.
#' @param e2 a \code{DifferentialOp} or missing.
#' @return A R6 class representing the difference of two differential operators.
#' @export 
`-.DifferentialOp` <-function(e1 ,e2){
            if(missing(e2)){
              e1$params[[1]] <- -e1$params[[1]]
              return(e1)
            }else if (!identical(e1$Function(), e2$Function())) {
              stop("operator arguments must be the same")
            }
            if (!identical(e1$Function(), e2$Function())) 
              stop("operator arguments must be the same")
            
            ## check for duplicated operator tokens
            if (any(duplicated(c(e1$tokens, e2$tokens)))) {
              stop("duplicated operator")
            }
            for(i in 1:length(e2$params)){
              e2$params[[i]] = -e2$params[[i]] 
            }
            .DifferentialOpCtr$new(
              tokens = c(e1$tokens, e2$tokens),
              params = c(e1$params, e2$params),
              Function = e1$Function()
            )
}

#setMethod("-", signature =c(e1 = "DifferentialOp", e2 = "missing"),

# minus (unary) operator for differential operators
# 
# @param e1 a DifferentialOp.
# @return A R6 class differential operator whose paramter has changed sign.
# @name DifferentialOps
# @export
# setMethod("-", signature =c(e1 = "DifferentialOp", e2 = "missing"),
#           function(e1 ,e2){
#             e1$params[[1]] <- -e1$params[[1]]
#             e1
# })

# @rdname minus_DiffOb_op
# `-.DifferentialOp` <- function(e1){
#             e1$params[[1]] <- -e1$params[[1]]
#             e1
# }

# setMethod("*", signature=c(e1="numeric", e2="DifferentialOp"),
# `*.DifferentialOp` <- function(e1,e2){  

#' product by scalar for differential operators
#'
#' @param e1 a numeric.
#' @param e2 a DifferentialOp.
#' @return A R6 class representing the product by scalar for a differential operator.
#' @name times_DifferentialOps
#' @export 
`*.DifferentialOp` <- function(e1,e2){
            if (!is.numeric(e1)) stop("bad product")
            e2$params[[1]] <- e1*e2$params[[1]]
            e2
}
  
## diffusion term
.DiffusionCtr <- R6Class("DiffusionOperator",
    inherit = .DifferentialOpCtr
)

setOldClass(c("DiffusionOperator", "DifferentialOp"))

## laplace() returns a special operator for the case of
## isotropic  and stationary diffusion

#' laplace operator for Function class
#'
#' @param Function a Function.
#' @return A R6 class representing the laplace operator applied to the function passed as parameter.
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
laplace <- function(Function) {
    if (!is(Function, "Function")) {
        stop("wrong argument type")
    }
    .DiffusionCtr$new(
        tokens = "diffusion",
        params = list(diffusion = 1),
        Function = Function
    )
}

## the general non-isotrpic, non-stationary diffusion operator

#' divergence operator FunctionGrad
#'
#' @param Function a Function.
#' @return A R6 class representing the diffusion term of a second order linear differential operator.
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
div <- function(Function) {
    if (is(Function, "FunctionGrad")) {
        if (!is.null(Function$K)) {
            return(.DiffusionCtr$new(
                       tokens = "diffusion",
                       params = list(diffusion = Function$K),
                       Function = Function$Function())
                   )
        }
    }
    stop("wrong argument type")
}

## transport term
.TransportCtr <- R6Class("TransportOperator",
    inherit = .DifferentialOpCtr
)

setOldClass(c("TransportOperator", "DifferentialOp"))

#' dot product between vector and FunctionGrad
#'
#' @param op1 a numeric vector.
#' @param op2 a FunctionGrad.
#' @return A R6 class representing the advection term of a second order linear differential operator.
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
                  Function = op2$Function()
              )
})

## reaction term
.ReactionCtr <- R6Class("ReactionOperator",
    inherit = .DifferentialOpCtr
)

setOldClass(c("ReactionOperator", "DifferentialOp"))

# setMethod("*", signature = c(e1="numeric", e2="Function"),
# `*.Function` <- function(e1,e2){

#' product by scalar for FunctionObejct
#'
#' @param e1 a numeric.
#' @param e2 a FunctioObject created by \code{Function}.
#' @return A R6 class representing the reaction term of a second order linear differential operator.
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
              Function = e2
            )
}

.TimeDerivativeCtr <- R6Class("TimeDerivative",
  inherit = .DifferentialOpCtr
)

setOldClass(c("TimeDerivative", "DifferentialOp"))

# overloading stats::dt function

#' time derivate of FunctionObejct
#'
#' @param x a Function.
#' @param df missing.
#' @param ncp missing.
#' @return A R6 class representing the time derivative of a Function.
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
               tokens = "time",
               params = list(time=1L),
               Function = x
            )
})
