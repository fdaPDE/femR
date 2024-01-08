# returns the private attributes of the R6 class passed as parameter 
setGeneric("extract_private", function(x) standardGeneric("extract_private"))
setMethod("extract_private", signature = "R6", 
          function(x){
            x$.__enclos_env__$private
})

setGeneric("set_geometry", function(x, value) standardGeneric("set_geometry"))
setMethod("set_geometry", signature = "Domain", 
          function(x, value){
            invisible(x$.__enclos_env__$private$geometry <- value)
})

setGeneric("set_time_interval", function(x, value) standardGeneric("set_time_interval"))
setMethod("set_time_interval", signature = "Domain", 
          function(x, value){
            invisible(x$.__enclos_env__$private$time_interval <- value)
})

setGeneric("set_crs", function(x, value) standardGeneric("set_crs"))
setMethod("set_crs", signature = "Domain", 
          function(x, value){
            invisible(x$.__enclos_env__$private$crs <- value)
})

setGeneric("set_times", function(x, value) standardGeneric("set_times"))
setMethod("set_times", signature = "Mesh", 
          function(x, value){
            invisible(x$.__enclos_env__$private$times <- value)
})

setGeneric("set_time_step", function(x, value) standardGeneric("set_time_step"))
setMethod("set_time_step", signature = "Mesh", 
          function(x, value){
            invisible(x$.__enclos_env__$private$time_step <- value)
})

# extract_private <- function(x){
#   if(!is(x,"R6")) stop("Input parameter must be an R6 class.")
#   x$.__enclos_env__$private
# }
