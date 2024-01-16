# returns the private attributes of the R6 class passed as parameter 
setGeneric("extract_private", function(x) standardGeneric("extract_private"))
setMethod("extract_private", signature = "R6", 
          function(x){
            x$.__enclos_env__$private
})

setGeneric("set_private", function(x, attribute, value) standardGeneric("set_private"))
setMethod("set_private", signature = c("R6", "character"), 
          function(x, attribute, value){
            invisible(x$.__enclos_env__$private[[attribute]] <- value)
})
