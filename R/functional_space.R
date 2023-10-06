.FunctionalSpaceCtr <- setRefClass(
  Class = "FunctionalSpaceObject",
  fields = c(
    mesh  = "ANY", 
    fe_order="integer"
  )
)

# constructor
setGeneric("FunctionalSpace", function(domain,fe_order) standardGeneric("FunctionalSpace"))
setMethod("FunctionalSpace", signature = c(domain="list", fe_order="numeric"),
          function(domain,fe_order){
            domain$elements <- domain$elements - 1
            storage.mode(domain$elements) <- "integer" 
            if(fe_order == 1){
              return(.FunctionalSpaceCtr(mesh=new(Mesh_2D, domain), fe_order=1L))
            }else if(fe_order == 2){
              return(.FunctionalSpaceCtr(mesh=new(Mesh_2D, domain), fe_order=2L))
            }
                                                  
})

## finite element function
.FunctionCtr <- setRefClass(
  Class = "FunctionObject",
  fields = c(
    FunctionalSpace  = "ANY", 
    coeff = "matrix",
    pde = "ANY"
  ),
  methods = list(
    eval_at = function(X) {
      M = dim(FunctionalSpace$mesh$nodes())[2]
      if(is.vector(X)) {
        pde$eval(FunctionalSpace$mesh, coeff, matrix(X, nrow=1,ncol=M))
      } else {           
        if(dim(X)[2] != M) {
          stop(paste("matrix of evaluation points should be an N x", M, "matrix"))
        }
        pde$eval(FunctionalSpace$mesh, coeff, as.matrix(X))
      }
    }
  )
)
## constructor
Function <- function(FunctionalSpace) {
  coeff = matrix(ncol = 1, nrow = 0)
  .FunctionCtr(coeff = coeff, FunctionalSpace = FunctionalSpace)
}

## FunctionObject plot overload
setMethod("plot", signature=c(x="FunctionObject"), function(x, ...){
  plot_data <- data.frame(X=x$FunctionalSpace$mesh$nodes()[,1], 
                          Y=x$FunctionalSpace$mesh$nodes()[,2],
                          coeff=x$coeff[1:nrow(x$FunctionalSpace$mesh$nodes())])
  I=x$FunctionalSpace$mesh$elements()[,1]
  J=x$FunctionalSpace$mesh$elements()[,2]
  K=x$FunctionalSpace$mesh$elements()[,3]
  plot_ly(plot_data, x=~X, y=~Y, z=~coeff,
          i = I, j = J, k = K,
          intensity=~coeff,color = ~coeff, type="mesh3d", ...) %>%
    layout(scene = list(
      aspectmode = "data", 
      xaxis = list(
        title = '',
        showgrid = F,
        zeroline = F,
        showticklabels = F),
      yaxis = list(
        title = '',
        showgrid = F,
        zeroline = F,
        showticklabels = F),
      zaxis = list(
        title = '',
        showgrid = F,
        zeroline = F,
        showticklabels = F)),
      camera = list(
        eye = list(x = 1.25, 
                   y = -1.25, 
                   z = 1.25))) %>%
    colorbar(len = 1, title="")
}
)

## FunctionObject contour overload
setMethod("contour", signature=c(x="FunctionObject"), function(x, ...){
  plot_data <- data.frame(X=x$FunctionalSpace$mesh$nodes()[,1], 
                          Y=x$FunctionalSpace$mesh$nodes()[,2],
                          coeff=x$coeff[1:nrow(x$FunctionalSpace$mesh$nodes())])
  I=x$FunctionalSpace$mesh$elements()[,1]
  J=x$FunctionalSpace$mesh$elements()[,2] 
  K=x$FunctionalSpace$mesh$elements()[,3]
  fig <- plot_ly(plot_data, type="contour", x=~X, y=~Y, z=~coeff, 
                 i = I, j = J, k = K,
                 intensity=~coeff, color = ~coeff,
                 contours=list(showlabels = TRUE),
                 colorbar=list(title=""), ...) %>%
    layout(xaxis = list(title = ""),
           yaxis = list(title = ""))
  fig
}
)

# gradient of Function
.FunctionGradCtr <- setRefClass(
  Class = "FunctionGradObject",
  fields = c(
    f = "FunctionObject", ## Function of which the gradient is taken
    K = "ANY" ## set by product operator
  )
)
