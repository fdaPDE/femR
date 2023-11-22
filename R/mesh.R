
## Mesh Class Definition
.MeshCtr <- setRefClass(
  Class = "MeshObject", contains = "DomainObject",
  fields = c(
    data  = "ANY",
    m = "integer",   # local dim
    n = "integer",   # embedding dim
    times = "numeric",
    deltaT = "numeric" # for parabolic problems
  ),
  methods = c(
    get_nodes = function(){
      data$nodes()
    },
    get_elements = function(){
      data$elements()
    },
    get_boundary = function(){
      data$boundary()
    },
    get_neighbors = function(){
      data$neighbors()
    },
    get_times = function(){
      return(times)
    },
    set_deltaT = function(deltaT){
      if(length(time_interval)==0)
        stop("Error!")
      deltaT <<- deltaT
      times <<- seq(time_interval[1], time_interval[2], by=deltaT)
    }
  )
)

#' Create mesh object
#'
#' @param domain could be a \code{triangulation} returned by \code{\link[RTriangle]{triangulate}} or a named list containing:
#' \itemize{
#'    \item{\code{nodes}, a #nodes-by-2 matrix containing the x and y coordinates of the mesh nodes;}
#'    \item{\code{elements}, a #elements-by-3 matrix specifiying the triangles giving the row's indices in \code{nodes} of the triangles' vertices;}
#'    \item{\code{boundary}, a #nodes-by-1 matrix, with entries either '1' or '0'. An entry '1' indicates that the corresponding node is a boundary node; 
#'           an entry '0' indicates that the corresponding node is not a boundary node.}
#' }
#' 
#' @return An S4 object representing a Mesh.
#' @rdname MeshObject
#' @export
#' @examples
#' \dontrun{
#' library(RTriangle)
#' library(femR)
#' p <- pslg(P=rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1)),
#' S=rbind(c(1, 2), c(2, 3), c(3, 4), c(4,1)))
#' unit_square <- triangulate(p, a = 0.00125, q=30)
#' mesh <- Mesh(unit_square)
#' }

setGeneric("Mesh", function(domain) standardGeneric("Mesh"))

#' @rdname MeshObject
setMethod("Mesh", signature = c(domain="list"),
          function(domain){
            domain$elements <- domain$elements - 1
            storage.mode(domain$elements) <- "integer"
            m <- ncol(domain$elements) - 1
            n <- ncol(domain$nodes) 
            if(m == 2 & n == 2)
              .MeshCtr(data=new(Mesh_2D, domain), m=as.integer(m),n=as.integer(n), 
                       times=vector(mode="double"))                                                  
            else
              stop("wrong input argument provided.")
})

# @importClassesFrom RTriangle triangulation

#' @rdname MeshObject
setMethod("Mesh", signature=c(domain="triangulation"),
          function(domain){
            elements <- domain$T - 1
            nodes <- domain$P
            boundary <- matrix(0, nrow=nrow(nodes), ncol=1)
            boundary[as.vector(domain$E[domain$EB == 1,]),] = 1
            
            storage.mode(elements) <- "integer"
            storage.mode(nodes) <- "numeric"
            storage.mode(boundary) <- "integer"
            
            m <- ncol(elements) - 1
            n <- ncol(nodes)
            
            domain <- list(elements = elements, nodes = nodes, boundary = boundary)
            if(m == 2 & n == 2)
              .MeshCtr(data=new(Mesh_2D, domain), m=as.integer(m),n=as.integer(n), 
                       times=vector(mode="double"))                                                  
            else
              stop("wrong input argument provided.")
            
})

# create spatio-temporal domain
#
# @param op1 A mesh object created by \code{Mesh}.
# @param op2 A numeric vector.
# @return An S4 object representing a spatio-temporal domain.
# @rdname DomainObject_times_vector
# @export 
#setGeneric("%X%", function(op1, op2) standardGeneric("%X%"))

# @rdname MeshObject_times_vector
# setMethod("%X%", signature=c(op1="MeshObject", op2="numeric"),
#           function(op1, op2){
#             if(op2[1] > op2[length(op2)])
#               stop("Error! First time instant is greater than last time instant.")
#             op1$times <- times
#             op1          
# })

## Mesh - auxiliary methods
unroll_edges_aux <- function(Mesh){
  mesh <- Mesh$data
  edges <- matrix(nrow=3*nrow(mesh$elements()), ncol=2)
  for(i in 1:nrow(mesh$elements())){
    edges[(3*(i-1) + 1),]   = mesh$elements()[i,c(1,2)] + 1
    edges[(3*(i-1) + 2),]   = mesh$elements()[i,c(2,3)] + 1
    edges[(3*(i-1) + 3),]   = mesh$elements()[i,c(3,1)] + 1
  }
  edges
}

setGeneric("unroll_edges", function(Mesh) standardGeneric("unroll_edges"))
setMethod("unroll_edges", "MeshObject", function(Mesh){
  unroll_edges_aux(Mesh)
})

plot_mesh_aux <- function(x, ...){
  edges <- unroll_edges(x)
  plot_ly(...) %>% 
    add_markers(x = x$data$nodes()[,1],
                y = x$data$nodes()[,2],
                color = I('black'), size = I(1),
                hoverinfo = 'text',
                text = paste('</br><b> Coordinates:', round(x$data$nodes()[,1],2),
                             round(x$data$nodes()[,2],2)),
                showlegend = T,
                visible = T) %>%
    add_segments(x = x$data$nodes()[edges[,1],1],
                 y = x$data$nodes()[edges[,1],2],
                 xend = x$data$nodes()[edges[,2],1],
                 yend = x$data$nodes()[edges[,2],2], 
                 color = I('black'), size = I(1),
                 showlegend = F) %>%
    layout(
      xaxis = list(
        title = '',
        showgrid = F,
        zeroline = F,
        showticklabels = F
      ),
      yaxis = list(
        title = '',
        showgrid = F,
        zeroline = F,
        showticklabels = F
      ))
}


# setMethod("plot", signature=c(x="MeshObject"), function(x, ...){
#   plot_mesh_aux(x, ...)  
# })

#' Plot a Mesh object
#'
#' @param x A \code{MeshObject} object defining the triangular mesh, as generated by \code{Mesh}
#' @param ... Arguments representing graphical options to be passed to \code{\link[plotly]{plot_ly}}.
#' @return A plotly object
#' @export
#' @examples
#' \dontrun{
#' library(femR)
#' data("unit_square")
#' mesh <- Mesh(unit_square)
#' plot(mesh)
#' }
plot.MeshObject <-function(x, ...){
  plot_mesh_aux(x, ...)
}