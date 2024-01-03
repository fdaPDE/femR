
## Mesh Class Definition
#' @export 
.MeshCtr <- R6Class("Mesh", inherit = .DomainCtr,
                    public = list(
                      cpp_handler = "ANY",       ## cpp backend
                      m = 0L,             ## local dim
                      n = 0L,             ## embedding dim
                      times = vector(mode="double", length = 0L),
                      deltaT = NA, # for parabolic problems
                      initialize = function(cpp_handler, m, n){
                        self$cpp_handler <- cpp_handler
                        self$m <- m
                        self$n <- n
                      },
                      get_nodes = function(){
                        self$cpp_handler$nodes()
                      },
                      get_elements = function(){
                        self$cpp_handler$elements()+1
                      },
                      get_boundary = function(){
                        self$cpp_handler$boundary()
                      },
                      get_neighbors = function(){
                        self$cpp_handler$neighbors()
                      },
                      get_times = function(){
                        return(self$times)
                      },
                      set_deltaT = function(deltaT){
                        if(length(self$time_interval)==0)
                          stop("Error!")
                        self$deltaT <- deltaT
                        self$times <- seq(self$time_interval[1], self$time_interval[2], by=deltaT)
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
#' @rdname Mesh
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

#' @rdname Mesh
setMethod("Mesh", signature = c(domain="list"),
          function(domain){
            domain$elements <- domain$elements - 1
            storage.mode(domain$elements) <- "integer"
            m <- ncol(domain$elements) - 1
            n <- ncol(domain$nodes) 
            if(m == 2 & n == 2)
              .MeshCtr$new(cpp_handler=new(cpp_2d_domain, domain), m=as.integer(m),n=as.integer(n))                                                  
            else
              stop("wrong input argument provided.")
})

setOldClass("triangulation")

#' @rdname Mesh
setMethod("Mesh", signature=c(domain="triangulation"),
          function(domain){
            nodes <- domain$P
            elements <- domain$T - 1
            boundary <- domain$PB
            
            storage.mode(elements) <- "integer"
            storage.mode(nodes) <- "numeric"
            storage.mode(boundary) <- "integer"
            
            m <- ncol(elements) - 1
            n <- ncol(nodes)
            
            domain <- list(elements = elements, nodes = nodes, boundary = boundary)
            if(m == 2 & n == 2)
              .MeshCtr$new(cpp_handler=new(cpp_2d_domain, domain), m=as.integer(m),n=as.integer(n))                                                  
            else
              stop("wrong input argument provided.")
})

# create spatio-temporal mesh
# 
# @param op1 A mesh object created by \code{Mesh}.
# @param op2 A numeric vector.
# @return An S4 object representing a spatio-temporal domain.
# @rdname Domain_times_vector
# @export
#setGeneric("%X%", function(op1, op2) standardGeneric("%X%"))

#' @rdname Domain_times_vector
setMethod("%X%", signature=c(op1="Mesh", op2="numeric"),
          function(op1, op2){
            if(op2[1] > op2[length(op2)])
              stop("Error! First time instant is greater than last time instant.")
            op1$times <- times
            op1$time_interval <- c(times[1], times[length(times)])
            op1$deltaT <- times[2] - times[1]
            op1
})

## Mesh - auxiliary methods
unroll_edges_aux <- function(Mesh){
  mesh <- Mesh$cpp_handler
  edges <- matrix(nrow=3*nrow(mesh$elements()), ncol=2)
  for(i in 1:nrow(mesh$elements())){
    edges[(3*(i-1) + 1),]   = mesh$elements()[i,c(1,2)] + 1
    edges[(3*(i-1) + 2),]   = mesh$elements()[i,c(2,3)] + 1
    edges[(3*(i-1) + 3),]   = mesh$elements()[i,c(3,1)] + 1
  }
  edges
}

setGeneric("unroll_edges", function(Mesh) standardGeneric("unroll_edges"))
setMethod("unroll_edges", "Mesh", function(Mesh){
  unroll_edges_aux(Mesh)
})

plot_mesh_aux <- function(x, ...){
  edges <- unroll_edges(x)
  plot_ly(...) %>% 
    add_markers(x = x$cpp_handler$nodes()[,1],
                y = x$cpp_handler$nodes()[,2],
                color = I('black'), size = I(1),
                hoverinfo = 'text',
                text = paste('</br><b> Coordinates:', round(x$cpp_handler$nodes()[,1],2),
                             round(x$cpp_handler$nodes()[,2],2)),
                showlegend = T,
                visible = T) %>%
    add_segments(x = x$cpp_handler$nodes()[edges[,1],1],
                 y = x$cpp_handler$nodes()[edges[,1],2],
                 xend = x$cpp_handler$nodes()[edges[,2],1],
                 yend = x$cpp_handler$nodes()[edges[,2],2], 
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


# setMethod("plot", signature=c(x="Mesh"), function(x, ...){
#   plot_mesh_aux(x, ...)  
# })

#' Plot a Mesh object
#'
#' @param x A \code{Mesh} object defining the triangular mesh, as generated by \code{Mesh}
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
plot.Mesh <-function(x, ...){
  plot_mesh_aux(x, ...)
}