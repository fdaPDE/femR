.MeshCtr <- R6Class("Mesh", inherit = .DomainCtr,
                    private = list(
                      cpp_handler_ = "ANY",                        ## cpp backend
                      m_ = 0L,                                     ## local dim
                      n_ = 0L,                                     ## embedding dim
                      time_nodes_ = vector(mode="double", length = 0L), ## time mesh for parabolic problems
                      time_step_ = NA                              ## time step for parabolic problems
                    ),
                    public = list(
                      initialize = function(cpp_handler, m, n){
                        private$cpp_handler_ <- cpp_handler
                        private$m_ <- m
                        private$n_ <- n
                      },
                      nodes = function(){
                        private$cpp_handler_$nodes()
                      },
                      elements = function(){
                        private$cpp_handler_$elements()+1
                      },
                      boundary = function(){
                        private$cpp_handler_$boundary()
                      },
                      neighbors = function(){
                        private$cpp_handler_$neighbors()
                      },
                      time_nodes = function(){
                        private$time_nodes_
                      },
                      set_time_step = function(time_step){
                        if(length(private$time_interval_)==0)
                          stop("Error!")
                        private$time_step_ <- time_step
                        private$time_nodes_ <- seq(private$time_interval_[1], private$time_interval_[2], by=time_step)
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
#' @return A R6 class representing a Mesh.
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
# @return A R6 object representing a spatio-temporal domain.
# @rdname Domain_times_vector
# @export
#setGeneric("%X%", function(op1, op2) standardGeneric("%X%"))

#' @rdname Domain_times_vector
setMethod("%X%", signature=c(op1="Mesh", op2="numeric"),
          function(op1, op2){
            if(op2[1] > op2[length(op2)])
              stop("Error! First time instant is greater than last time instant.")
            set_private(op1, "time_nodes_", op2)
            set_private(op1, "time_interval_", c(op2[1], op2[length(op2)]))
            set_private(op1, "time_step_", (op2[2] - op2[1]))
            op1
})

## Mesh - auxiliary methods
# unroll_edges_aux <- function(mesh){
#   edges <- matrix(nrow=3*nrow(mesh$elements()), ncol=2)
#   for(i in 1:nrow(mesh$elements())){
#     edges[(3*(i-1) + 1),]   = mesh$elements()[i,c(1,2)] 
#     edges[(3*(i-1) + 2),]   = mesh$elements()[i,c(2,3)] 
#     edges[(3*(i-1) + 3),]   = mesh$elements()[i,c(3,1)] 
#   }
#   edges
# }

setGeneric("unroll_edges", function(mesh) standardGeneric("unroll_edges"))
setMethod("unroll_edges", "Mesh", function(mesh){
  #unroll_edges_aux(mesh)
  edges <- matrix(nrow=3*nrow(mesh$elements()), ncol=2)
  for(i in 1:nrow(mesh$elements())){
    edges[(3*(i-1) + 1),]   = mesh$elements()[i,c(1,2)] 
    edges[(3*(i-1) + 2),]   = mesh$elements()[i,c(2,3)] 
    edges[(3*(i-1) + 3),]   = mesh$elements()[i,c(3,1)] 
  }
  edges
})

plot_mesh_aux <- function(x, ...){
  edges <- unroll_edges(x)
  plot_ly(...) %>% 
    add_markers(x = x$nodes()[,1],
                y = x$nodes()[,2],
                color = I('black'), size = I(1),
                hoverinfo = 'text',
                text = paste('</br><b> Coordinates:', round(x$nodes()[,1],2),
                             round(x$nodes()[,2],2)),
                showlegend = T,
                visible = T) %>%
    add_segments(x = x$nodes()[edges[,1],1],
                 y = x$nodes()[edges[,1],2],
                 xend = x$nodes()[edges[,2],1],
                 yend = x$nodes()[edges[,2],2], 
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