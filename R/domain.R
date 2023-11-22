.DomainCtr <- setRefClass("DomainObject", 
                          fields=c(
                            geometry = "ANY",
                            time_interval = "numeric",
                            crs ="ANY"
                          ),
                          methods=c(
                            set_crs = function(crs){
                              crs <<- crs
                            }
                          )
                          
)

#' S4 class representing a spatial (spatio-temporal) domain
#' 
#' @param x could be a \code{pslg} returned by \code{\link[RTriangle]{triangulate}}, 
#' a \code{sfc} returned by \code{\link[sf]{st_geometry}} or a named list containing:
#' \itemize{
#'    \item{\code{nodes}, a #nodes-by-2 matrix containing the x and y coordinates of the boundary nodes;}
#' }
#' @rdname DomainObject
#' @export 
setGeneric("Domain", function(x) standardGeneric("Domain"))

#' @rdname DomainObject
setMethod("Domain", signature = "list", function(x){
  geometry <- list(nodes=x$nodes, nodes_group=rep(1, times=nrow(x$nodes)),
                   edges=x$edges, edges_group=rep(1, times=nrow(x$edges))) 
  
  crs = NA
  if(!is.null(x$crs)) crs <- x$crs
  .DomainCtr(geometry = geometry, time_interval = vector(mode="numeric", length = 0), crs =crs)
})

#' @rdname DomainObject
setMethod("Domain", signature = "pslg", function(x){
  nodes <- x$P
  nodes_group <- as.integer(x$PB)
  edges <- x$S
  edges_group <- as.integer(x$SB)
  geometry <- list(nodes = nodes, edges = edges,
                   nodes_group = nodes_group, edges_group = edges_group)
  .DomainCtr(geometry = geometry, time_interval = vector(mode="numeric", length = 0), crs = NA)
})

#' @rdname DomainObject
setMethod("Domain", signature = "sfc", function(x){
  if(!any(class(x) %in% c("sfc_POLYGON", "sfc_MULTIPOLYGON")))
    stop("Error!")
  x_boundary <- x %>% st_union()             # return a sfc_POLYGON with 1 feature
  nodes <- unique(st_coordinates(x_boundary))[,1:2]
  edges <- cbind(1:nrow(nodes),c(2:nrow(nodes),1))
  nodes_group <- rep(1, times=nrow(nodes))
  edges_group <- rep(1, times=nrow(edges))
  geometry <- list(nodes = nodes, edges = edges,
                   nodes_group = nodes_group, edges_group = edges_group)
  crs <- NA
  if(!is.na(st_crs(x))) crs <- st_crs(x)$input
  .DomainCtr(geometry = geometry, time_interval = vector(mode="numeric", length = 0), crs = crs)
})

#' Delaunay triangulation of the spatial domain
#' 
#' @param Domain 
#' @param maximum_area maximum triangle area
#' @param minumum_angle minimum triangle angle in degrees
#' @rdname Triangulate
#' @export
setGeneric("Triangulate", function(Domain, maximum_area, minimum_angle) standardGeneric("Triangulate"))

#' @rdname Triangulate
setMethod("Triangulate", signature=c("DomainObject","numeric", "numeric"), 
          function(Domain, maximum_area, minimum_angle){
            pslg <- RTriangle::pslg(P = Domain$geometry$nodes, PB = Domain$geometry$nodes_group,
                                    S = Domain$geometry$edges, SB = Domain$geometry$edges_group) 
            triangulation <- RTriangle::triangulate(p = pslg, a = maximum_area, q = minimum_angle)
            res <- Mesh(triangulation)
            res$geometry <- Domain$geometry
            res$time_interval <- Domain$time_interval
            res$set_crs(Domain$crs)
            return(res)
})

#' create spatio-temporal domain
#'
#' @param op1 A mesh object created by \code{Mesh}.
#' @param op2 A numeric vector.
#' @return An S4 object representing a spatio-temporal domain.
#' @rdname DomainObject_times_vector
#' @export 
setGeneric("%X%", function(op1, op2) standardGeneric("%X%"))

#' @rdname DomainObject_times_vector
setMethod("%X%", signature=c(op1="DomainObject", op2="numeric"),
          function(op1, op2){
            if(op2[1] > op2[length(op2)])
              stop("Error! First time instant is greater than last time instant.")
            op1$time_interval <- op2
            op1          
})