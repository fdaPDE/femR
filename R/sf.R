#' sf methods for DomainObject and MeshObject
#'
#' \code{\link[sf]{sf}} methods for \code{\link{DomainObject}} and \code{\link{MeshObject}} objects.
#'
#' @param x An object of class \code{\link{DomainObject}} or \code{\link{MeshObject}}.
#'
#' @param y An object of class of \code{\link[sf:st]{sfg}}.
#' 
#' @param ... Arguments passed on the corresponding \code{sf} function.
#'
#' @param value The value to be assigned. See the documentation of the
#' corresponding sf function for details.
#'
#' @return The \code{femR} method for \code{\link[sf]{st_as_sf}} returns
#' ............ \code{\link[sf]{sf}}.
#' 
#' @details See the \code{\link[sf]{sf}} documentation.
#' @importFrom sf st_as_sf
#' @importFrom sf st_polygon
#' @importFrom sf st_sfc
#' @export
st_as_sf.DomainObject <- function(x, ...){
  groups_id <- unique(x$geometry$edges_group)
  n_physical_lines <- length(groups_id)
  path_list = vector(mode="list", length=n_physical_lines)
  for(i in 1:n_physical_lines){
    path  <- t(x$geometry$edges[x$geometry$edges_group==groups_id[i],])[1,]
    path  <- c(path,path[1]) 
    nodes <- x$geometry$nodes[path,]
    path_list[[i]] <- nodes
  }
  st_sfc( st_polygon(path_list) )
}

#' @importFrom sf st_as_sf
#' @importFrom sf st_polygon
#' @importFrom sf st_sfc
#' @name sf
#' @export
st_as_sf.MeshObject <- function(x, ...){
  polygon_list <- list()
  for(e in 1:nrow(x$get_elements())){
    element_sf <- st_polygon(list(rbind(x$get_nodes()[x$get_elements()[e,],],
                                        x$get_nodes()[x$get_elements()[e,1],])))
    polygon_list[[e]]  <- element_sf
    
  }
  
  mesh_sf <- st_sfc(polygon_list, crs = st_crs(x))
  mesh_sf
}

#' @name sf
#' @importFrom sf st_crs
#' @export
st_crs.DomainObject <- function(x, ...){
  st_crs(x$crs, ...)
}

#' @name sf
#' @importFrom sf st_crs<- st_crs
#' @export
`st_crs<-.DomainObject` = function(x, value) {
  st_crs(x) <- value
}

