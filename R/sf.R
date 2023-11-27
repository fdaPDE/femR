#' sf methods for DomainObject and MeshObject
#'
#' \code{\link[sf]{sf}} methods for \code{\link{DomainObject}} and \code{\link{MeshObject}} objects.
##
#' @param x An object of class \code{\link{DomainObject}} or \code{\link{MeshObject}}.
#' @param ... Arguments passed on the corresponding \code{sf} function.
#' @param value The value to be assigned. See the documentation of the
#' corresponding sf function for details.
#' @return The \code{femR} method for \code{\link[sf]{st_as_sf}} returns
#' ............ \code{\link[sf]{sf}}.
#' @importFrom sf st_as_sf
#' @importFrom sf st_polygon
#' @importFrom sf st_as_sfc
#' @export
st_as_sfc.DomainObject <- function(x, ...){
  polygon_list <- list(mode="list", length=length(x$geometry))
  for(sub_id in 1:length(x$geometry)){
    groups_id <- unique(x$geometry[[sub_id]]$edges_group)
    n_physical_lines <- length(groups_id)
    path_list = vector(mode="list", length=n_physical_lines)
    for(i in 1:n_physical_lines){
      path  <- t(x$geometry[[sub_id]]$edges[x$geometry[[sub_id]]$edges_group==groups_id[i],])[1,]
      path  <- c(path,path[1]) 
      nodes <- as.matrix(x$coords[path,1:2]) #x$geometry[[sub_id]]$nodes[path,]
      path_list[[i]] <- nodes
    }
    polygon_list[[sub_id]] <- st_cast(st_polygon(path_list),
                                      to="MULTIPOLYGON")
  }
  if(length(x$geometry)==1){
    result <- st_sfc(polygon_list[[1]], crs=st_crs(x))
  }else{
    result <- st_sfc(polygon_list, crs=st_crs(x))
  }
  return(result)
}


#' @name sf
#' @importFrom sf st_as_sf
#' @importFrom sf st_polygon
#' @importFrom sf st_linestring
#' @importFrom sf st_point
#' @importFrom sf st_sfc
#' @export
st_as_sf.DomainObject <- function(x, ...){
  geom_sf <- list()
  for(sub_id in 1:length(x$geometry)){
    label <- vector(mode="character")
    sub <- vector(mode="integer")
    local_id <- vector(mode="integer")
    group <- vector(mode="numeric")
    boundary <- vector(mode="integer")
    bc = vector(mode="character")
    path_list = vector(mode="list")
    edge_list = vector(mode="list")
    pts_list = vector(mode="list")
    
    geometry <- x$geometry[[sub_id]]
    path_id <- unique(geometry$edges_group)
    n_physical_lines <- length(path_id)
    
    path_list_sub = vector(mode="list", length=n_physical_lines)
    for(i in 1:n_physical_lines){
      path  <- t(geometry$edges[geometry$edges_group==path_id[i],][,1]) # )[1,]
      path  <- c(path,path[1]) 
      nodes <- as.matrix(x$coords[path,1:2]) #geometry$nodes[path,]
      path_list_sub[[i]] <- st_linestring(nodes)
    }
    edge_list_sub = vector(mode="list", length=n_physical_lines)
    for(i in 1:nrow(geometry$edges)){
      edge  <- geometry$edges[i,]
      nodes <- as.matrix(x$coords[edge,1:2]) #geometry$nodes[path,]
      edge_list_sub[[i]] <- st_linestring(nodes)
    }
    
    groups_pts_id <- unique(geometry$nodes_group)
    n_physical_pts <- length(groups_pts_id)
    pts_list_sub = vector(mode="list", length=n_physical_pts)
    for(i in 1:nrow(geometry$nodes)){
      node  <- geometry$nodes[i,]
      pts_list_sub[[i]] <- st_point(t(as.matrix(node)))
    }
    
    label= append(label, c(rep("path", times=n_physical_lines),
                           rep("edge", times=nrow(geometry$edges)),
                           rep("node", times=nrow(geometry$nodes))))
    group = append(group, c(path_id , geometry$edges_group, geometry$nodes_group))
    local_id = append(local_id, 
                      c(1:n_physical_lines,1:nrow(geometry$edges), 1:nrow(geometry$nodes)))
    bc = append(bc, 
                c(rep(NA, times = n_physical_lines), geometry$BC, rep(NA, times=nrow(geometry$nodes))))
    # sub = append(sub, rep(sub_id, 
    #                       times=(n_physical_lines + nrow(geometry$edges) + nrow(geometry$nodes))))  
    path_list = append(path_list, path_list_sub) 
    edge_list = append(edge_list, edge_list_sub)
    pts_list  = append(pts_list, pts_list_sub)
    
    group <- as.factor(group)
    geom_sfc <- st_sfc( append(append(path_list, edge_list), pts_list))
    geom_sf[[sub_id]] <- st_as_sf(data.frame(label= label,
                                             group = group, local_id=local_id,
                                             bc=bc),
                                  geometry = geom_sfc)
    
  }
  if(length(x$geometry)==1) geom_sf <- geom_sf[[1]]
  return(geom_sf)
}
# st_as_sf.DomainObject <- function(x, ...){
#   groups_id <- unique(x$geometry$edges_group)
#   n_physical_lines <- length(groups_id)
#   path_list = vector(mode="list", length=n_physical_lines)
#   for(i in 1:n_physical_lines){
#     path  <- t(x$geometry$edges[x$geometry$edges_group==groups_id[i],])[1,]
#     path  <- c(path,path[1]) 
#     nodes <- x$geometry$nodes[path,]
#     path_list[[i]] <- nodes
#   }
#   st_sfc( st_polygon(path_list) )
# }

#' @importFrom sf st_as_sf
#' @importFrom sf st_polygon
#' @importFrom sf st_sfc
#' @name sf
#' @export
st_as_sfc.MeshObject <- function(x, ...){
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

