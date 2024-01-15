#' sf methods for Domain and Mesh
#'
#' \code{\link[sf]{sf}} methods for \code{Domain} and \code{Mesh} objects.
##
#' @param x An object of class \code{Domain} or \code{Mesh}.
#' @param ... Arguments passed on the corresponding \code{sf} function.
#' @param value The value to be assigned. See the documentation of the
#' corresponding sf function for details.
#' @return The \code{femR} method for \code{\link[sf]{st_as_sf}} returns
#' ............ \code{\link[sf]{sf}}.
#' @importFrom sf st_as_sf
#' @importFrom sf st_polygon
#' @importFrom sf st_as_sfc
#' @rdname sf
#' @export
st_as_sf.Domain <- function(x, ...){
  geom_sf <- list()
  for(sub_id in 1:length(extract_private(x)$geometry_)){
    
    geometry <- extract_private(x)$geometry_[[sub_id]]
    path_id <- unique(geometry$edges_ring)
    
    path_list = vector(mode="list", length=length(path_id))
    for(i in 1:length(path_id)){
      path  <- t(geometry$edges[geometry$edges_ring==path_id[i],][,1]) 
      path  <- c(path,path[1]) 
      nodes <- as.matrix(extract_private(x)$coords_[path,1:2]) 
      path_list[[i]] <- st_linestring(nodes)
    }
    
    # edges to st_linestring ! (FAST)
    edge_list <- lapply(as.list(as.data.frame(t(geometry$edges))), FUN=
                          function(edge){
                            st_linestring(as.matrix(extract_private(x)$coords_[edge,1:2]))
                          })
    
    # nodes to st_point ! (FAST)
    pts_list <- lapply(as.list(as.data.frame(t(geometry$nodes))), st_point)
    
    label= c(rep("path", times=length(path_id)),
             rep("edge", times=length(edge_list)),
             rep("node", times=length(pts_list)))
    group = c(path_id , geometry$edges_group, geometry$nodes_group)
    id =  c(1:length(path_id),
            1:nrow(geometry$edges), 
            1:nrow(geometry$nodes))
    
    group <- as.factor(group)
    crs <- extract_private(x)$crs_
    if(is.na(crs)) crs <- NA_crs_
    geom_sfc <- st_sfc( append(append(path_list, edge_list), pts_list), crs=crs)
    geom_sf[[sub_id]] <- st_as_sf(data.frame(label= label,
                                             group = group, id=id),
                                  geometry = geom_sfc)
    
  }
  if(length(extract_private(x)$geometry_)==1) geom_sf <- geom_sf[[1]]
  return(geom_sf)
}

#' @importFrom sf st_as_sf
#' @importFrom sf st_polygon
#' @importFrom sf st_linestring
#' @importFrom sf st_point
#' @importFrom sf st_sfc
#' @importFrom sf st_cast
#' @rdname sf
#' @export
st_as_sfc.Domain <- function(x, ...){
  polygon_list <- list(mode="list", length=length(extract_private(x)$geometry_))
  for(sub_id in 1:length(extract_private(x)$geometry_)){
    geometry <- extract_private(x)$geometry_[[sub_id]]
    
    path_id <- unique(geometry$edges_ring)
    path_list = vector(mode="list", length=length(path_id))
    for(i in 1:length(path_id)){
      path  <- t(geometry$edges[geometry$edges_ring==path_id[i],])[1,]
      path  <- c(path,path[1]) 
      nodes <- as.matrix(extract_private(x)$coords_[path,1:2]) 
      path_list[[i]] <- nodes
    }
    polygon_list[[sub_id]] <- st_cast(st_polygon(path_list),
                                      to="MULTIPOLYGON")
  }
  crs <- extract_private(x)$crs_
  if(is.na(crs)) crs <- NA_crs_
  
  if(length(extract_private(x)$geometry_)==1){
    result <- st_sfc(polygon_list[[1]], crs=crs)
  }else{
    result <- st_sfc(polygon_list, crs=crs)
  }
  return(result)
}

#' @importFrom sf st_as_sf
#' @importFrom sf st_polygon
#' @importFrom sf st_sfc
#' @importFrom sf st_cast
#' @rdname sf
#' @export
st_as_sfc.Mesh <- function(x, ...){

  polygon_list <- apply(x$elements(), MARGIN=1, FUN=function(elem){
    st_cast(st_linestring(x$nodes()[elem,]),
            to ="POLYGON"
    )
  })
  
  crs <- extract_private(x)$crs_
  if(is.na(crs)) crs <- NA_crs_
  mesh_sf <- st_sfc(polygon_list, crs = crs)
  mesh_sf
}

#' @importFrom sf st_crs
#' @rdname sf
#' @export
st_crs.Domain <- function(x, ...){
  st_crs(st_as_sfc(x), ...)
}

#' @importFrom sf st_crs<- st_crs
#' @rdname sf
#' @export
`st_crs<-.Domain` = function(x, value) {
  sfc <- st_as_sfc(x)
  st_crs(sfc) <- value
  set_private(x, "crs_", value)
}

