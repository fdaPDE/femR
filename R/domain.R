.DomainCtr <- setRefClass("DomainObject", 
                          fields=c(
                            geometry = "ANY",
                            coords = "data.frame",         # (x,y,boundary)
                            time_interval = "numeric",
                            crs ="ANY"
                           )#,
                          # methods=c(
                          #    set_crs = function(crs){
                          #      crs <<- crs
                          #    }
                          # )
                          
)

#' S4 class representing a spatial (spatio-temporal) domain
#' 
#' @param x could be a \code{pslg} returned by \code{\link[RTriangle]{triangulate}}, 
#' a \code{sfc} returned by \code{\link[sf]{st_geometry}} or a named list containing:
#' \itemize{
#'    \item{\code{nodes}, a #nodes-by-2 matrix containing the x and y coordinates of the boundary nodes;}
#' }
#' @rdname DomainObject
#' @importFrom sf NA_crs_
#' @importFrom sf st_crs
#' @export 
setGeneric("Domain", function(x) standardGeneric("Domain"))

#' @rdname DomainObject
setMethod("Domain", signature = "list", function(x){
    if(is.null(x$nodes_group)) x$nodes_group <- rep(1, times=nrow(x$nodes))
    if(is.null(x$edges_group)) x$edges_group <- rep(1, times=nrow(x$edges))
    if(is.null(x$nodes_boundary)) x$nodes_boundary <- rep(1, times=nrow(x$nodes))
    if(is.null(x$edges_boundary)) x$edges_boundary <- rep(1, times=nrow(x$edges))
    
    geometry <- vector(mode="list", length=1L)
    geometry[[1]] <- list(nodes=x$nodes, nodes_group=x$nodes_group, nodes_boundary=x$nodes_boundary,
                          edges=x$edges, edges_group=x$edges_group, edges_boundary=x$edges_boundary) 
    
    coords <- data.frame(x=x$nodes[,1], y=x$nodes[,2], boundary=x$nodes_boundary)
  crs = NA_crs_
  if(!is.null(x$crs)) crs <- st_crs(x$crs)
  .DomainCtr(geometry = geometry, time_interval = vector(mode="numeric", length = 0), 
             coords=coords, crs = crs)
})

#' @rdname DomainObject
setMethod("Domain", signature = "pslg", function(x){
  nodes <- x$P
  nodes_group <- as.integer(x$PB)
  edges <- x$S
  edges_group <- as.integer(x$SB)
  geometry <- list(nodes = nodes, edges = edges,
                   nodes_group = nodes_group, edges_group = edges_group)
  .DomainCtr(geometry = geometry, time_interval = vector(mode="numeric", length = 0),
             coords=geometry$nodes, crs = NA_crs_)
})

#' @rdname DomainObject
setMethod("Domain", signature="sfc", function(x){
  if(length(x) == 1) x <- st_cast(x, to="MULTIPOLYGON") 
  coords <- st_coordinates(x)
  #coords <- coords[-nrow(coords),]
  coords <- cbind(coords, 
                  rep(0, times=nrow(coords)), # id 
                  rep(0, times=nrow(coords))) # boundary (1,0)
  coords <- as.data.frame(coords)
  colnames(coords)[(ncol(coords)-1):ncol(coords)] <- c("id", "boundary")
  storage.mode(coords$id) <- "integer"; storage.mode(coords$boundary) <- "integer" 
  pts_list <- vector(mode="list", length=nrow(coords))
  for(i in 1:nrow(coords)){
    pts_list[[i]] <- st_point(as.numeric(coords[i,1:2]))
  }
  pts_sfc <- st_sfc(pts_list)
  
  # setto indici dei nodi al bordo
  coords_bd <- st_coordinates(st_cast(st_union(st_geometry(x)),
                                      to="LINESTRING"))
  pts_bd_list <- vector(mode="list", length=nrow(coords_bd))
  for(i in 1:nrow(coords_bd)){
    pts_bd_list[[i]] <- st_point(coords_bd[i,1:2])
  }
  pts_bd_sfc <- st_sfc(pts_bd_list)
  dist_bd <- st_distance(pts_sfc, pts_bd_sfc)
  
  num_nodes <- 1
  for(i in 1:nrow(coords_bd)){
    id_bd <- which(dist_bd[,i] < 100 *.Machine$double.eps )
    if(any( coords[id_bd,(ncol(coords)-1)] != 0)){
      global_id <- which(coords[id_bd,(ncol(coords)-1)] != 0)[1]
      coords[id_bd,(ncol(coords)-1)] <-  rep(global_id, times=length(id_bd))
    }else{
      coords[id_bd,(ncol(coords)-1)] <- rep(num_nodes, times=length(id_bd))
      num_nodes <- num_nodes + 1  
    }
    coords[id_bd, ncol(coords)] <- rep(1, times=length(id_bd))
  }
  
  for(sub_id in 1:length(st_geometry(x))){
    coords_sub_id <- st_coordinates(st_geometry(x)[sub_id])
    pts_sub_id_list <- vector(mode="list", length=nrow(coords_sub_id))
    for(j in 1:nrow(coords_sub_id)){
      pts_sub_id_list[[j]] <- st_point(coords_sub_id[j,1:2])
    }
    pts_sub_id_sfc <- st_sfc(pts_sub_id_list)
    dist_matrix_sub_id <- st_distance(pts_sfc, pts_sub_id_sfc)
    for(j in 1:nrow(coords_sub_id)){
      id <- which(dist_matrix_sub_id[,j] < 100* .Machine$double.eps)
      if(all(coords[id,(ncol(coords)-1)] != 0 )){ next }
      if( any( coords[id,(ncol(coords)-1)] != 0) ){
        global_id <- which(coords[id,(ncol(coords)-1)] != 0)[1]
        coords[id,(ncol(coords)-1)] <-  rep(global_id, times=length(id))
      }else{
        coords[id,(ncol(coords)-1)] <- rep(num_nodes, times = length(id))
        num_nodes = num_nodes + 1
      }
    }
  }
  
  geometry <- vector(mode="list", length =length(st_geometry(x)))
  
  for(sub_id in 1:length(st_geometry(x))){
    idx <- which(coords[,5]==sub_id)  
    nodes <- unique(coords[idx,])
    n_main_ring <- diff(range(nodes[,4])) + 1 # number of main ring!!! see st_coordinates
    edge <- matrix(nrow=0, ncol=2)
    edge_sub <- matrix(nrow=0, ncol=2)
    edge_group <- matrix(nrow=0, ncol=1)
    edge_bd <- matrix(nrow=0, ncol=1)
    node_group <- matrix(nrow=0, ncol=1)
    for(ring in 1:n_main_ring){
      idx_ring <- which(nodes[,4] == ring)
      nodes_ring <- nodes[idx_ring,]
      # in the current main ring
      n_physical_lines <- diff(range(nodes_ring[,3])) + 1 
      for(i in 1:n_physical_lines){
        id_edge_sub <- which(nodes_ring[,3] == i)
        # local
        edge_sub_i <-  cbind(id_edge_sub,
                             c(id_edge_sub[2:length(id_edge_sub)],id_edge_sub[1]))
        # global 
        edge_i <- cbind(nodes_ring[edge_sub_i[,1],(ncol(coords)-1)], nodes_ring[edge_sub_i[,2],(ncol(coords)-1)])           
        
        edge_bd_i <- matrix(0, nrow=nrow(edge_sub_i),ncol=1)
        node_group_i <-matrix(0, nrow(nodes[id_edge_sub,]),ncol=1)
        for(e in 1:nrow(edge_sub_i)){
          if( nodes_ring[edge_sub_i[e,1],ncol(coords)] & nodes_ring[edge_sub_i[e,2],ncol(coords)] ){ # sono entrambi sul bordo
            edge_bd_i[e] <- 1
            node_group_i[edge_sub_i[e,1]] <- 1
            node_group_i[edge_sub_i[e,2]] <- 1
          } 
        }
        if(i == 1){
          edge_group_i <- as.matrix(rep(ring, times=nrow(edge_sub_i)))
        }else{
          edge_group_i <- as.matrix(rep((-i+1), times=nrow(edge_sub_i)))
        }
        edge <- rbind(edge, edge_i)
        edge_group <- rbind(edge_group, edge_group_i)
        edge_bd <- rbind(edge_bd, edge_bd_i)
        edge_sub <- rbind(edge_sub, edge_sub_i)
        node_group <- rbind(node_group, node_group_i)
      }
    }
    storage.mode(edge) <-"integer"
    storage.mode(edge_bd) <-"integer"
    storage.mode(nodes[,ncol(coords)]) <- "integer" 
    geometry[[sub_id]] <- list(nodes = nodes[,1:2], nodes_group = node_group, nodes_boundary = nodes[,ncol(coords)],
                               edges = edge, edges_group = edge_group, edges_boundary = edge_bd)
    
  }
  # da salvare
  NODES <- unique(coords[order(coords[,(ncol(coords)-1)]),][,c(1:2,ncol(coords))])
  NODES <- as.data.frame(NODES)
  storage.mode(NODES$boundary) <- "integer"
  
  crs <- NA_crs_
  if(!is.na(st_crs(x))) crs <- st_crs(x)
      .DomainCtr(geometry = geometry, time_interval = vector(mode="numeric", length = 0), 
                coords=NODES, crs = crs)
})

# setMethod("Domain", signature = "sfc", function(x){
#   if(!any(class(x) %in% c("sfc_POLYGON", "sfc_MULTIPOLYGON")))
#     stop("Error!")
#   x_boundary <- x %>% st_union()             # return a sfc_POLYGON with 1 feature
#   nodes <- unique(st_coordinates(x_boundary))[,1:2]
#   edges <- cbind(1:nrow(nodes),c(2:nrow(nodes),1))
#   nodes_group <- rep(1, times=nrow(nodes))
#   edges_group <- rep(1, times=nrow(edges))
#   geometry <- list(coords = nodes, edges = edges,
#                    nodes_group = nodes_group, edges_group = edges_group)
#   crs <- NA_crs_
#   if(!is.na(st_crs(x))) crs <- st_crs(x)
#   .DomainCtr(geometry = geometry, time_interval = vector(mode="numeric", length = 0), 
#              coords=geometry$nodes, crs = crs)
# })

#' Delaunay triangulation of the spatial domain
#' 
#' @param Domain 
#' @param maximum_area maximum triangle area
#' @param minumum_angle minimum triangle angle in degrees
#' @rdname build_mesh
#' @importFrom RTriangle pslg
#' @importFrom RTriangle triangulate
#' @importFrom sf NA_crs_
#' @export
setGeneric("build_mesh", function(domain, maximum_area, minimum_angle) standardGeneric("build_mesh"))

#' @rdname build_mesh
setMethod("build_mesh", signature=c("DomainObject", "numeric", "numeric"),
          function(domain, maximum_area, minimum_angle){
  S = matrix(nrow=0, ncol=2)
  SB = matrix(nrow=0, ncol=1)
  H = matrix(nrow=0, ncol=2)
  
  for( sub_id in 1:length(domain$geometry)){
    groups_id <- unique(domain$geometry[[sub_id]]$edges_group)
    holes_id <- groups_id[ groups_id < 0 ]
    num_holes <- length(holes_id)
    #holes <- matrix(0, nrow=num_holes, ncol=2)
    if(num_holes > 0){
      for(i in 1:num_holes){
        edges_id <- which(domain$geometry[[sub_id]]$edges_group==holes_id[i])
        if( ! all(as.integer(domain$geometry[[sub_id]]$edges_boundary[edges_id])) ){ # NON VERO buco
          next
        } 
        path  <- t(domain$geometry[[sub_id]]$edges[edges_id,])[1,]
        path  <- c(path,path[1]) 
        nodes <- as.matrix(domain$coords[path,1:2])
        holes <- st_point_on_surface(st_polygon(list(nodes)))
        H = rbind(H, holes)
      }
    }
    S = rbind(S, domain$geometry[[sub_id]]$edges)
    SB = rbind(SB, as.matrix(domain$geometry[[sub_id]]$edges_boundary))
  }
  
  if(nrow(H)==0){
    H = NA
  }
  pslg <- pslg(P = as.matrix(domain$coords[,1:2]), PB = as.matrix(domain$coords[,3]),
               S = S, SB = SB, H = H) 
  triangulation <- triangulate(p = pslg, a = maximum_area, q = minimum_angle)
  res <- Mesh(triangulation)
  res$geometry <- domain$geometry
  res$time_interval <- domain$time_interval
  res$crs <- domain$crs
  return(res)
})
# setMethod("build_mesh", signature=c("DomainObject","numeric", "numeric"), 
#           function(domain, maximum_area, minimum_angle){
#             groups_id <- unique(domain$geometry$edges_group)
#             holes_id <- groups_id[ groups_id < 0 ]
#             num_holes <- length(holes_id)
#             holes <- matrix(0, nrow=num_holes, ncol=2)
#             if(num_holes > 0){
#               for(i in 1:num_holes){
#                 path  <- t(domain$geometry$edges[domain$geometry$edges_group==holes_id[i],])[1,]
#                 path  <- c(path,path[1]) 
#                 nodes <- domain$geometry$nodes[path,]
#                 holes[i,] <- st_point_on_surface(st_polygon(list(nodes)))
#               }
#             }else{
#               holes <- NA
#             }
#             pslg <- pslg(P = domain$geometry$nodes, PB = domain$geometry$nodes_group,
#                          S = domain$geometry$edges, SB = domain$geometry$edges_group,
#                          H = holes) 
#             triangulation <- triangulate(p = pslg, a = maximum_area, q = minimum_angle)
#             res <- Mesh(triangulation)
#             res$geometry <- domain$geometry
#             res$time_interval <- domain$time_interval
#             res$crs <- domain$crs
#             return(res)
# })

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