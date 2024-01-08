.DomainCtr <- R6Class("Domain", 
                      private = list(
                        geometry = "ANY",
                        coords = "matrix",         # (x,y)
                        boundary = "matrix",        # 
                        time_interval = vector(mode="numeric", length = 0),
                        crs = "ANY"
                      ),
                      public = list(
                        initialize = function(geometry, coords, boundary, crs){
                                      private$geometry <- geometry 
                                      private$coords <- coords 
                                      private$boundary <- boundary
                                      private$crs <- crs
                                    }
                      )
)

#' @name Domain
#'
#' @exportClass Domain
setOldClass(c("Domain", "R6"))

#' Mesh 
#'
#' @name Mesh
#'
#' @exportClass Mesh
setOldClass(c("Mesh","Domain"))

#' S4 class representing a spatial (spatio-temporal) domain
#' 
#' @param x could be a \code{pslg} returned by \code{\link[RTriangle]{triangulate}}, 
#' a \code{sfc} returned by \code{\link[sf]{st_geometry}} or a named list containing:
#' \itemize{
#'    \item{\code{nodes}, a #nodes-by-2 matrix containing the x and y coordinates of the boundary nodes;}
#' }
#' @rdname Domain
#' @importFrom sf NA_crs_
#' @importFrom sf st_crs
#' @importFrom sf st_within
#' @export 
setGeneric("Domain", function(x) standardGeneric("Domain"))

#' @rdname Domain
setMethod("Domain", signature = "list", function(x){
    if(is.null(x$nodes_group)) x$nodes_group <- rep(1, times=nrow(x$nodes))
    if(is.null(x$edges_group)) x$edges_group <- rep(1, times=nrow(x$edges))
    if(is.null(x$nodes_boundary)) x$nodes_boundary <- rep(1, times=nrow(x$nodes))
    if(is.null(x$edges_boundary)) x$edges_boundary <- rep(1, times=nrow(x$edges))
    if(is.null(x$edges_ring))     x$edges_ring <- rep(1, times = nrow(x$edges))
    
    geometry <- vector(mode="list", length=1L)
    geometry[[1]] <- list(nodes=x$nodes, nodes_group=x$nodes_group, nodes_boundary=x$nodes_boundary,
                          edges=x$edges, edges_group=x$edges_group, edges_boundary=x$edges_boundary,
                          edges_ring = x$edges_ring) 
    
    coords <- cbind(x$nodes[,1], x$nodes[,2])
    boundary <- x$nodes_boundary
    crs = NA
  .DomainCtr$new(geometry = geometry, coords=coords, boundary = boundary, crs = crs)
})

setOldClass("pslg")

#' @rdname Domain
setMethod("Domain", signature = "pslg", function(x){
  nodes <- x$P
  edges <- x$S
  nodes_boundary <- x$PB
  edges_boundary <- x$SB
  
  if(!nrow(edges))  edges <- cbind(1:nrow(nodes), c(2:nrow(nodes),1))
  if(!any(nodes_boundary)) nodes_boundary <- rep(1,times=nrow(nodes))
  if(!length(edges_boundary)) edges_boundary <- rep(1, times=nrow(edges))
  
  nodes_group <- rep(1, times=nrow(nodes))
  edges_group <- rep(1, times=nrow(edges))
  edges_ring <- rep(1, times=nrow(edges))
  
  geometry <- vector(mode="list", length=1L)
  geometry[[1]] <- list(nodes=nodes, nodes_group=nodes_group, nodes_boundary=nodes_boundary,
                        edges=edges, edges_group=edges_group, edges_boundary=edges_boundary,
                        edges_ring=edges_ring) 
  
  coords <- cbind(nodes[,1], nodes[,2])
  boundary <- nodes_boundary
  
  .DomainCtr$new(geometry = geometry, coords=coords, boundary=boundary, crs = NA)
})

#' @rdname Domain
#' @importFrom sf st_cast
#' @importFrom sf st_coordinates
#' @importFrom sf st_union
#' @importFrom sf st_geometry
#' @importFrom sf st_distance
#' @importFrom sf st_centroid
setMethod("Domain", signature="sfc", function(x){
  if(length(x) == 1) x <- st_cast(x, to="MULTIPOLYGON") 
  st_coords <- st_coordinates(x)
  
  st_coords <- cbind(st_coords, 
                  rep(0, times=nrow(st_coords)), # id 
                  rep(0, times=nrow(st_coords))) # boundary (1,0)
  st_coords <- as.data.frame(st_coords)
  colnames(st_coords)[(ncol(st_coords)-1):ncol(st_coords)] <- c("id", "boundary")
  storage.mode(st_coords$id) <- "integer"; storage.mode(st_coords$boundary) <- "integer" 
  
  # nodes to st_point ! (FAST)
  pts_list <- lapply(as.list(as.data.frame(t(st_coords[,1:2]))), st_point)
  pts_sfc <- st_sfc(pts_list)
  
  # setto indici dei nodi al bordo
  coords_bd <- st_coordinates(st_cast(st_union(st_geometry(x)),
                                      to="LINESTRING"))
  
  pts_bd_list <- lapply(as.list(as.data.frame(t(coords_bd[,1:2]))), st_point)
  pts_bd_sfc <- st_sfc(pts_bd_list)
  dist_bd <- st_distance(pts_sfc, pts_bd_sfc)
  
  num_nodes <- 1
  for(i in 1:nrow(coords_bd)){
    id_bd <- which(dist_bd[,i] < 100 *.Machine$double.eps )
    if(any( st_coords[id_bd,(ncol(st_coords)-1)] != 0)){
      global_id <- which(st_coords[id_bd,(ncol(st_coords)-1)] != 0)[1]
      st_coords[id_bd,(ncol(st_coords)-1)] <-  rep(global_id, times=length(id_bd))
    }else{
      st_coords[id_bd,(ncol(st_coords)-1)] <- rep(num_nodes, times=length(id_bd))
      num_nodes <- num_nodes + 1  
    }
    st_coords[id_bd, ncol(st_coords)] <- rep(1, times=length(id_bd))
  }
  
  for(sub_id in 1:length(st_geometry(x))){
    coords_sub_id <- st_coordinates(st_geometry(x)[sub_id])
    pts_sub_id_list <- lapply(as.list(as.data.frame(t(coords_sub_id[,1:2]))), st_point)
    
    pts_sub_id_sfc <- st_sfc(pts_sub_id_list)
    dist_matrix_sub_id <- st_distance(pts_sfc, pts_sub_id_sfc)
    for(j in 1:nrow(coords_sub_id)){
      id <- which(dist_matrix_sub_id[,j] < 100* .Machine$double.eps)
      if(all(st_coords[id,(ncol(st_coords)-1)] != 0 )){ next }
      if( any( st_coords[id,(ncol(st_coords)-1)] != 0) ){
        global_id <- which(st_coords[id,(ncol(st_coords)-1)] != 0)[1]
        st_coords[id,(ncol(st_coords)-1)] <-  rep(global_id, times=length(id))
      }else{
        st_coords[id,(ncol(st_coords)-1)] <- rep(num_nodes, times = length(id))
        num_nodes = num_nodes + 1
      }
    }
  }
  
  geometry <- vector(mode="list", length =length(st_geometry(x)))
  
  for(sub_id in 1:length(st_geometry(x))){
      
    nodes <- unique(st_coords[which(st_coords[,5]==sub_id),])
    n_rings <- diff(range(nodes[,4])) + 1 # number of "rings"! see st_coordinates,
                                          # 1 main ring (external boundary), 
                                          #>1 interior rings (holes!)
    
    edges <- matrix(nrow=0, ncol=2)        
    
    edges_group <- matrix(NA, nrow=0, ncol=1)
    edges_bd <- matrix(nrow=0, ncol=1)
    edges_ring <- matrix(nrow=0, ncol=1)
    nodes_group <- matrix(nrow=0, ncol=1)
    for(ring in 1:n_rings){

      nodes_ring <- nodes[which(nodes[,4] == ring),]
      # in the current main ring
      n_lines <- diff(range(nodes_ring[,3])) + 1 
      for(i in 1:n_lines){
        id_edge_sub <- which(nodes_ring[,3] == i)
        # local
        edges_sub_i <-  cbind(id_edge_sub,
                             c(id_edge_sub[2:length(id_edge_sub)],id_edge_sub[1]))
        # global 
        edges_i <- cbind(nodes_ring[edges_sub_i[,1],(ncol(st_coords)-1)], nodes_ring[edges_sub_i[,2],(ncol(st_coords)-1)])           
        
        edges_bd_i <- matrix(0, nrow=nrow(edges_sub_i),ncol=1)
        edges_group_i <- matrix(NA, nrow=nrow(edges_sub_i), ncol=1)
        
        nodes_group_i <-matrix(0, nrow(nodes[id_edge_sub,]),ncol=1)
        
        centroids_list <- lapply( as.list(as.data.frame(t(edges_sub_i))), FUN=function(edge){
          st_centroid(st_linestring(as.matrix(nodes_ring[edge,1:2])))
        })
        
        is_inside <- st_within(st_sfc(centroids_list, crs=st_crs(x)$input), 
                                    st_union(st_geometry(x)), sparse = F)
        
        for(e in 1:nrow(edges_sub_i)){
          if( nodes_ring[edges_sub_i[e,1],ncol(st_coords)] & nodes_ring[edges_sub_i[e,2],ncol(st_coords)] ){ 
            if(!is_inside[e]) edges_bd_i[e] <- 1 
            
            edges_group_i[e] <- sub_id
            nodes_group_i[edges_sub_i[e,1]] <- e
            nodes_group_i[edges_sub_i[e,2]] <- e
          } 
        }
        
        if(i == 1){ 
          edge_ring_i <- as.matrix(rep(ring, times=nrow(edges_sub_i)))
        }else{  
          edge_ring_i <- as.matrix(rep((-i+1), times=nrow(edges_sub_i)))
        }
        edges <- rbind(edges, edges_i)
        edges_ring <- rbind(edges_ring, edge_ring_i)
        edges_bd <- rbind(edges_bd, edges_bd_i)
        edges_group <- rbind(edges_group, edges_group_i)
        #edge_sub <- rbind(edge_sub, edges_sub_i)
        nodes_group <- rbind(nodes_group, nodes_group_i)
      }
    }
    storage.mode(edges) <-"integer"
    storage.mode(edges_bd) <-"integer"
    storage.mode(nodes[,ncol(st_coords)]) <- "integer" 
    geometry[[sub_id]] <- list(nodes = as.matrix(nodes[,1:2]), nodes_group = nodes_group, nodes_boundary = nodes[,ncol(st_coords)],
                               edges = edges, edges_group = edges_group, edges_boundary = edges_bd,
                               edges_ring = edges_ring)
    
  }
  
  coords <- unique(st_coords[order(st_coords[,(ncol(st_coords)-1)]),][,c(1:2,ncol(st_coords))])
  
  boundary <- coords[,3]
  storage.mode(boundary) <- "integer"
  coords <- coords[,1:2]
  
  crs <- NA
  if(!is.na(st_crs(x))) crs <- st_crs(x)$input
      .DomainCtr$new(geometry = geometry, coords=coords, boundary=boundary, crs = crs)
})

#' Delaunay triangulation of the spatial domain
#' 
#' @param domain a domain object returned by \code{Domain}.
#' @param maximum_area maximum triangle area.
#' @param minimum_angle minimum triangle angle in degrees.
#' @rdname build_mesh
#' @importFrom RTriangle pslg
#' @importFrom RTriangle triangulate
#' @importFrom sf NA_crs_
#' @importFrom sf st_point_on_surface
#' @export
setGeneric("build_mesh", function(domain, maximum_area, minimum_angle) standardGeneric("build_mesh"))

#' @rdname build_mesh
setMethod("build_mesh", signature=c("Domain", "numeric", "numeric"),
          function(domain, maximum_area, minimum_angle){
  S = matrix(nrow=0, ncol=2)
  SB = matrix(nrow=0, ncol=1)
  H = matrix(nrow=0, ncol=2)
  
  for( sub_id in 1:length(extract_private(domain)$geometry)){
    groups_id <- unique(extract_private(domain)$geometry[[sub_id]]$edges_ring)
    holes_id <- groups_id[ groups_id < 0 ]
    num_holes <- length(holes_id)
    
    if(num_holes > 0){
      for(i in 1:num_holes){
        edges_id <- which(extract_private(domain)$geometry[[sub_id]]$edges_ring==holes_id[i])
        if( ! all(as.integer(extract_private(domain)$geometry[[sub_id]]$edges_boundary[edges_id])) ){ # NON VERO buco
          next
        } 
        path  <- t(extract_private(domain)$geometry[[sub_id]]$edges[edges_id,])[1,]
        path  <- c(path,path[1]) 
        nodes <- as.matrix(extract_private(domain)$coords[path,1:2])
        holes <- st_point_on_surface(st_polygon(list(nodes)))
        H = rbind(H, holes)
      }
    }
    S = rbind(S, extract_private(domain)$geometry[[sub_id]]$edges)
    SB = rbind(SB, as.matrix(extract_private(domain)$geometry[[sub_id]]$edges_boundary))
  }
  
  if(nrow(H)==0){
    H = NA
  }
  pslg <- pslg(P = extract_private(domain)$coords[,1:2], PB = as.matrix(extract_private(domain)$boundary),
               S = S, SB = SB, H = H) 
  triangulation <- triangulate(p = pslg, a = maximum_area, q = minimum_angle)
  
  res <- Mesh(triangulation)
  set_geometry(res, extract_private(domain)$geometry)
  set_time_interval(res, extract_private(domain)$time_interval)
  set_crs(res, extract_private(domain)$crs)
  return(res)
})

#' create spatio-temporal domain
#'
#' @param op1 A mesh object created by \code{Mesh}.
#' @param op2 A numeric vector.
#' @return An S4 object representing a spatio-temporal domain.
#' @rdname Domain_times_vector
#' @export 
setGeneric("%X%", function(op1, op2) standardGeneric("%X%"))

#' @rdname Domain_times_vector
setMethod("%X%", signature=c(op1="Domain", op2="numeric"),
          function(op1, op2){
            if(op2[1] > op2[length(op2)])
              stop("Error! First time instant is greater than last time instant.")
            set_time_interval(op1, op2)
            op1          
})