# reading list -----------------------------------------------------------------
# 1.
library(femR)
library(sf)
library(dplyr)
nodes_ext <- rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
nodes_int <- rbind(c(0.25, 0.25), c(0.75, 0.75), c(0.75, 0.25), c(0.25, 0.75))
edges_ext <- rbind(c(1, 2), c(2, 3), c(3, 4), c(4,1))
edges_int <- rbind(c(5,7), c(7,6), c(6,8), c(8,5))
# hole
input <- list(nodes = rbind(nodes_ext, nodes_int), 
              edges = rbind(edges_ext, edges_int), 
              edges_group = c(rep(1, times=nrow(edges_ext)),
                              rep(-1, times=nrow(edges_int))))
# no hole
# input <- list(nodes = nodes_ext,
#               edges = edges_ext,
#               edges_group = rep(1, times=nrow(edges_ext)))
domain <- Domain(input)

domain_sf <- st_as_sfc(domain)
plot(domain_sf, col="red")
st_centroid(domain_sf)
st_point_on_surface(domain_sf) # :)
st_crs(domain_sf)

domain_sf <- st_as_sf(domain) # + dataset
plot(domain_sf)
domain_sf %>% filter(label=="edge") %>% 
              select(local_id)
plot(domain_sf %>% filter(label=="edge") %>% 
       select(local_id), lwd=3)

plot(domain_sf %>% filter(label=="edge" & local_id >3) %>% 
       select(local_id), lwd=3)
#plot(domain_sf %>% select(bc) %>% filter(!is.na(bc)))
#plot(domain_sf %>% select() %>% filter(!is.na(bc)))

mesh <- build_mesh(domain, maximum_area = 0.05, minimum_angle = 20)
plot(mesh)
mesh_sf <- st_as_sfc(mesh)
plot(mesh_sf, col="red")

# espolare st_triangulate / st_triangulate_constrained
plot( st_triangulate(mesh_sf), col="red" )

# 2.
nodes_ext <- rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
nodes_int1 <- rbind(c(0.125, 0.125), c(0.25, 0.125), c(0.25, 0.25), c(0.125, 0.25))
nodes_int2 <- rbind(c(0.75, 0.75), c(0.975, 0.75), c(0.975, 0.975), c(0.75, 0.975))
edges_ext <- rbind(c(1, 2), c(2, 3), c(3, 4), c(4,1))
edges_int1 <- rbind(c(5,6), c(6,7), c(7,8), c(8,5))
edges_int2 <- rbind(c(9,10), c(10,11), c(11,12), c(12,9))
# hole
input <- list(nodes = rbind(nodes_ext, nodes_int1,  nodes_int2), 
              edges = rbind(edges_ext, edges_int1, edges_int2), 
              edges_group = c(rep(1, times=nrow(edges_ext)),
                              rep(-1, times=nrow(edges_int1)),
                              rep(-2, times=nrow(edges_int2))))

domain <- Domain(input)
domain_sf <- st_as_sfc(domain)
plot(domain_sf,col="red")

mesh <- build_mesh(domain, maximum_area = 0.00125, minimum_angle = 20)
plot(mesh)
mesh_sf <- st_as_sfc(mesh)
plot(mesh_sf, col="red")



# reading sf object ------------------------------------------------------------
library(femR)
library(sf)
library(mapview)

data("franconia",package = "mapview")
class(franconia)
franconia_bd <- st_geometry(franconia)
domain <- Domain(franconia_bd)

domain_sfc <- st_as_sfc(domain)
plot(domain_sfc)

domain_sf <- st_as_sf(domain) # + dataset
domain_sf[[2]] %>% 
        select()

mesh <- build_mesh(domain, maximum_area = 0.001, minimum_angle = 10)
plot(mesh) # :)

# vecchia versione 
domain <- Domain(st_union(st_geometry(franconia)))
plot(st_as_sf(domain))
mesh <- build_mesh(domain, maximum_area = 0.001, minimum_angle = 10)
plot(mesh)
#domain <- Domain(franconia_bd)
#mesh <- build_mesh(domain, maximum_area = 0.01, minimum_angle = 20) 
#plot(mesh)
#st_crs(mesh) # :)

franconia_bd <- st_cast(st_union(franconia_bd), to="LINESTRING")
plot(franconia_bd)

sub_domain_on_boundary <- st_intersects(franconia_bd, st_geometry(franconia), sparse = F)
plot(franconia_bd)
plot(st_geometry(franconia)[which(sub_domain_on_boundary==T)], add=T, col="red")

adj_matrix <- st_intersects(st_geometry(franconia),
                            st_geometry(franconia))

# esempio 
plot(st_geometry(franconia))
plot(st_geometry(franconia)[1], col="red", add=T)
plot(st_geometry(franconia)[adj_matrix[[1]]], col="blue", add=T)
data <- franconia 
data$sub_domain_id <- rep(1:nrow(franconia))
data <- data %>% dplyr::select(sub_domain_id)
plot(data, reset=FALSE)
plot(data[1,], col="red", add=T)
plot(data[2,], col="blue", add=T)
plot(data[3,], col="red", add=T)
plot(data[5,], col="gray", add=T) #5 ok
plot(data[6,], col="yellow", add=T)
plot(data[7,], col="gray", add=T)
plot(data[8,], col="green3", add=T) # 8 ok

data <- data %>% dplyr::filter(sub_domain_id %in% c(5,8)) # pazzesko
tmp <- st_intersection(st_geometry(data))

data_sfc <- st_geometry(data)
st_coordinates(data_sfc)
plot(data_sfc)
coords <- data_sfc[1]

# ------------------------------------------------------------------------------

coords <- st_coordinates(franconia)
coords_unique <- unique(coords[,1:2]) # global numbering

tmp <- unlist(apply(coords[,1:2], MARGIN=1, FUN=function(x){ 
  which( apply(coords_unique-x, MARGIN=1, FUN=function(y){Matrix::norm(y, type="2")})  
         < 10* .Machine$double.eps) })) 

tmp <- which(apply(coords, 1, function(x) all(x %in% coords_unique))) 
st_cast(franconia_boundary, to="POLYGON", group_or_split = F)
tmp2 <- st_cast(st_cast(franconia_boundary, to="POLYGON", group_or_split = F), 
                to="MULTILINESTRING", group_or_split=T)
tmp3 <- st_cast(tmp2, to="LINESTRING")

tmp <- st_cast(st_cast(st_cast(franconia_boundary, 
                              to="POLYGON"), 
                              to="MULTILINESTRING"), 
                              to="LINESTRING")

franconia_intersection <- st_intersection(franconia_boundary)

plot(st_triangulate_constrained(franconia_boundary))
tmp <- st_triangulate_constrained(franconia_boundary)
tmp
tmp[1]
# reading pslg object ----------------------------------------------------------
library(RTriangle)
library(femR)
library(sf)
nodes <- rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
edges <- rbind(c(1, 2), c(2, 3), c(3, 4), c(4,1))
pslg <- pslg(P=nodes, S=edges, PB=rep(1,nrow(nodes)), SB=rep(1,nrow(edges)))
domain <- Domain(pslg)
mesh <- build_mesh(domain, maximum_area = 0.1, minimum_angle = 20)
plot(mesh)

# pslg <- pslg(P=nodes, PB=c(1,0,0,1))
# domain <- Domain(pslg)
# plot(st_as_sf(domain))
# points(nodes[domain$geometry[[1]]$nodes_boundary == 1,], pch=16, col="red")

#mesh <- build_mesh(domain, maximum_area = 0.0125, minimum_angle = 20)
#plot(st_as_sf(mesh))
#points(mesh$get_nodes()[mesh$get_boundary()==1,], pch=16, col="red")
#points(mesh$get_nodes()[mesh$get_boundary()==0,], pch=16, col="blue")

# st_as_sf(c) 
#x<- domain
# results <- data.frame(id=)
TMP_FUNCTION<-function(x){

  geom_sf <- list()
  for(sub_id in 1:length(x$geometry)){
    label <- vector(mode="character")
    sub <- vector(mode="integer")
    local_id <- vector(mode="integer")
    bc <- vector(mode="character")
    group <- vector(mode="numeric")
    boundary <- vector(mode="integer")
    path_list = vector(mode="list")
    edge_list = vector(mode="list")
    pts_list = vector(mode="list")
    
    geometry <- x$geometry[[sub_id]]
    path_id <- unique(geometry$edges_group)
    n_physical_lines <- length(path_id)
    
    path_list_sub = vector(mode="list", length=n_physical_lines)
    for(i in 1:n_physical_lines){
      path  <- t(geometry$edges[geometry$edges_group==path_id[i],][,1]) # )[1,]
      #if(is.vector(path)) path <- t(as.matrix(path))
      #nodes <- matrix(nrow=0,ncol=2)
      #for(j in 1:nrow(path)){
      #  nodes <- rbind(nodes, as.matrix(rbind(x$coords[path[j,1],1:2],
      #                                        x$coords[path[j,2],1:2])))
      #}
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
      pts_list_sub[[i]] <- st_point(as.matrix(node))
    }
  # nodes_sf <- st_sf(data.frame(nodes_group=as.factor(nodes_group)), geometry=(pts_list))
  # plot(nodes_sf)
    label= append(label, c(rep("path", times=n_physical_lines),
                           rep("edge", times=nrow(geometry$edges)),
                           rep("node", times=nrow(geometry$nodes))))
    group = append(group, c(path_id , geometry$edges_group, geometry$nodes_group))
    local_id = append(local_id, 
                      c(1:n_physical_lines,1:nrow(geometry$edges), 1:nrow(geometry$nodes)))
    
    bc = append(bc, 
                c(rep(NA, times = n_physical_lines), geometry$BC, rep(NA, times=nrow(geometry$nodes))))
    sub = append(sub, rep(sub_id, 
                          times=(n_physical_lines + nrow(geometry$edges) + nrow(geometry$nodes))))  
    path_list = append(path_list, path_list_sub) 
    edge_list = append(edge_list, edge_list_sub)
    pts_list  = append(pts_list, pts_list_sub)
  
    group <- as.factor(group)
    geom_sfc <- st_sfc( append(append(path_list, edge_list), pts_list))
    geom_sf[[sub_id]] <- st_as_sf(data.frame(label= label,
                                             group = group, local_id=local_id, bc),
                                  geometry = geom_sfc)
    
  }
  return(geom_sf)
}

geom_sf <- TMP_FUNCTION(domain)
geom_sf[[1]]
edge_sf <- geom_sf[[1]] %>% dplyr::filter(label == "edge") #%>% dplyr::select(group) 
plot(edge_sf)
nodes_sf <- geom_sf %>% dplyr::filter(label == "node") %>% dplyr::select(group)

# ------------------------------------------------------------------------------

set_group <- function(domain, local_id, value, sub_id=1){
    domain$geometry[[sub_id]]$edges_group[local_id] <- rep(value, times=length(local_id))
    nodes_id <- unique(as.vector(domain$geometry[[sub_id]]$edges[local_id,]))
    domain$geometry[[sub_id]]$nodes_group[nodes_id] <- rep(value, length(nodes_id))
    return(domain)
}
domain<-x
domain <- set_group(domain, local_id = c(1,2,3,4), value=0) 

domain <- set_group(domain, local_id = 4, value = 1)
geom_sfc <- TMP_FUNCTION(domain)
plot(geom_sfc %>% dplyr::filter(label=="path"))
plot(geom_sfc)

# prova franconia
library(mapview)
data("franconia")

domain<- Domain(st_geometry(franconia))

geom_franconia <- TMP_FUNCTION(domain)
sub_id <- 1
sub_domain_1 <- geom_franconia[[sub_id]] 
plot(sub_domain_1)

plot( geom_franconia %>% filter(sub_id == 1) )
# - dplyr --------------------------------------------------
library(femR)
library(sf)
library(dplyr)
data("franconia", package="mapview")
x <- st_geometry(franconia)

n_sub_domains <- length(x)
nodes <- vector(mode="list", length=n_sub_domains)
edges <- vector(mode="list", length=n_sub_domains)
nodes_group <- vector(mode="list", length=n_sub_domains)
edges_group <- vector(mode="list", length=n_sub_domains)
geometry <- vector(mode="list", length=n_sub_domains)
for(n in 1:n_sub_domains){
  nodes <- unique(st_coordinates(x[1]))[,1:2]
  edges <- cbind(1:nrow(nodes),c(2:nrow(nodes),1))
  nodes_group <- rep(n, times=nrow(nodes))
  edges_group <- rep(n, times=nrow(edges))
  geometry <- list(nodes = nodes, edges = edges,
                   nodes_group = nodes_group, edges_group = edges_group)
  crs <- NA_crs_
}




groups_id <- unique(x$geometry$edges_group)
n_physical_lines <- length(groups_id)

groups_pts_id <- unique(x$geometry$nodes_group)
n_physical_pts <- length(groups_pts_id)

path_list = vector(mode="list", length=n_physical_lines)
for(i in 1:n_physical_lines){
  path  <- t(x$geometry$edges[x$geometry$edges_group==groups_id[i],])[1,]
  path  <- c(path,path[1]) 
  nodes <- x$geometry$nodes[path,]
  path_list[[i]] <- st_linestring(nodes)
}
# edges_sf <- st_sf(data.frame(edges_group=as.factor(edges_group)), geometry=(path_list))
# plot(edges_sf)

pts_list = vector(mode="list", length=n_physical_pts)
for(i in 1:nrow(x$geometry$nodes)){
  node  <- x$geometry$nodes[i,]
  pts_list[[i]] <- st_point(node)
}
# nodes_sf <- st_sf(data.frame(nodes_group=as.factor(nodes_group)), geometry=(pts_list))
# plot(nodes_sf)

geom_sfc <- st_sfc( append(path_list, pts_list))
geom_sf <- st_as_sf(data.frame(label= 
                                 c(rep("edge", times=nrow(domain$geometry$edges)),
                                   rep("node", times=nrow(domain$geometry$nodes))),
                               physical_group = c(as.factor(edges_group), 
                                                  as.factor(nodes_group))),
                    geometry = geom_sfc)

edge_sf <- geom_sf %>% dplyr::filter(label == "edge") %>% dplyr::select(physical_group) 
nodes_sf <- geom_sf %>% dplyr::filter(label == "node") %>% dplyr::select(physical_group)

plot(edge_sf, lwd=3)
plot(nodes_sf, cex=2, pch=16)



seg <- st_linestring(matrix(c(0,0,1,0),nrow=2,ncol=2,byrow=T))

st_intersects(seg, st_sfc(list( st_point(c(0,0)), st_point(c(1,0)))),sparse = FALSE)
