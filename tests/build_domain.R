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
              edges_ring = c(rep(1, times=nrow(edges_ext)),
                              rep(-1, times=nrow(edges_int))),
              edge_group = c(rep(0, times=nrow(edges_ext)),
                              rep(1, times=nrow(edges_int)))) 
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
              edges_ring = c(rep(1, times=nrow(edges_ext)),
                              rep(-1, times=nrow(edges_int1)),
                              rep(-2, times=nrow(edges_int2))))

domain <- Domain(input)
domain_sf <- st_as_sfc(domain)
plot(domain_sf,col="red")

mesh <- build_mesh(domain, maximum_area = 0.00125, minimum_angle = 20)
plot(mesh)
mesh_sf <- st_as_sfc(mesh)
plot(mesh_sf, col="red")

# 3. reading sf object ---------------------------------------------------------
library(femR)
library(sf)
library(mapview)

data("franconia",package = "mapview")
class(franconia)
franconia_bd <- st_geometry(franconia)
plot(franconia_bd)
domain <- Domain(franconia_bd)

domain_sfc <- st_as_sfc(domain)
plot(domain_sfc)

domain_sf <- st_as_sf(domain) # + dataset

mesh <- build_mesh(domain, maximum_area = 0.0015, minimum_angle = 20)
plot(mesh) # :)

mesh_sf <- st_as_sfc(mesh)
plot(mesh_sf)
# solo boundary  
domain <- Domain(st_union(st_geometry(franconia)))
plot(st_as_sf(domain))
mesh <- build_mesh(domain, maximum_area = 0.001, minimum_angle = 10)
plot(mesh)

franconia_bd <- st_cast(st_union(franconia_bd), to="LINESTRING")
plot(franconia_bd)

sub_domain_on_boundary <- st_intersects(franconia_bd, st_geometry(franconia), sparse = F)
plot(franconia_bd)
plot(st_geometry(franconia)[which(sub_domain_on_boundary==T)], add=T, col="red")

adj_matrix <- st_intersects(st_geometry(franconia),
                            st_geometry(franconia))

# esempio ----------------------------------------------------------------------
nodes <- matrix(c(0,0,1,0,1,1,0,1,0.5,0,0.5,1), ncol=2, byrow=T)

poly1 <- st_cast(st_cast(st_linestring(nodes[c(1,5,6,4,1),]),
                         to="POLYGON"), to="MULTIPOLYGON")
poly2 <- st_cast(st_cast(st_linestring(nodes[c(5,2,3,6,5),]),
                         to="POLYGON"), to="MULTIPOLYGON")
geom_sfc <- st_sfc(list(poly1, poly2)) 
domain <- Domain(geom_sfc)
plot(st_as_sfc(domain))
st_as_sf(domain)

tmp <- st_as_sf(domain)
mesh <- build_mesh(domain, maximum_area = 0.0025, minimum_angle = 20)
plot(st_as_sfc(mesh))
points(mesh$get_nodes()[mesh$get_boundary()==1,], pch=16, col="red") # da sistemare ?
# ------------------------------------------------------------------------------

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

geom_sf <- st_as_sf(domain)
geom_sf[[1]]
edge_sf <- geom_sf[[1]] %>% dplyr::filter(label == "edge") #%>% dplyr::select(group) 
plot(edge_sf)
nodes_sf <- geom_sf %>% dplyr::filter(label == "node") %>% dplyr::select(group)


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
