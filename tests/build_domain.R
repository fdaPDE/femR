# reading list -----------------------------------------------------------------
# 1.
library(femR)
library(sf)
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

domain_sf <- st_as_sf(domain)
plot(domain_sf, col="red")
st_centroid(domain_sf)
st_point_on_surface(domain_sf) # :)

mesh <- build_mesh(domain, maximum_area = 0.05, minimum_angle = 20)
plot(mesh)
mesh_sf <- st_as_sf(mesh)
plot(mesh_sf, col="red")

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
domain_sf <- st_as_sf(domain)
plot(domain_sf,col="red")

mesh <- build_mesh(domain, maximum_area = 0.00125, minimum_angle = 20)
plot(mesh)
mesh_sf <- st_as_sf(mesh)
plot(mesh_sf, col="red")

# reading sf object ------------------------------------------------------------
library(femR)
library(sf)
library(mapview)

data("franconia",package = "mapview")
class(franconia)
franconia_boundary <- st_geometry(franconia)
plot(franconia_boundary)
domain <- Domain(franconia_boundary)
mesh <- build_mesh(domain, maximum_area = 0.01, minimum_angle = 20) 
plot(mesh)
st_crs(mesh) # :)

# reading pslg object ----------------------------------------------------------
library(RTriangle)
nodes <- rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
edges <- rbind(c(1, 2), c(2, 3), c(3, 4), c(4,1))
pslg <- pslg(P=nodes, S=edges, PB=rep(1,nrow(nodes)), SB=rep(1,nrow(edges)))
domain <- Domain(pslg)
domain$crs
domain %X% c(0,1)

mesh <- build_mesh(domain, maximum_area = 0.1, minimum_angle = 20)
plot(mesh)
mesh$crs
mesh$times
mesh$set_deltaT(0.1)
mesh$times
