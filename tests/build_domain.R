# reading sf object -----------------------------------------------------------
library(femR)
library(sf)
library(mapview)

data("franconia",package = "mapview")
class(franconia)
franconia_boundary <- st_geometry(franconia)
plot(franconia_boundary)
domain <- Domain(franconia_boundary)
mesh <- Triangulate(domain, maximum_area = 0.01, minimum_angle = 20) 
mesh$set_crs(4326)
plot(mesh)

# reading pslg object ----------------------------------------------------------
library(RTriangle)
pslg <- pslg(P=rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1)),
          S=rbind(c(1, 2), c(2, 3), c(3, 4), c(4,1)))
domain <- Domain(pslg)
domain$crs
domain %X% c(0,1)

mesh <- Triangulate(domain, maximum_area = 0.1, minimum_angle = 20)
plot(mesh)
mesh$crs
mesh$times
mesh$set_deltaT(0.1)
