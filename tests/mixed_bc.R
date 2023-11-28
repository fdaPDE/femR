library(RTriangle)
library(femR)
library(sf)
options(warn = -1)
p <- pslg(P=rbind(c(0, 0),  # first node
                  c(1, 0),  # ...   node
                  c(1, 1),  # ...   node
                  c(0, 1)), # ...   node
          S=rbind(c(1, 2),  # first physical edge 
                  c(2, 3),  # ...   physical edge
                  c(3, 4),  # ...   physical edge
                  c(4,1)),  # ...   physical edge
          SB = cbind(c(1,1,1,1)), 
          PB = cbind(c(1,1,1,1))) 

domain <- Domain(p)
plot(st_as_sfc(domain))
domain_sf <- st_as_sfc(domain)
plot(domain_sf)

triangulation <- triangulate(pslg, a = 0.0025, q=20)
plot(triangulation)

num_bd_nodes <- sum(triangulation$PB==1)
bd_id <- which(triangulation$PB==1)
pts_list <- list()
for(i in 1:num_bd_nodes){
  pts_list[[i]] <- st_point( t(as.matrix(triangulation$P[bd_id[i],1:2])))
}
 
triangulation$PB <- matrix(0, nrow=nrow(triangulation$P),ncol=1) 
edge <- st_linestring(as.matrix(triangulation$P[pslg$S[4,],1:2]))
dirichlet_nodes <- st_intersects(edge, 
                                 st_sfc(pts_list),sparse = FALSE)
dirichlet_id <- which(dirichlet_nodes[1,] == T)
triangulation$PB[ bd_id[dirichlet_id] ] = 1

plot(triangulation)
points(triangulation$P[triangulation$PB==1,], pch=16, col="red")

## 1. building mesh
times <- seq(0,3,by=0.05)
nodes <- triangulation$P; edges <- triangulation$T; boundary  <- triangulation$PB
mesh <- Mesh(domain=list(nodes=nodes, elements=edges, boundary=boundary))
mesh <- mesh %X% times
plot(mesh)

## 2. Defining the solution of the PDE
Vh <- FunctionSpace(mesh, fe_order=1) 
u <- Function(Vh)
## 3. Defining the differential operator
Lu <- dt(u) -laplace(u) + dot(c(0.25,0), grad(u)) 
## 4. Defining the forcing term and Dirichlet boundary conditions
##    as standard R functions
# forcing term
f <- function(points, times){
  matrix(0, nrow=nrow(points), ncol=length(times))
}
Q <- 10
dirichlet_id <- which(mesh$get_boundary()==1)

plot(st_as_sfc(mesh))
points(mesh$get_nodes()[dirichlet_id,], pch=16)

# Dirichlet boundary conditions 
dirichletBC <- function(points, times){
  res <- matrix(0, nrow=nrow(points), ncol=length(times))
  for( t in 1:length(times)){
    res[dirichlet_id,t] <- rep(Q*exp(-times[t]), times=length(dirichlet_id))
  }
  return(res)
}

initialCondition <- function(points){
  res <- matrix(0, nrow=nrow(points))
  res[dirichlet_id,] <- rep(Q, times=length(dirichlet_id))
  return(res)
}
## 5. Building the PDE object
pde <- Pde(Lu, f)
# setting initial conditions
pde$set_initialCondition(initialCondition)

# setting boundary conditions
pde$set_dirichletBC(dirichletBC)
## 7. computing the discrete solution
pde$solve()

## 8. Plots
contour(u)
