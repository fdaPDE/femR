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
domain_sf <- st_as_sf(domain)
plot(domain_sf)

domain_sf %>% filter(label=="edge") %>%
  select(local_id)

plot(domain_sf %>% filter(label=="edge") %>%
       select(local_id), lwd=3)

domain$set_boundary_type(on = list(edges_id=c(4)), 
                         type = "dirichlet")
domain_sf <- st_as_sf(domain)
plot(domain_sf)

## 1. building mesh
times <- seq(0,3,by=0.05)
mesh <- build_mesh(domain, maximum_area = 0.00125, minimum_angle = 20)
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
