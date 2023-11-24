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
          SB = cbind(c(1,2,3,4))) # use SB to specify "physical boundaries"
# 1 bottom, 2 right, 3 up, 4 left 
unit_square <- triangulate(p, a = 0.00125, q=30)

plot(p)
points(unit_square$P[unit_square$PB==1,],pch=16,col="red") # dirichlet
points(unit_square$P[unit_square$PB==2,],pch=16, col="blue")
points(unit_square$P[unit_square$PB==3,],pch=16, col="green3")
points(unit_square$P[unit_square$PB==4,],pch=16, col="black")
points(unit_square$P[1,1],unit_square$P[1,2],
       pch=8,col="black") # dirichlet
points(unit_square$P[4,1],unit_square$P[4,2],
       pch=8,col="black") # dirichlet

# select nodes which do not belong to 1 and 4 physical edges
mask <- unit_square$PB != 4  

# Homogeneous Neumann BC will be imposed on nodes belonging to 1,2,3
unit_square$PB[mask,] = 0 

# Dirichlet BC will be imposed on nodes belonging to 4
unit_square$PB[!mask,] = 1
unit_square$PB[1,] = 1 # basso a sx
unit_square$PB[4,] = 1 # alto  a sx

plot(unit_square)
points(unit_square$P[unit_square$PB==1,],pch=16,col="black") # dirichlet

# create list to be passed to Mesh()
domain <- list( elements = unit_square$T, 
                nodes = unit_square$P, boundary = unit_square$PB)

times <- seq(0,3,by=0.05)
## 1. building mesh
mesh <- Mesh(domain) %X% times
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

plot(st_as_sf(mesh))
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
