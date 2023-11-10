library(RTriangle)
library(femR)

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

plot(p)
unit_square <- triangulate(p, a = 0.00125, q=30)
# select nodes which do not belong to 1 and 4 physical edges
mask <- unit_square$PB != 1 & unit_square$PB != 4  

# Homogeneous Neumann BC will be imposed on nodes belonging to 2 and 3
unit_square$PB[mask,] = 0 

# Dirichlet BC will be imposed on nodes belonging to 1 and 4
unit_square$PB[!mask,] = 1

plot(unit_square)
points(unit_square$P[unit_square$PB==1,],pch=16,col="red") # dirichlet

# create list to be passed to Mesh()
domain <- list( elements = unit_square$T, 
                nodes = unit_square$P, boundary = unit_square$PB)

## 1. building mesh
mesh <- Mesh(domain)
plot(mesh)

## 2. Defining the solution of the PDE
Vh <- FunctionSpace(mesh, fe_order=1) 
u <- Function(Vh)
## 3. Defining the differential operator
Lu <- -laplace(u) 
## 4. Defining the forcing term and Dirichlet boundary conditions
##    as standard R functions
# forcing term
f <- function(points){
  return(4*((points[,1]^2+points[,2]^2) * sin(points[,1]^2 + points[,2]^2 - 1) 
            - cos(points[,1]^2 + points[,2]^2 - 1))) 
}
# Dirichlet boundary conditions 
dirichletBC <- function(points){
  return(matrix(0,nrow=nrow(points), ncol=1))
}
## 5. Building the PDE object
pde <- Pde(Lu, f)
# setting boundary conditions
pde$set_dirichletBC(dirichletBC)
## 7. computing the discrete solution
pde$solve()

## 8. Plots
plot(u)
