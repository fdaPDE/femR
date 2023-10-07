library(femR)

## 1. Defining domain
data("unit_square", package="femR")
mesh <- Mesh(unit_square)

# Functional Space
Vh <- FunctionSpace(mesh, fe_order=2)

## 2. Defining the solution of the PDE
f <- Function(Vh)
## 3. Defining the differential operator
L <- -laplace(f) 
## 4. Defining the forcing term and the Dirichlet boundary conditions as standard R functions
# forcing term
u <- function(points){
  return(8.*pi^2* sin(2.*pi*points[,1])*sin(2.*pi*points[,2])) 
}
# Dirichlet BC
dirichletBC <- function(points){
  return(matrix(0,nrow=nrow(points), ncol=1))
}
## 5. Building the PDE object
pde <- Pde(L, u, dirichletBC)
## 6. computing the discrete solution
pde$solve()
## 7. Plots
plot(f)

## perform evaluation at single point
point = c(0.2, 0.5)
f$eval_at(point)

## evaluate over a 10x10 grid
x <- seq(0, 1, length.out = 50)
y <- x
points <- expand.grid(x, y)
tmp <- f$eval_at(points)

max(f$eval_at(points) - pde$eval(Vh$mesh, as.matrix(f$coeff), as.matrix(points)))
