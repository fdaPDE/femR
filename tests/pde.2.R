library(femR)

## load domain data and generate mesh object
data("unit_square", package="femR")
fe_order = 2

# create Functional Space
Vh <- FunctionalSpace(unit_square, fe_order = 2)
plot(Vh$mesh)

exact_solution <- function(points){
  return( sin(2. * pi * points[,1]) * sin(2. * pi * points[,2]) )
}

## define differential operator in its strong formulation
f <- feFunction(Vh)
L <- -laplace(f)
## forcing term
u <- function(points){
    return(8.*pi^2* sin( 2.* pi * points[,1]) * sin(2.*pi* points[,2]) ) 
}
## Dirichlet BC
dirichletBC <- function(points){
  return(rep(0, times=nrow(points)))
}
## create pde
pde <- pde(L, u, dirichletBC)

## solve problem
pde$solve()

## compute L2 norm of the error
u_ex <- as.matrix(exact_solution(pde$get_dofs_coordinates()))
error.L2 <- sqrt(sum(pde$get_mass() %*% (u_ex - pde$solution())^2))
cat("L2 error = ", error.L2, "\n")

## perform evaluation at single point
point = c(0.2, 0.5)
f$eval_at(point)

## evaluate over a 10x10 grid
x <- seq(0, 1, length.out = 50)
y <- x
points <- expand.grid(x, y)
f$eval_at(points)

## plot solution 
options(warn=-1)
plot(f) %>% layout(scene=list(aspectmode="cube")) %>% hide_colorbar()

contour(f)
