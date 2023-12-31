library(femR)

## load domain data and generate mesh object
data("unit_square", package="femR")

mesh = Mesh(unit_square)
class(mesh)
plot(mesh)

# create Functional Space
fe_order = 1
Vh <- FunctionSpace(mesh, fe_order)

exact_solution <- function(points){
  return( sin(2. * pi * points[,1]) * sin(2. * pi * points[,2]) )
}

## define differential operator in its strong formulation
f <- Function(Vh)
L <- -laplace(f)
## forcing term
u <- function(points){
    return(8.*pi^2* sin( 2.* pi * points[,1]) * sin(2.*pi* points[,2]) ) 
}
## Dirichlet BC
g <- function(points){
  return(rep(0, times=nrow(points)))
}
## Pde constructor
pde <- Pde(L,u)
pde$set_boundary_condition(fun=0.,
                           type = "dirichlet")

## solve problem
pde$solve()

## compute L2 norm of the error
u_ex <- as.matrix(exact_solution(pde$get_dofs_coordinates()))
error.L2 <- sqrt(sum(pde$get_mass() %*% (u_ex - pde$get_solution())^2))
cat("L2 error = ", error.L2, "\n")

## perform evaluation at single point
point = c(0.2, 0.5)
f$eval_at(point)

## evaluate over a 10x10 grid
x <- seq(0, 1, length.out = 50)
y <- x
points <- expand.grid(x, y)
max(abs( f$eval_at(points) - as.matrix(exact_solution(points))))

evaluations <- max(abs(f$eval_at(mesh$get_nodes()) - exact_solution(mesh$get_nodes())))
## plot solution 
options(warn=-1)
plot(f) %>% layout(scene=list(aspectmode="cube")) %>% hide_colorbar()

plot(f)
contour(f)
