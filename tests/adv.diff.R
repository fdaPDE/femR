library(femR)

## load domain data and generate mesh object
data("unit_square", package="femR")

mesh = Mesh(unit_square)
class(mesh)
plot(mesh)

# create Functional Space
fe_order = 1
Vh <- FunctionSpace(mesh, fe_order)

W_ <- 1.
R_ <- 1.
H_ <- 1.
beta_ <- 1.

alpha_ <- H_ * beta_ / R_
gamma_ <- pi * W_ / R_  

lambda1 <- -alpha_/2 - sqrt((alpha_/2)^2 + pi^2)
lambda2 <- -alpha_/2 + sqrt((alpha_/2)^2 + pi^2)

p_ <- (1-exp(lambda2))/(exp(lambda1)-exp(lambda2))

exact_solution <- function(points){
  return( -gamma_/pi^2 * ( p_ * exp(lambda1 * points[,1]) + (1 - p_) * exp(lambda2 * points[,1]) - 1. ) * sin(pi * points[,2] ) )
}

## define differential operator in its strong formulation
f <- Function(Vh)
L <- -laplace(f) + dot(c(alpha_,0), grad(f))
## forcing term
u <- function(points){
  return(gamma_ * sin(pi * points[,2])) 
}
## Dirichlet BC
dirichletBC <- function(points){
  return(matrix(0,nrow=nrow(points), ncol=1))
}
## create pde
pde <- Pde(L, u, dirichletBC)

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
max(abs( f$eval_at(points) - as.matrix(exact_solution(points))))
# plot solution 
options(warn=-1)

plot(f)
plot(f) %>% layout(scene=list(aspectmode="cube"))
plot(f, colorscale="Rainbow") 
plot(f, colorscale="Jet") %>% layout(scene=list(aspectmode="cube"), 
                                     title="Solution")
plot(f, colorscale="Electric") %>% 
  layout(scene=list(aspectmode="cube"), 
         title="Solution") %>%
  hide_colorbar()

contour(f)
contour(f, colorscale="Jet")
# see: https://plotly.com/r/reference/mesh3d/