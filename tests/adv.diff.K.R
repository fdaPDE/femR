## load domain data and generate mesh object
library(femR)

rm(list=ls())
graphics.off()

data("unit_square", package="femR")
unit_square = create_mesh(data = unit_square)
K <- matrix(c(1,0,0,1), byrow=T, nrow=2,ncol=2)

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
f <- Function(domain = unit_square)
#L <- -laplace(f) + dot(c(alpha_,0), grad(f))
L <- -div(K*grad(f)) + dot(c(alpha_,0), grad(f))
## forcing term
u <- function(points){
  return(gamma_ * sin(pi * points[,2])) 
}
## create pde
pde <- pde(L, u, fe_order = 1)

## set boundary conditions
nodes <- pde$get_dofs_coordinates()
dirichletBC <- as.matrix(rep(0., times = dim(nodes)[1]))
pde$set_dirichlet_bc(dirichletBC)

## solve problem
pde$solve()

u_ex <- as.matrix(exact_solution(nodes))
error.L2 <- sqrt(sum(pde$get_mass() %*% (u_ex - pde$solution())^2))
cat("L2 error = ", error.L2, "\n")