library(femR)

## load domain data and generate mesh object
data("unit_square", package="femR")

mesh = Mesh(unit_square)
class(mesh)
plot(mesh)

# create Functional Space
fe_order = 2
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
pde$set_boundary_condition(boundary_condition =0.,
                           type = "dirichlet")

## solve problem
pde$solve()

## compute L2 norm of the error
u_ex <- as.matrix(exact_solution(pde$dofs_coordinates()))
error.L2 <- sqrt(sum(pde$mass() %*% (u_ex - f$coefficients())^2))
cat("L2 error = ", error.L2, "\n")

## perform evaluation at single point
point = c(0.2, 0.5)
f$eval(point)

## evaluate over a 10x10 grid
x <- seq(0, 1, length.out = 50)
y <- x
points <- expand.grid(x, y)
max(abs( f$eval(points) - as.matrix(exact_solution(points))))

evaluations <- max(abs(f$eval(mesh$nodes()) - exact_solution(mesh$nodes())))
## plot solution 
options(warn=-1)
plot(f) %>% layout(scene=list(aspectmode="cube")) %>% hide_colorbar()

plot(f)
contour(f)

### 

basis <- Vh$get_basis()
Psi <- basis$eval(basis$get_dofs_coordinates())
dim(Psi)

A  <- basis$eval(as.matrix(points))
dim(A)

R0 <- pde$get_mass()
R1 <- pde$get_stiff()
n <- basis$size()
lambda <- 1.
row_1 <- cbind( 1/n * t(Psi)%*%Psi, lambda*t(R1))
row_2 <- cbind(lambda*R1, R0)

M <- rbind(row_1,
           row_2)
dim(M)

y <- u_ex()