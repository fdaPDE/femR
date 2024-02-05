library(femR, quietly = TRUE)

## load domain data and generate mesh object
data("unit_square", package="femR")

t0 = 0.
t_max = 1
dT = 0.05
M = (t_max - t0)/dT + 1
times = seq(t0,t_max,length=M)

# Spatio-temporal domain ( Mesh inherits from DOMAIN :)  )
mesh = Mesh(unit_square) %X% c(t0, t_max)  #times

mesh$set_time_step(dT) # NUOVO
class(mesh)
plot(mesh)

# create Functional Space
Vh <- FunctionSpace(mesh)

exact_solution <- function(points,times){
  res <- matrix(0, nrow=nrow(points), ncol=length(times))
  for( t in 1:length(times)){
    res[,t] = sin( 2.* pi * points[,1]) * sin(2.*pi* points[,2]) * exp(-times[t])
  }
  return(res)
}

## define differential operator in its strong formulation
u <- Function(Vh)

## define differential operator in its strong formulation
L <- dt(u) -laplace(u) 

#L <- dt(u) - laplace(u) + dot(c(0,1),grad(u))

## forcing term
f <- function(points,times){
  res <- matrix(0, nrow=nrow(points), ncol=length(times))
  for( t in 1:length(times)){
    res[,t] = (8*pi^2 -1) * sin( 2.* pi * points[,1]) * sin(2.*pi* points[,2]) * exp(-times[t])
  }
  return(res) 
}
## Dirichlet BC
g <- function(points, times){
  return(matrix(0, nrow=nrow(points), ncol=length(times)))
}

# initial condition 
u0 <- function(points){
 return(sin( 2.* pi * points[,1]) * sin(2.*pi* points[,2]))
}
  
## create pde
pde <- Pde(L, f)
pde$set_boundary_condition(g)
pde$set_initial_condition(u0)

pde$set_boundary_condition(g(pde$dofs_coordinates(),times))
pde$set_initial_condition(u0(pde$dofs_coordinates()))

pde$solve()

## compute L2 norm of the error
u_ex <- exact_solution(pde$dofs_coordinates(), times)
 
error.L2 <- matrix(0, nrow=length(times), ncol=1)
for( t in 1:length(times)){
  error.L2[t] <- sqrt(sum(pde$mass() %*% (u_ex[,t] - u$coefficients()[,t])^2))
  cat(paste0("L2 error at time ",times[t]," ", error.L2[t], "\n"))
}
# ------------------------------------------------------------------------------
## perform evaluation at single point
point = c(0.2, 0.5)
u$eval(point)
# 
# ## evaluate over a 10x10 grid 
x <- seq(0, 1, length.out = 50)
y <- x
points <- expand.grid(x, y)
dim( u$eval(points) )
max(abs( u$eval(points) - as.matrix(exact_solution(points))))
# 
# ## plot solution 
options(warn=-1)
plot(u) 
plot(u) %>% layout(scene =list(camera=list(eye=list(z=1))))
plot(u, showscale=FALSE) # no color bar :)

# ------------------------------------------------------------------------------

u_exact <- Function(Vh)
u_exact$set_coefficients(exact_solution(Vh$basis()$dofs_coordinates(), times))
plot(u_exact) %>% layout(scene =list(camera=list(eye=list(z=1))))
plot(u_exact, showscale=FALSE) # no color bar :)

plot(u)

# contour
contour(u)
