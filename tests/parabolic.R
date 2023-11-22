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

mesh$set_deltaT(dT) # NUOVO
class(mesh)
plot(mesh)

# create Functional Space
Vh <- FunctionSpace(mesh)

exact_solution <- function(points,t){
  res <- matrix(0, nrow=nrow(points), ncol=length(times))
  for( t in 1:length(times)){
    res[,t] = sin( 2.* pi * points[,1]) * sin(2.*pi* points[,2]) * exp(-times[t])
  }
  return(res)
}

## define differential operator in its strong formulation
u <- Function(Vh)

## define differential operator in its strong formulation
#L <- dt(f) + (-1)*laplace(f) 

L <- dt(u) - laplace(u)

## forcing term
f <- function(points,times){
  res <- matrix(0, nrow=nrow(points), ncol=length(times))
  for( t in 1:length(times)){
    res[,t] = (8*pi^2 -1) * sin( 2.* pi * points[,1]) * sin(2.*pi* points[,2]) * exp(-times[t])
  }
  return(res) 
}
## Dirichlet BC
dirichletBC <- function(points, times){
  return(matrix(0, nrow=nrow(points), ncol=length(times)))
}

# initial condition 
initialCondition <- function(points){
 return(sin( 2.* pi * points[,1]) * sin(2.*pi* points[,2]))
}
  
## create pde
pde <- Pde(L, f)
pde$set_dirichletBC(dirichletBC)
pde$set_initialCondition(initialCondition)
pde$solve()

## compute L2 norm of the error
u_ex <- exact_solution(pde$get_dofs_coordinates(), times)
 
error.L2 <- matrix(0, nrow=length(times), ncol=1)
for( t in 1:length(times)){
  error.L2[t] <- sqrt(sum(pde$get_mass() %*% (u_ex[,t] - pde$solution()[,t])^2))
  cat(paste0("L2 error at time ",times[t]," ", error.L2[t], "\n"))
}
# ------------------------------------------------------------------------------
## perform evaluation at single point
point = c(0.2, 0.5)
u$eval_at(point)
# 
# ## evaluate over a 10x10 grid 
x <- seq(0, 1, length.out = 50)
y <- x
points <- expand.grid(x, y)
dim( u$eval_at(points) )
max(abs( u$eval_at(points) - as.matrix(exact_solution(points))))
# 
# ## plot solution 
# options(warn=-1)
plot(u) %>% hide_colorbar()

# scene = list(camera = list(eye = list(x = 0, y = 0, z = 1.25)),
#              xaxis = list(autorange="reversed"),
#              yaxis = list(autorange="reversed"))
# 
# plot_tmp <- plot(u) %>% 
#   layout(dragmode=FALSE, scene=scene,aspectmode="cube") %>%
#   hide_colorbar() %>% 
#   config(displayModeBar =FALSE)
# 
# plot_tmp

# contour
contour(u)
