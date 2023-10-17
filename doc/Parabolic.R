## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, fig.width= 6., 
                                   fig.height= 5., 
                                   fig.align= "center",
                                   collapse = TRUE,
                                   comment = "#>")

library(femR)
options(warn = -1)

## ----echo=FALSE---------------------------------------------------------------
exact_solution <- function(points,t){
  res <- matrix(0, nrow=nrow(points), ncol=length(times))
  for( t in 1:length(times)){
    res[,t] = sin( 2.* pi * points[,1]) * sin(2.*pi* points[,2]) * exp(-times[t])
  }
  return(res)
}

t0 = 0.
t_max = 1
dT = 1e-2
M = (t_max - t0)/dT + 1
times = seq(t0,t_max,length=M)

data("unit_square", package = "femR")
domain <- Mesh(unit_square) %X% times
eval_points <- unit_square$nodes
data_plot <- data.frame(X=rep(eval_points[,1], times=length(times)), Y=rep(eval_points[,2], times=length(times)),
                        coeff=as.vector(exact_solution(eval_points, times)),
                        t=rep(times, times=nrow(eval_points)))
I=(domain$elements()[,1]); J=(domain$elements()[,2]); K=(domain$elements()[,3]) 
limits = c(min(data_plot$coeff), max(data_plot$coeff))

fig_exact <- plot_ly(data_plot, type="contour", x=~X, y=~Y, z=~coeff, intensity=~coeff,color = ~coeff,
                     frame =~t,
                 i = I, j = J, k = K, zmin=limits[1], zmax=limits[2],
          contours=list(showlabels = TRUE),
          colorbar=list(title="")) %>%
   layout(xaxis = list(title = "", showgrid=F, zeroline=F, ticks="", showticklabels=F),
          yaxis = list(title = "", showgrid=F, zeroline=F, ticks="", showticklabels=F)) %>%
      animation_opts(frame=5) %>% 
      animation_slider(currentvalue = list(prefix ="t = "))

fig_exact
# fig2 <- plot(domain) %>% layout(scene=list(aspectmode="cube"))
# subplot(fig_exact %>% hide_colorbar(), fig2, margin=0.05) %>%
#         layout(annotations = list(
#             list(x = 0.155 , y = 1.05, text = "exact solution", showarrow = F, xref='paper', yref='paper'),
#             list(x = 0.8 , y = 1.05, text = "mesh", showarrow = F, xref='paper', yref='paper'))
# )

## ----code block---------------------------------------------------------------
## load domain data
data("unit_square", package="femR")

## time steps
t0 = 0.
t_max = 1
dT = 1e-2
M = (t_max - t0)/dT + 1
times = seq(t0,t_max,length=M)

## Spatio-temporal domain
domain = Mesh(unit_square) %X% times

## create Functional Space
Vh <- FunctionSpace(domain)

## define differential operator in its strong formulation
f <- Function(Vh)

## define differential operator in its strong formulation
L <- dt(f) - laplace(f)

## forcing term
u <- function(points,times){
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

## initial condition 
initialCondition <- function(points){
 return(sin( 2.* pi * points[,1]) * sin(2.*pi* points[,2]))
}
  
## create pde
pde <- Pde(L, u, dirichletBC, initialCondition)

## solve problem
pde$solve()

## plot solution
contour(f)


