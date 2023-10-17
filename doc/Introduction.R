## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, fig.width= 7., 
                                   fig.height= 4., 
                                   fig.align= "center",
                                   collapse = TRUE,
                                   comment = "#>")

library(femR)
options(warn = -1)

## ----echo=FALSE---------------------------------------------------------------
exact_solution <- function(points){
    return( sin(2.*pi*points[,1])*sin(2.*pi*points[,2]))
}
data("unit_square", package = "femR")
domain <- Mesh(unit_square)
eval_points <- unit_square$nodes
data_plot <- data.frame(X=eval_points[,1], Y=eval_points[,2],
                        Z=exact_solution(eval_points), coeff=exact_solution(eval_points))
I=(domain$elements()[,1]); J=(domain$elements()[,2]); K=(domain$elements()[,3])  
fig_exact <- plot_ly(data_plot, type="contour", x=~X, y=~Y, z=~Z, intensity=~Z,color = ~Z,
                 i = I, j = J, k = K,
          contours=list(showlabels = TRUE),
          colorbar=list(title="")) %>%
   layout(xaxis = list(title = "", showgrid=F, zeroline=F, ticks="", showticklabels=F),
          yaxis = list(title = "", showgrid=F, zeroline=F, ticks="", showticklabels=F)) %>% layout(scene=list(aspectmode="cube"))
fig2 <- plot(domain) %>% layout(scene=list(aspectmode="cube"))
subplot(fig_exact %>% hide_colorbar(), fig2, margin=0.05) %>%
        layout(annotations = list(
            list(x = 0.155 , y = 1.05, text = "exact solution", showarrow = F, xref='paper', yref='paper'),
            list(x = 0.8 , y = 1.05, text = "mesh", showarrow = F, xref='paper', yref='paper'))
)

## ----code block---------------------------------------------------------------
## 1. Defining domain
data("unit_square", package="femR")
mesh <- Mesh(unit_square)
## 2. Defining the solution of the PDE
Vh <- FunctionSpace(mesh, fe_order=2) 
f <- f <- Function(Vh)
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

## -----------------------------------------------------------------------------
fig_h <- contour(f)
subplot(fig_exact %>%hide_colorbar(), fig_h, margin = 0.05) %>%layout(annotations = list(
            list(x=0.155, y=1.05, text = "exact solution",    showarrow = F, 
                 xref='paper', yref='paper'),
            list(x=0.875,  y=1.05, text = "discrete solution", showarrow = F, 
                 xref='paper', yref='paper')))

## ----Creating mesh, fig.width=4, fig.height=4---------------------------------
data("unit_square", package="femR")
unit_square <- Mesh(unit_square)
plot(unit_square)

## ----FreeFem++, eval=FALSE----------------------------------------------------
#  filename <- "path/to/domain.mesh" #domain.mesh stored by the FreeFem++ utility savemesh()
#  domain <- read_mesh(filename)
#  mesh <- Mesh(domain)

## ----Defining Differential Operator-------------------------------------------
## solution of the PDE
Vh <- FunctionSpace(unit_square, fe_order=2)
f <- Function(Vh)

## ----Differentail Operator----------------------------------------------------
L <- -laplace(f) # L <- -div(I*grad(f))
                 # In general:
                 # L <- (-mu)*laplace(f) + dot(b,grad(f)) + c*f

## -----------------------------------------------------------------------------
## forcing term
u <- function(points){
    return(8.*pi^2* sin(2.*pi*points[,1])*sin(2.*pi*points[,2])) 
}

## Dirichlet BC
dirichletBC <- function(points){
  return(matrix(0,nrow=nrow(points), ncol=1))
}

## ----Defining PDE object------------------------------------------------------
## create pde
pde <- Pde(L, u, dirichletBC)

## ----Solving the problem------------------------------------------------------
## solve problem
pde$solve()

## ----Compute error------------------------------------------------------------
exact_solution <- function(points){
    return( sin(2.*pi*points[,1])*sin(2.*pi*points[,2]))
}
u_ex <- as.matrix(exact_solution(pde$get_dofs_coordinates()))
error.L2 <- sqrt(sum(pde$get_mass()%*%(u_ex-pde$solution())^2))
cat("L2 error = ", error.L2, "\n")

## ----Plots--------------------------------------------------------------------
## plot solution 
plot(f)

## ----Contour, fig.width=5-----------------------------------------------------
# contour 
contour(f)

## ----Using Plotly-------------------------------------------------------------
plot(f, colorscale="Jet") %>% layout(scene=list(aspectmode="cube"))

