## Mesh - auxiliary methods
setGeneric("unroll_edges", function(mesh) 0)
setMethod("unroll_edges", "Rcpp_Mesh_2D", function(mesh){
  edges <- matrix(nrow=3*nrow(mesh$elements()), ncol=2)
  for(i in 1:nrow(mesh$elements())){
    edges[(3*(i-1) + 1),]   = mesh$elements()[i,c(1,2)] + 1
    edges[(3*(i-1) + 2),] = mesh$elements()[i,c(2,3)] + 1
    edges[(3*(i-1) + 3),] = mesh$elements()[i,c(3,1)] + 1
  }
  edges
  }
)

setMethod("plot", "Rcpp_Mesh_2D", function(x, ...){
  edges <- unroll_edges(x)
  plot_ly(...) %>% 
    add_markers(x = x$nodes()[,1],
                y = x$nodes()[,2],
                color = I('black'), size = I(1),
                hoverinfo = 'text',
                text = paste('</br><b> Coordinates:', round(x$nodes()[,1],2),
                             round(x$nodes()[,2],2)),
                showlegend = T,
                visible = T) %>%
    add_segments(x = x$nodes()[edges[,1],1],
                 y = x$nodes()[edges[,1],2],
                 xend = x$nodes()[edges[,2],1],
                 yend = x$nodes()[edges[,2],2], 
                 color = I('black'), size = I(1),
                 showlegend = F) %>%
    layout(
      xaxis = list(
        title = '',
        showgrid = F,
        zeroline = F,
        showticklabels = F
      ),
      yaxis = list(
        title = '',
        showgrid = F,
        zeroline = F,
        showticklabels = F
      ))
})
