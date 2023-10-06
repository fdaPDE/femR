## Mesh2D Class Definition
.Mesh2DCtr <- setRefClass(
  Class = "Mesh2DObject",
  fields = c(
    data  = "ANY"
  ),
  methods = c(
    nodes = function(){
      data$nodes()
    },
    elements = function(){
      data$elements()
    },
    boundary = function(){
      data$boundary()
    },
    neighbors = function(){
      data$neighbors()
    }
  )
)

setGeneric("Mesh2D", function(domain) standardGeneric("Mesh2D"))
setMethod("Mesh2D", signature = c(domain="list"),
          function(domain){
            domain$elements <- domain$elements - 1
            storage.mode(domain$elements) <- "integer" 
            .Mesh2DCtr(data=new(Mesh_2D, domain))                                                  
})

## Mesh - auxiliary methods
unroll_edges_aux <- function(Mesh2D){
  mesh <- Mesh2D$data
  edges <- matrix(nrow=3*nrow(mesh$elements()), ncol=2)
  for(i in 1:nrow(mesh$elements())){
    edges[(3*(i-1) + 1),]   = mesh$elements()[i,c(1,2)] + 1
    edges[(3*(i-1) + 2),]   = mesh$elements()[i,c(2,3)] + 1
    edges[(3*(i-1) + 3),]   = mesh$elements()[i,c(3,1)] + 1
  }
  edges
}

setGeneric("unroll_edges", function(Mesh2D) standardGeneric("unroll_edges"))
setMethod("unroll_edges", "Mesh2DObject", function(Mesh2D){
  unroll_edges_aux(Mesh2D)
})

plot_mesh_aux <- function(x, ...){
  edges <- unroll_edges(x)
  plot_ly(...) %>% 
    add_markers(x = x$data$nodes()[,1],
                y = x$data$nodes()[,2],
                color = I('black'), size = I(1),
                hoverinfo = 'text',
                text = paste('</br><b> Coordinates:', round(x$data$nodes()[,1],2),
                             round(x$data$nodes()[,2],2)),
                showlegend = T,
                visible = T) %>%
    add_segments(x = x$data$nodes()[edges[,1],1],
                 y = x$data$nodes()[edges[,1],2],
                 xend = x$data$nodes()[edges[,2],1],
                 yend = x$data$nodes()[edges[,2],2], 
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
}

setMethod("plot", "Mesh2DObject", function(x, ...){
  plot_mesh_aux(x, ...)  
})

