create_mesh <- function(data, fe_order = 1) {
  if(fe_order == 1)  
    return(new(Mesh_2D_ORDER_1, data))
  else(fe_order == 2)
    return(new(Mesh_2D_ORDER_2, data))
}
