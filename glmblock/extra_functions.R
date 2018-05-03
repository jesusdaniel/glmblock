uptri <- function(network, diag = F)  {
  network[upper.tri(network, diag = diag)]
}
to_matrix <- function(vec) {
  Nodes = (sqrt(8*length(vec) +1) + 1)/2
  network = matrix(0,nrow = Nodes, ncol = Nodes)
  network[upper.tri(network)] = as.double(vec)
  network = network + t(network)
  return(network)
}

