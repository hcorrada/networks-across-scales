library(tidygraph)
library(ggraph)

graph <- highschool %>% 
  as_tbl_graph() %>%
  activate(edges) %>%
  filter(year == 1958) %>%
  convert(to_undirected) %>%
  convert(to_simple)

graph %>%
  ggraph(layout="kk") +
  geom_edge_fan() +
  geom_node_text(aes(label=name)) +
  theme_graph(foreground=NA)

plot(graph)

## FW iterative version
fw_dists_iter <- function(amat) {
  nvertices <- nrow(amat)
  dmat <- matrix(Inf, nvertices, nvertices)
  diag(dmat) <- 0
  dmat[amat == 1] <- 1
  
  for (k in seq_len(nvertices)) {
    for (i in seq_len(nvertices)) {
      for (j in seq_len(nvertices)) {
        dk <- dmat[i,k] + dmat[k,j]
        if (dk < d[i,j]) {
          d[i,j] <- dk
        }
      }
    }
  }
  dmat
}

## FW using vectorized operations
fw_dists <- function(amat) {
  nvertices <- nrow(amat)
  dmat <- matrix(Inf, nvertices, nvertices)
  diag(dmat) <- 0
  dmat[amat == 1] <- 1
  
  for (k in seq_len(nvertices)) {
    #    cat(".")
    dk <- outer(dmat[,k],dmat[k,],FUN="+")
    switch <- dmat > dk
    dmat[switch] <- dk[switch]
  }
  #  cat("\n")
  dmat
}

