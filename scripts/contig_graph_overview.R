library(tidyverse)
library(igraph)

contig_lengths <- read_tsv(snakemake@input[["reference_index"]],
                           col_names = c("contig", "length"),
                           col_types = "cd---")

plot_contig_graph <- function(e, plot_filename, ...) {
  g <- graph_from_data_frame(e, directed = FALSE)
  layout <- layout_in_circle(g, order = names(sort(degree(g))))

  d <- degree(g)
  scaled_d <- 3 + (d - min(d)) * 8 / (max(d) - min(d))

  scaled_edge_width <- with(e, 1 + (n_edges - min(n_edges)) * 4 / (max(n_edges) - min(n_edges)))

  vertex_palette <- RColorBrewer::brewer.pal(9, "BuGn")

  node_meta <- tibble(contig = names(V(g))) %>%
    left_join(contig_lengths, by = "contig") %>%
    mutate(length_bin = cut(length, length(vertex_palette)),
           colour = vertex_palette[as.integer(length_bin)])

  png(plot_filename, ...)
  plot(g,
       layout = layout,
       vertex.size = scaled_d,
       vertex.color = node_meta$colour,
       vertex.label = NA,
       edge.color = rgb(0, 0, 0, 0.3),
       edge.width = scaled_edge_width)

  legend("bottomleft",
         legend = levels(node_meta$length_bin),
         fill = vertex_palette,
         bty = "n",
         box.lwd = 0,
         title = "Contig length",
         title.adj = 0.1)

  mtext(paste0(length(V(g)), " nodes, ", length(E(g)), " edges"),
        side = 3, cex = 1.5, adj = 0.5)
  dev.off()

  invisible(g)
}

full_edge_df <- read_tsv(snakemake@input[["edgelist"]],
    col_names = c("source", "target", "kmer",
        "source_position", "target_position",
        "source_direction", "target_direction"),
    col_types = "cccddcc") %>%
    group_by(source, target) %>%
    summarise(n_edges = n())

edge_df_kmer_density <- full_edge_df %>%
  inner_join(contig_lengths %>% rename(source_length = length), by = c(source = "contig")) %>%
  inner_join(contig_lengths %>% rename(target_length = length), by = c(target = "contig")) %>%
  group_by(source, target, source_length, target_length) %>%
  mutate(edges_per_kb_of_shortest = 1000 * n_edges / pmin(source_length, target_length))

edge_df_filtered <- edge_df_kmer_density %>%
  filter(edges_per_kb_of_shortest > 1, source != target)

g <- plot_contig_graph(edge_df_filtered, snakemake@output[[1]],
		       width = 12, height = 8, units = "in", res = 300)
