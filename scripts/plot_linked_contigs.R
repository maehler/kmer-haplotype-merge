library(tidyverse)

contig1 <- snakemake@wildcards[["contig1"]]
contig2 <- snakemake@wildcards[["contig2"]]

contig_lengths <- read_tsv(snakemake@input[["reference_index"]],
                           col_names = c("contig", "length"),
                           col_types = "cd---")

kmer_blocks <- read_tsv(snakemake@input[["kmer_blocks"]]) %>% 
  inner_join(contig_lengths %>% rename(source_length = length), by = c(source = "contig")) %>% 
  inner_join(contig_lengths %>% rename(target_length = length), by = c(target = "contig")) %>% 
  mutate(source_block_length = source_end - source_start,
         target_block_length = target_end - target_start)

subset_df <- kmer_blocks %>%
    filter(source %in% c(contig1, contig2),
           target %in% c(contig1, contig2)) %>%
  mutate(poly_end = NA) %>% 
  filter(source_block_length > 1000)

x_coords <- subset_df %>%
  rename(tmp_target_start = target_start,
         tmp_target_end = target_end) %>% 
  mutate(target_start = ifelse(same_direction, tmp_target_start, tmp_target_end),
         target_end = ifelse(same_direction, tmp_target_end, tmp_target_start)) %>%
  select(source_start, source_end, target_start, target_end, poly_end) %>%
  as.matrix %>% t %>% as.vector

png(snakemake@output[[1]], height = 8, width = 10, units = "in", res = 300)
plot.new()
plot.window(xlim = c(0, max(subset_df$source_length, subset_df$target_length)),
            ylim = c(0, 1))

polygon(x_coords, rep(c(0, 0, 1, 1, NA), length(x_coords) / 5),
        col = c("gold", "forestgreen"), border = NA)
lines(x = c(0, subset_df$source_length[1]), y = rep(0, 2), lwd = 3, col = "steelblue")
lines(x = c(0, subset_df$target_length[1]), y = rep(1, 2), lwd = 3, col = "steelblue")
mtext(paste0(contig1, "\n", format(subset_df$source_length[1], big.mark = ","), " basepairs"),
      side = 1, line = 0.5, at = subset_df$source_length[1]/2)
mtext(paste0(contig2, "\n", format(subset_df$target_length[1], big.mark = ","), " basepairs"),
      side = 3, line = -0.5, at = subset_df$target_length[1]/2)
dev.off()
