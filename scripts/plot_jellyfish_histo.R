library(tidyverse)

kmer_counts <- read_delim(snakemake@input[[1]],
        delim = " ",
        col_names = c("freq", "count"),
        col_types = "dd")  %>%
    mutate(direction = sign(diff(c(count[1], count)))) %>%
    mutate(peak = c(diff(direction) < 0, FALSE))

round_up_nice <- function(x, nice=c(1,2,4,5,6,8,10)) {
  if(length(x) != 1) stop("'x' must be of length 1")
  10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
}

# Find a good coverage cutoff, ignoring any peaks at 0 or 1.
max_index <- kmer_counts %>%
    filter(freq > 3) %>%
    with(which.max(count) + sum(kmer_counts$freq <= 3))

if (max_index == 4) {
    # If there is no second peak, show everything.
    max_index <- 1
}
n_bins <- nrow(kmer_counts)
cutoff_freq <- with(kmer_counts, round_up_nice(freq[max_index:n_bins][count[max_index:n_bins] < 0.001 * max(count)][1]))

ggplot(kmer_counts, aes(freq, count)) +
    geom_line() +
    coord_cartesian(xlim = c(0, cutoff_freq), y = c(0, kmer_counts$count[max_index])) +
    labs(x = "k-mer depth", y = "Count")

ggsave(snakemake@output[["hist_plot"]])

kmer_counts %>%
    filter(peak, freq > 1) %>%
    distinct(freq) %>%
    top_n(2, -freq) %>%
    with(round(freq[2] * c(0.9, 1.1))) %>%
    write_lines(path = snakemake@output[["hom_limits"]])
