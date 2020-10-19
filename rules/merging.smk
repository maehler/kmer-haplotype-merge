rule contig_graph:
    input: 'alignments/{dataset}/short_reads_assembly_{k}mer_alignments.bam'
    output: 'results/{dataset}/short_reads_assembly_{k}mer_contig_graph.edgelist',
    conda: '../envs/merging.yaml'
    script: '../scripts/build_contig_graph.py'

rule contig_graph_overview:
    input: 'results/{dataset}/short_reads_assembly_{k}mer_contig_graph.edgelist'
    output: 'results/{dataset}/short_reads_assembly_{k}mer_contig_graph_overview.png'
    conda: '../envs/r_graphs.yaml'
    script: '../scripts/contig_graph_overview.R'

rule extract_kmer_blocks:
    input: 'results/{dataset}/short_reads_assembly_{k}mer_contig_graph.edgelist'
    output: 'results/{dataset}/short_reads_assembly_{k}mer_blocks_maxdist{max_distance}.tsv'
    conda: '../envs/merging.yaml'
    script: '../scripts/extract_kmer_blocks.py'
