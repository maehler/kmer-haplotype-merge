rule contig_graph:
    input: 'alignments/{dataset}/short_reads_assembly_{k}mer_alignments.bam'
    output: temp('results/{dataset}/short_reads_assembly_{k}mer_contig_graph.edgelist')
    conda: '../envs/merging.yaml'
    script: '../scripts/build_contig_graph.py'

rule sort_contig_graph:
    input: 'results/{dataset}/short_reads_assembly_{k}mer_contig_graph.edgelist'
    output: 'results/{dataset}/short_reads_assembly_{k}mer_contig_graph_sorted.edgelist'
    shell: 'sort -k1,2 {input} > {output}'

rule contig_graph_overview:
    input:
        edgelist='results/{dataset}/short_reads_assembly_{k}mer_contig_graph_sorted.edgelist',
        reference_index=lambda wildcards: '{0}.fai'.format(config[wildcards.dataset]['assembly'])
    output: 'results/{dataset}/short_reads_assembly_{k}mer_contig_graph_overview.png'
    conda: '../envs/r_graphs.yaml'
    script: '../scripts/contig_graph_overview.R'

rule extract_kmer_blocks:
    input: 'results/{dataset}/short_reads_assembly_{k}mer_contig_graph_sorted.edgelist'
    output: 'results/{dataset}/short_reads_assembly_{k}mer_blocks_maxdiff{max_difference}.tsv'
    conda: '../envs/merging.yaml'
    script: '../scripts/extract_kmer_blocks.py'

rule plot_linked_contigs:
    input:
        kmer_blocks='results/{dataset}/short_reads_assembly_{k}mer_blocks_maxdiff{max_difference}.tsv',
        reference_index=lambda wildcards: '{0}.fai'.format(config[wildcards.dataset]['assembly'])
    output: 'results/{dataset}/short_reads_assembly_{k}mer_maxdiff{max_difference}_plots/{contig1}_{contig2}.png'
    conda: '../envs/r_graphs.yaml'
    script: '../scripts/plot_linked_contigs.R'
