import networkx as nx
import pysam
import sys

def construct_edge(aln):
    # Set the lexicographically smaller contig name as source.
    if aln[0].reference_name <= aln[1].reference_name:
        source_index = 0
        target_index = 1
    else:
        source_index = 1
        target_index = 0

    return (
        aln[source_index].reference_name,
        aln[target_index].reference_name,
        {
            'label': aln[source_index].query_name,
            'source_position': aln[source_index].reference_start,
            'target_position': aln[target_index].reference_start,
            'source_direction': '+' if aln[source_index].flag & 16 == 0 else '-',
            'target_direction': '+' if aln[target_index].flag & 16 == 0 else '-'
        }
    )

def walk_bamfile(fname, batch_size=1000):
    bam = pysam.AlignmentFile(fname, 'rb', require_index=False)

    g = nx.MultiGraph()
    edge_batch = []

    current_kmer = None
    current_kmer_aln = []

    for aln in bam:
        if len(edge_batch) >= batch_size:
            g.add_edges_from(edge_batch)
            edge_batch = []
        if current_kmer != aln.query_name:
            if current_kmer is not None and len(current_kmer_aln) != 2:
                print(f'error: found {len(current_kmer_aln)} alignments '
                      f'for {current_kmer}, should be 2', file=sys.stderr)
                sys.exit(1)
            if current_kmer is not None:
                edge_batch.append(construct_edge(current_kmer_aln))
            current_kmer = aln.query_name
            current_kmer_aln = []
        current_kmer_aln.append(aln)

    # Add the last batch of edges
    g.add_edges_from(edge_batch)

    return g

def main():
    bamfile = snakemake.input[0]
    graphfile = snakemake.output[0]

    contig_graph = walk_bamfile(bamfile)
    # nx.write_graphml(contig_graph, graphfile)
    nx.write_edgelist(contig_graph, graphfile,
        data=['label', 'source_position', 'target_position',
            'source_direction', 'target_direction'])
    print(f'{contig_graph.number_of_nodes()} nodes, {contig_graph.number_of_edges()} edges')

if __name__ == '__main__':
    main()
