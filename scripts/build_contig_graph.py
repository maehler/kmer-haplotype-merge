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

def walk_bamfile(fname):
    bam = pysam.AlignmentFile(fname, 'rb', require_index=False)

    current_kmer = None
    current_kmer_aln = []

    for aln in bam:
        if current_kmer != aln.query_name:
            if current_kmer is not None and len(current_kmer_aln) != 2:
                print(f'error: found {len(current_kmer_aln)} alignments '
                      f'for {current_kmer}, should be 2', file=sys.stderr)
                sys.exit(1)
            if current_kmer is not None:
                yield construct_edge(current_kmer_aln)
            current_kmer = aln.query_name
            current_kmer_aln = []
        current_kmer_aln.append(aln)

def main():
    bamfile = snakemake.input[0]
    graphfile = snakemake.output[0]

    edges = walk_bamfile(bamfile)
    with open(graphfile, 'w') as f:
        for source, target, attr in edges:
            print('\t'.join(map(str, (source, target, attr['label'],
                             attr['source_position'],
                             attr['target_position'],
                             attr['source_direction'],
                             attr['target_direction']))), file=f)

if __name__ == '__main__':
    main()
