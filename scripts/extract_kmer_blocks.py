from collections import defaultdict

def process_pairs(source, target, alignments, k, max_distance=None):
    sorted_alignments = sorted(alignments, key=lambda x: x['source_position'])
    target_positions = [x['target_position'] for x in sorted_alignments]

    prev_increasing = None
    prev_same_strand = None
    prev_target_end = None
    block = 1

    blocks = defaultdict(dict)

    for i, (x, y) in enumerate(zip(target_positions[:-1], target_positions[1:])):
        increasing = y - x > 0
        if prev_increasing is None:
            prev_increasing = increasing

        source_position = sorted_alignments[i]['source_position']
        source_direction = sorted_alignments[i]['source_direction']
        target_position = sorted_alignments[i]['target_position']
        target_direction = sorted_alignments[i]['target_direction']

        source_start = source_position if source_direction == '+' else source_position - k
        source_end = source_position + k if source_direction == '+' else source_position
        target_start = target_position if target_direction == '+' else target_position - k
        target_end = target_position + k if target_direction == '+' else target_position

        same_strand = source_direction == target_direction
        if prev_same_strand is None:
            prev_same_strand = same_strand

        if prev_target_end is None:
            prev_target_end = target_end

        if increasing != prev_increasing or \
                same_strand != prev_same_strand or \
                target_start - prev_target_end > max_distance:
            block += 1

        blocks[block]['source_start'] = min(blocks[block].get('source_start', source_start), source_start)
        blocks[block]['source_end'] = max(blocks[block].get('source_end', source_end), source_end)
        blocks[block]['target_start'] = min(blocks[block].get('target_start', target_start), target_start)
        blocks[block]['target_end'] = max(blocks[block].get('target_end', target_end), target_end)
        blocks[block]['same_direction'] = increasing

        prev_increasing = increasing
        prev_same_strand = same_strand
        prev_target_end = target_end

    return blocks

def extract_candidates(graphfile, k, max_distance):
    current_pair = None
    pair_alignments = []

    kmer_blocks = {}

    with open(graphfile) as f:
        for line in f:
            contig_pair = line.strip().split(' ')[:2]
            if contig_pair != current_pair:
                if current_pair is not None:
                    kmer_blocks[tuple(current_pair)] = process_pairs(current_pair[0], current_pair[1], pair_alignments, k, max_distance)
                current_pair = contig_pair
                pair_alignments = []
            metadata = line.strip().split(' ')[2:]
            pair_alignments.append(
                {
                    'name': metadata[0],
                    'source_position': int(metadata[1]),
                    'target_position': int(metadata[2]),
                    'source_direction': metadata[3],
                    'target_direction': metadata[4]
                }
            )

    return kmer_blocks

def main():
    graphfile = snakemake.input[0]
    output = snakemake.output[0]
    k = int(snakemake.wildcards['k'])
    max_distance = int(snakemake.wildcards['max_distance'])

    kmer_blocks = extract_candidates(graphfile, k, max_distance)

    with open(output, 'w') as f:
        print('\t'.join(('source', 'target', 'source_start', 'source_end',
            'target_start', 'target_end', 'same_direction')),
            file=f)
        for (source, target), blocks in kmer_blocks.items():
            for i, b in blocks.items():
                print('{source}\t{target}\t{source_start}\t{source_end}\t{target_start}\t{target_end}\t{same_direction}' \
                    .format(source=source, target=target, **b), file=f)

if __name__ == '__main__':
    main()
