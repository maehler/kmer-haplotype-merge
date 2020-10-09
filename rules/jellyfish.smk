rule jellyfish_count:
    input: 'data/{species}/{prefix}.{ext}'
    output: 'counts/{species}/{prefix}.{ext}_{k}mer-counts.jf'
    conda: '../envs/jellyfish.yaml'
    threads: 4
    shell:
        """
        jellyfish count \\
            --canonical \\
            -m {wildcards.k} \\
            -s 6G \\
            -t {threads} \\
            -o {output} \\
            {input}
        """

rule jellyfish_histo:
    input: 'counts/{species}/{prefix}.{ext}_{k}mer-counts.jf'
    output: 'results/{species}/{prefix}.{ext}_{k}mer-hist.txt'
    conda: '../envs/jellyfish.yaml'
    threads: 4
    shell:
        """
        jellyfish histo \\
            --output {output} \\
            {input}
        """

rule plot_jellyfish_histo:
    input: 'results/{species}/{prefix}.{ext}_{k}mer-hist.txt'
    output:
        hist_plot='results/{species}/{prefix}.{ext}_{k}mer-hist.png',
        hom_limits='results/{species}/{prefix}.{ext}_{k}mer-hom-limits.txt'
    script: '../scripts/plot_jellyfish_histo.R'

rule executables:
    output: 'bin/identify_homozygous_kmers'
    conda: '../envs/jellyfish.yaml'
    shell:
        """
        export PKG_CONFIG_PATH=$CONDA_PREFIX/lib/pkgconfig
        cd src/identify_homozygous_kmers
        make install
        """

rule homozygous_kmers:
    input:
        executable='bin/identify_homozygous_kmers',
        read_jf='counts/{species}/{read_prefix}.{read_ext}_{k}mer-counts.jf',
        assembly_jf='counts/{species}/{asm_prefix}.{asm_ext}_{k}mer-counts.jf',
        hom_limits='results/{species}/{read_prefix}.{read_ext}_{k}mer-hom-limits.txt'
    output: 'results/{species}/{read_prefix}.{read_ext}_{asm_prefix}.{asm_ext}_{k}mer_homozygous.fasta'
    shell:
        """
        {input.executable} $(cat {input.hom_limits}) {input.read_jf} {input.assembly_jf} > {output}
        """

rule index_bowtie:
    input:
        fasta='data/{species}/{prefix}.{ext}'
    output: multiext('alignments/reference/{species}_{prefix}.{ext}', \
        '.1.ebwt', '.2.ebwt', '.rev.1.ebwt', '.rev.2.ebwt')
    conda: '../envs/bowtie.yaml'
    threads: 4
    shell:
        """
        INDEX={output[0]}
        INDEX=${{INDEX%%.1.ebwt}}
        bowtie-build \\
            --threads {threads} \\
            {input.fasta} \\
            ${{INDEX}}
        """

rule align_kmers:
    input:
        kmers='results/{species}/{read_prefix}.{read_ext}_{asm_prefix}.{asm_ext}_{k}mer_homozygous.fasta',
        reference=multiext('alignments/reference/{species}_{asm_prefix}.{asm_ext}',
            '.1.ebwt', '.2.ebwt', '.rev.1.ebwt', '.rev.2.ebwt')
    output: 'alignments/{species}/{read_prefix}.{read_ext}_{asm_prefix}.{asm_ext}_{k}mer_alignments.bam'
    conda: '../envs/bowtie.yaml'
    threads: 4
    shell:
        """
        INDEX={input.reference[0]}
        INDEX=${{INDEX%%.1.ebwt}}
        bowtie \\
            -f \\
            -v 0 \\
            --seedlen {wildcards.k} \\
            -k 20 \\
            --sam \\
            --threads {threads} \\
            ${{INDEX}} \\
            {input.kmers} | \\
            samtools view -bo {output}
        """
