rule jellyfish_count:
    input: lambda wildcards: config[wildcards.dataset][wildcards.data_type]
    output: 'counts/{dataset}/{data_type}_{k}mer-counts.jf'
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
    input: 'counts/{dataset}/{data_type}_{k}mer-counts.jf'
    output: 'results/{dataset}/{data_type}_{k}mer-hist.txt'
    conda: '../envs/jellyfish.yaml'
    threads: 4
    shell:
        """
        jellyfish histo \\
            --output {output} \\
            {input}
        """

rule plot_jellyfish_histo:
    input: 'results/{dataset}/{data_type}_{k}mer-hist.txt'
    output:
        hist_plot='results/{dataset}/{data_type}_{k}mer-hist.png',
        hom_limits='results/{dataset}/{data_type}_{k}mer-hom-limits.txt'
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
        read_jf='counts/{dataset}/short_reads_{k}mer-counts.jf',
        assembly_jf='counts/{dataset}/assembly_{k}mer-counts.jf',
        hom_limits='results/{dataset}/short_reads_{k}mer-hom-limits.txt'
    output: 'results/{dataset}/short_reads_assembly_{k}mer_homozygous.fasta'
    shell:
        """
        {input.executable} \\
            $(cat {input.hom_limits}) \\
            {input.read_jf} \\
            {input.assembly_jf} > {output}
        """
