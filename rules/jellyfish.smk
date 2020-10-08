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

rule executables:
    output: 'bin/identify_homozygous_kmers'
    conda: '../envs/jellyfish.yaml'
    shell:
        """
        export PKG_CONFIG_PATH=$CONDA_PREFIX/lib/pkgconfig
        cd src/identify_homozygous_kmers
        make install
        """
