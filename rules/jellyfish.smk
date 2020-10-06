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
