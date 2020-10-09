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
            samtools view -b | \\
            samtools sort -n -o {output}
        """
