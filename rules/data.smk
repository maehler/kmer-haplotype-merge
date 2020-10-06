rule download_data:
    output:
        asm='data/athaliana/genome.p-ctg.fa',
        reads='data/athaliana/illumina_reads.fastq'
    script: '../scripts/init.sh'
