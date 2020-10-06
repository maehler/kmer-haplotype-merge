#!/bin/bash

set -euo pipefail

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
PROJDIR="$(dirname ${DIR})"
ATDIR="${PROJDIR}/data/athaliana"

echo "Checking assembly..."
if [ ! -f "${ATDIR}/genome.p-ctg.fa" ] || \
		! md5sum --quiet --status -c "${ATDIR}/genome.p-ctg.fa.md5" ; then
	echo "  Downloading assembly..."
	curl https://zenodo.org/record/1419699/files/athaliana.tar.gz \
		-o "${ATDIR}/athaliana.tar.gz"

	tar -zxf "${ATDIR}/athaliana.tar.gz" -C "${ATDIR}"

	mv "${ATDIR}/athaliana/genome.p-ctg.fa" "${ATDIR}"
	rm -r "${ATDIR}/athaliana"
	rm "${ATDIR}/athaliana.tar.gz"
else
	echo "  Assembly already downloaded"
fi

echo "Checking Illumina reads..."
if [ ! -f "${ATDIR}/illumina_reads.fastq" ] || \
		! md5sum --quiet --status -c "${ATDIR}/illumina_reads.fastq.md5"; then
	echo "  Downloading read data..."
	fastq-dump --stdout SRR3703081 SRR3703082 SRR3703105 > "${ATDIR}/illumina_reads.fastq"
else
	echo "  Read data already downloaded"
fi
