#!/bin/bash

set -euo pipefail

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
PROJDIR="$(dirname ${DIR})"
ATDIR="${PROJDIR}/data/athaliana"

mkdir -p "${ATDIR}"

echo "Downloading data..."
curl https://zenodo.org/record/1419699/files/athaliana.tar.gz \
	-o "${ATDIR}/athaliana.tar.gz"

tar -zxf "${ATDIR}/athaliana.tar.gz" -C "${ATDIR}"

mv "${ATDIR}/athaliana/genome.p-ctg.fa" "${ATDIR}"
rm -r "${ATDIR}/athaliana"
rm "${ATDIR}/athaliana.tar.gz"