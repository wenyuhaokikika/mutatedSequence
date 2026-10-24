#!/usr/bin/env bash
###
# @Description: download Ensembl reference files used by mutatedSequence.
# Usage: ./downLoad.sh 109
#
# For a human genome, the pipeline can read the gzipped FASTA directly but
# will load it into memory.  For low-memory random access, decompress the
# genome and index it:
#   gunzip -c resource/Homo_sapiens.GRCh38.dna.toplevel.fa.gz \
#       > resource/Homo_sapiens.GRCh38.dna.toplevel.fa
#   samtools faidx resource/Homo_sapiens.GRCh38.dna.toplevel.fa
###
set -euo pipefail
filepath=$(cd "$(dirname "$0")" && pwd)
source "$filepath/log.sh"

release=${1:-109}
base=${ENSEMBL_FTP:-https://ftp.ensembl.org/pub}/release-"$release"
outdir=${2:-"$filepath/resource"}
mkdir -p "$outdir"

download() {
  local url=$1
  local dest=$2
  log info "download $url"
  wget -c "$url" -O "$dest"
  if [ $? -ne 0 ]; then
    log error "download failed: $url"
    exit 1
  fi
}

download "$base/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz" \
  "$outdir/Homo_sapiens.GRCh38.cdna.all.fa.gz"
download "$base/fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz" \
  "$outdir/Homo_sapiens.GRCh38.cds.all.fa.gz"
download "$base/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz" \
  "$outdir/Homo_sapiens.GRCh38.pep.all.fa.gz"
download "$base/gtf/homo_sapiens/Homo_sapiens.GRCh38.${release}.chr_patch_hapl_scaff.gtf.gz" \
  "$outdir/Homo_sapiens.GRCh38.${release}.chr_patch_hapl_scaff.gtf.gz"
download "$base/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz" \
  "$outdir/Homo_sapiens.GRCh38.dna.toplevel.fa.gz"

log info "downloaded files are in $outdir"
log info "tip: gunzip the genome and run 'samtools faidx' to avoid loading 3 GB into memory"
