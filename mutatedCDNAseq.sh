#!/bin/bash
# DEPRECATED: the bcftools consensus / per-exon CDS pipeline was replaced by
# the Python implementation.  It mishandles negative-strand alleles, variants
# spanning exon junctions, multiallelic loci and symbolic SVs.
#
# Use instead:
#   python -m mutatedseq run -g genome.fa -a annotation.gtf -v sample.vcf -o out \
#       --cds Homo_sapiens.GRCh38.cds.all.fa.gz \
#       --cdna Homo_sapiens.GRCh38.cdna.all.fa.gz \
#       --pep Homo_sapiens.GRCh38.pep.all.fa.gz
#
# The original script is kept in git history and below as a reference.
echo "DEPRECATED: use 'python -m mutatedseq run ...' (see README.md)." >&2
exit 1

#从bed文件获取fasta文件
filepath=$(cd "$(dirname "$0")"; pwd) # 获取当前文件夹的绝对路径
cat $filepath/loger.txt;
check_status () {
    if [ $? -ne 0 ]; then
        log error "bcftools consensus ERROR."
        exit 1
    else
        log info "bcftools consensus successfully."
    fi
}
while getopts m:d:o: flag
do
    case "${flag}" in
        m) mutation=${OPTARG};;
        d) database=${OPTARG};;
        o) out=${OPTARG};;
    esac
done
source $filepath/log.sh;
log info "mutation: $mutation, database: $database, out: $out"
if [ ! -d $out ]; then
  mkdir $out
fi
if [ ! -f $database/cds.fa ]
  then
    log info 'cds is bed file, convert it to cds.fa';
    awk -F '\t' '{print ">"$1":"$2"-"$3" transcript_id:"$5";rank:"$4";strand:"$6";stopWithCodon:"$7";addHead:"$8"\n"$9}' $database/cds.bed  | seqkit seq -w 60 > $database/cds.fa;
  else
    log info 'cds is fasta file';
fi
log info 'extract mutations from vcf file'
if [ ! -f $mutation.tbi ];then tabix $mutation;fi
bcftools view -R $database/cds.bed $mutation | bcftools view -f 'PASS,.' | bgzip -c > $out/all.vcf.gz && tabix $out/all.vcf.gz
log info 'get mutated dna sequence'
bcftools consensus -f $database/cds.fa $out/all.vcf.gz -o $out/mutated.fa > $out/consensusLog.txt 2>&1
check_status
log info 'get mutated dna sequence'
seqkit fx2tab $out/mutated.fa > $out/mutated.bed
