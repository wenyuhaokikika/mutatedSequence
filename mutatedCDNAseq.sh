#从bed文件获取fasta文件
#!/bin/bash
while getopts m:d:o: flag
do
    case "${flag}" in
        m) mutation=${OPTARG};;
        d) database=${OPTARG};;
        o) out=${OPTARG};;
    esac
done
source ./log.sh;
log info "mutation: $mutation, database: $database, out: $out"
# 如果cds输入未bed文件将他转化为fa文件，如果是fa文件则不改变。
#if [[ $cdsfa =~ 'bed' ]]
if[ ! -f $database/cds.fa ];
  then
    log info 'cds is bed file, convert it to cds.fa'
    awk -F '\t' '{print ">"$1":"$2"-"$3" transcript_id:"$5";rank:"$4";strand:"$6";stopWithCodon:"$7";addHead:"$8"\n"$8}' $database/cds.bed  | seqtk seq -l 60 > $database/cds.fa
  else
    log info 'cds is fasta file'
fi

# 只选择在bed范围内的突变，并设置过滤条件为filter=PASS，并建立tabix索引
log info 'extract mutations from vcf file'
bcftools view -R $database/cds.bed $mutation | bcftools view -f 'PASS,.' | bgzip -c > $out/all.vcf.gz && tabix $out/all.vcf.gz

log info 'get mutated dna sequence'
bcftools consensus -f $database/cds.fa $out/all.vcf.gz -o $out/mutated.fa

log info 'get mutated dna sequence'
seqkit fx2tab $out//mutated.fa > $out/mutated.bed
