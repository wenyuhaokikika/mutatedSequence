#!/bin/bash
###
 # @Description: 
 # @version: 
 # @Author: wenyuhao
 # @Date: 2023-05-07 03:05:42
 # @LastEditors: wenyuhao
 # @LastEditTime: 2023-05-16 09:44:58
### 
#从bed文件获取fasta文件
filepath=$(cd "$(dirname "$0")"; pwd) # 获取当前文件夹的绝对路径
#echo 'mutatedSequence' | figlet -t -c > /data/wenyuhao/55/ensembl/ensembl109/db/mutatedSequence/log.txt
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
# 如果cds输入未bed文件将他转化为fa文件，如果是fa文件则不改变。
#if [[ $cdsfa =~ 'bed' ]]

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

# 只选择在bed范围内的突变，并设置过滤条件为filter=PASS，并建立tabix索引
log info 'extract mutations from vcf file'
if [ ! -f $mutation.tbi ];then tabix $mutation;fi
#tabix -f $mutation
bcftools view -R $database/cds.bed $mutation | bcftools view -f 'PASS,.' | bgzip -c > $out/all.vcf.gz && tabix $out/all.vcf.gz

log info 'get mutated dna sequence'
bcftools consensus -f $database/cds.fa $out/all.vcf.gz -o $out/mutated.fa > $out/consensusLog.txt 2>&1
check_status

log info 'get mutated dna sequence'
seqkit fx2tab $out/mutated.fa > $out/mutated.bed
