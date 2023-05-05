###
 # @Description: download bio sequence from ensembl
 # @version: 
 # @Author: wenyuhao
 # @Date: 2023-05-05 11:15:40
 # @LastEditors: wenyuhao
 # @LastEditTime: 2023-05-05 11:15:40
### 

source ./log.sh; #https://github.com/Zordrak/bashlog
check_status () {
    if [ $? -ne 0 ]; then
        log error "Error: download ERROR."
        exit 1
    else
        log info "download successfully."
    fi
}
release=$1
mkidr resource
log info 'download cdna ...'
wget https://ftp.ensembl.org/pub/release-"$release"/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz -P ./resource
check_status
log info 'download cds ...'
wget https://ftp.ensembl.org/pub/release-"$release"/fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz -P ./resource
check_status 
log info 'download proteins ...'
wget https://ftp.ensembl.org/pub/release-"$release"/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz -P ./resource
check_status
log info 'download gtf ...'
wget https://ftp.ensembl.org/pub/release-"$release"/gtf/homo_sapiens/Homo_sapiens.GRCh38."$release".chr_patch_hapl_scaff.gtf.gz -P ./resource
check_status 
log info 'download dna ...'
wget https://ftp.ensembl.org/pub/release-"$release"/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz -P ./resource
check_status

