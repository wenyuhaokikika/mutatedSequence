source ./log.sh;

release=$1
mkidr source
wget https://ftp.ensembl.org/pub/release-"$release"/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz -P ./source
wget https://ftp.ensembl.org/pub/release-"$release"/fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz -P ./source
wget https://ftp.ensembl.org/pub/release-"$release"/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz -P ./source
wget https://ftp.ensembl.org/pub/release-"$release"/gtf/homo_sapiens/Homo_sapiens.GRCh38."$release".chr_patch_hapl_scaff.gtf.gz -P ./source
wget https://ftp.ensembl.org/pub/release-"$release"/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz -P ./source
