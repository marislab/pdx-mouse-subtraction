# This script creates data directories and downloads reference data
# set your rootdirectory and reference directory by refering to the yaml files
rootdir="${ROOT_DIR}" 
refdir="${REF_DIR}"

wget --output-document='$refdir/star_fusion_db/hg19/GRCh37_gencode_v19_CTAT_lib_July192017/ctat_genome_lib_build_dir' \
https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/GRCh37_gencode_v19_CTAT_lib_July192017.source_data.tar.gz

wget --output-document='$rootdir/data/ref/mm10_hg19/hg19.bed' \
https://ndownloader.figshare.com/files/14515394

wget --output-document='$rootdir/data/ref/mm10_hg19/mm10.bed' \
https://ndownloader.figshare.com/files/14515889

wget --output-document='$rootdir/data/capture_designs/HG19_vcrome2.1.bed' \
https://ndownloader.figshare.com/files/14515391

wget --output-document='$rootdir/data/ref/hg19_genome/Homo_sapiens.GRCh37.71.hap.ERCC.sm.fa' \
https://ndownloader.figshare.com/files/14515886

wget --output-document='$rootdir/data/ref/mm10_hg19/hg19_mm10.fa' \
https://ndownloader.figshare.com/files/14515310

wget --output-document='$rootdir/data/gtf/hg19_mm10.gtf' \
https://ndownloader.figshare.com/files/14515388

wget --output-document='$rootdir/data/gtf/Homo_sapiens.GRCh37.71.hap.ERCC.gtf' \
https://ndownloader.figshare.com/files/14570963

wget --output-document='$rootdir/data/gtf/Homo_sapiens.GRCh37.71.hap.gene_name.gtf' \
https://ndownloader.figshare.com/files/14570981

wget --output-document='$rootdir/data/ERCC92/ERCC_conc.tsv' \
https://ndownloader.figshare.com/files/14570921

wget --output-document='$rootdir/data/ref/transcripts/protein_coding_canonical.T_chr.fa' \
https://ndownloader.figshare.com/files/14570957
