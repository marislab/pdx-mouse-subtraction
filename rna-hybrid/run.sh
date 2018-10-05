export PROJDIR="/mnt/isilon/maris_lab/target_nbl_ngs/PPTC-PDX-genomics/mouse_subtraction_pipeline/rna-hybrid/results"
export SCRIPTSDIR="/mnt/isilon/maris_lab/target_nbl_ngs/PPTC-PDX-genomics/mouse_subtraction_pipeline/scripts"
export FQDIR="/mnt/isilon/maris_lab/target_nbl_ngs/PPTC-PDX-genomics/mouse_subtraction_pipeline/data/fastq"
export sample="PPTC-AF02-XTP1-B-1-0-R"
export FQ1="PPTC-AF02-XTP1-B-1-0-R_1_sequence.txt.bz2"
export FQ2="PPTC-AF02-XTP1-B-1-0-R_2_sequence.txt.bz2"
export VALENCE="/mnt/isilon/maris_lab/target_nbl_ngs/PPTC-PDX-genomics/mouse_subtraction_pipeline/rna-hybrid/valence.txt"

sh $SCRIPTSDIR/GenerateRNASeq.py \
dir $FQDIR \
-sample $sample \
-library $LIBNAME \
-fastq $FQ1 \
-fastqpair $FQ2 \
-masterxml $SCRIPTSDIR/RNASeq_STAR.xml  \
-mastercsv $SCRIPTSDIR/RNASeq_STAR.csv \
-outdir $PROJDIR/raw -valence $VALENCE


