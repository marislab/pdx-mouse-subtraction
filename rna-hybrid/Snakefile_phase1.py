from util.varsub import varsub

configfile: "config.yaml"
varsub(config)

shell.prefix("source ~/.bash_profile")

with open(config['samples']) as f:
	SAMPLES = f.read().splitlines()
	print(SAMPLES)

rule all:
	input:
		config['data']['ref']['genomedir_hybrid'] + "/" + "SAindex",
		expand(config['dirs']['outdirs']['alignstatsdir'] + "{file}" + "/" + "{file}_hybrid_alignstats.txt", file = SAMPLES),
		expand(config['dirs']['outdirs']['bamhybrid'] + "{file}" + "/" + "{file}_star_hg19_final.bam.bai", file = SAMPLES),
		expand(config['dirs']['outdirs']['alignstatsdir'] + "{file}" + "/" + "{file}_hybrid_hg19_alignstats.txt", file = SAMPLES),
		expand(config['dirs']['outdirs']['bamhybrid'] + "{file}" + "/" + "{file}_star_mm10_final.bam.bai", file = SAMPLES),
		expand(config['dirs']['outdirs']['alignstatsdir'] + "{file}" + "/" + "{file}_hybrid_mm10_alignstats.txt", file = SAMPLES),
		expand(config['dirs']['outdirs']['bamhybrid'] + "{file}" + "/" + "{file}_hybrid_hg19_{rep}_sequence.txt.bz2", file = SAMPLES, rep = [1,2])
		
# This pipeline will automatically remove intermediate files that are not required

# STAR genome
rule star_genome:
	input:
		genomefile_hybrid = config['data']['ref']['genomefile_hybrid']
	output:
		genomedir_hybrid = config['data']['ref']['genomedir_hybrid'],
		genomedir_index = config['data']['ref']['genomedir_hybrid'] + "/" + "SAindex"
	log:
		out = config['dirs']['logdir'] + "star_genome.log",
		err = config['dirs']['logdir'] + "star_genome.err"
	params:
		star = config['tools']['star'],
		gtf_hybrid = config['data']['annotation']['gtf_hybrid']
	threads: 4
	shell:
		"""
		mkdir -p {output.genomedir_hybrid}
		
		{params.star} \
		--runMode genomeGenerate \
		--genomeDir {output.genomedir_hybrid} \
		--outFileNamePrefix {output.genomedir_hybrid} \
		--runThreadN 4 \
		--genomeFastaFiles {input.genomefile_hybrid} \
		--sjdbOverhang 100 \
		--sjdbGTFfile {params.gtf_hybrid} 2> {log.err} 1> {log.out}
		"""

# STAR
rule star_align:
	input:
		fq1 = config['dirs']['fastqdir'] + "{file}" + "/" + "{file}_1_sequence.txt.bz2",
		fq2 = config['dirs']['fastqdir'] + "{file}" + "/" + "{file}_2_sequence.txt.bz2",
		genomedir_hybrid = config['data']['ref']['genomedir_hybrid'],
		genomedir_index = config['data']['ref']['genomedir_hybrid'] + "/" + "SAindex"
	output:
		bam = config['dirs']['outdirs']['bamhybrid'] + "{file}" + "/" + "{file}_star_hybrid.bam"
	log:
		out = config['dirs']['logdir'] + "{file}" + "_star_align.log",
		err = config['dirs']['logdir'] + "{file}" + "_star_align.err"
	params:
		bamhybrid = config['dirs']['outdirs']['bamhybrid'] + "{file}" + "/",
		outprefix = config['dirs']['outdirs']['bamhybrid'] + "{file}" + "/" + "{file}_",
		star = config['tools']['star'],
		samtools = config['tools']['samtools']
	threads: 4
	shell:  
		"""
		mkdir -p {params.bamhybrid}
		
		{params.star} \
		--genomeDir {input.genomedir_hybrid} \
		--outFileNamePrefix {params.outprefix} \
		--readFilesIn {input.fq1} {input.fq2} \
		--readFilesCommand "bzcat" \
		--chimSegmentMin 18 \
		--chimScoreMin 12 \
		--runThreadN 8 \
		--outFilterMultimapNmax 20 \
		--outFilterMismatchNoverLmax 0.04 \
		--outFilterIntronMotifs RemoveNoncanonicalUnannotated \
		--alignIntronMax 200000 \
		--outSAMstrandField intronMotif \
		--outStd SAM \
		--outSAMunmapped Within | {params.samtools} view - -b -S -o {output.bam} 2> {log.err} 1> {log.out}
		"""

# alignstats
rule alignstats:
	input:
		bam = config['dirs']['outdirs']['bamhybrid'] + "{file}" + "/" + "{file}_star_hybrid.bam"
	output:
		alignstatsdir = config['dirs']['outdirs']['alignstatsdir'] + "{file}" + "/" + "{file}_hybrid_alignstats.txt"
	log:
		out = config['dirs']['logdir'] + "{file}" + "_alignstats.log",
		err = config['dirs']['logdir'] + "{file}" + "_alignstats.err"
	params:
		alignstatsdir = config['dirs']['outdirs']['alignstatsdir'] + "{file}" + "/",
		alignstats = config['tools']['alignstats'],
		vcrome_bed = config['data']['bed']['vcrome_bed']
	threads: 2
	shell:
		"""
		mkdir -p {params.alignstatsdir}

		{params.alignstats} -v -i {input.bam} -o {output.alignstatsdir} -t {params.vcrome_bed} -C -W 2> {log.err} 1> {log.out}
		"""

# process_bam1
rule process_bam1:
	input:
		bam = config['dirs']['outdirs']['bamhybrid'] + "{file}" + "/" + "{file}_star_hybrid.bam"
	output:
		bam = temp(config['dirs']['outdirs']['bamhybrid'] + "{file}" + "/" + "{file}_star_hg19.bam")
	log:
		# out = config['dirs']['logdir'] + "{file}" + "process_bam1.log",
		err = config['dirs']['logdir'] + "{file}" + "_process_bam1.err"
	params:
		samtools = config['tools']['samtools'],
		hg19_bed = config['data']['bed']['hg19_bed'],
		hg19_fasta = config['data']['ref']['hg19_fasta']
	threads: 2
	shell:
		"""
		{params.samtools} view -L {params.hg19_bed} {input.bam} | grep human | grep -v mouse | sed "s/human//g" | \
		{params.samtools} view -bt {params.hg19_fasta} - > {output.bam} 2> {log.err}
		"""

# process_bam2
rule process_bam2:
	input:
		bam1 = config['dirs']['outdirs']['bamhybrid'] + "{file}" + "/" + "{file}_star_hybrid.bam",
		bam2 = config['dirs']['outdirs']['bamhybrid'] + "{file}" + "/" + "{file}_star_hg19.bam"
	output:
		bam = temp(config['dirs']['outdirs']['bamhybrid'] + "{file}" + "/" + "{file}_star_hg19_header_correct.bam")
	log:
		# out = config['dirs']['logdir'] + "{file}" + "process_bam2.log",
		err = config['dirs']['logdir'] + "{file}" + "_process_bam2.err"
	params:
		samtools = config['tools']['samtools'],
	threads: 2
	shell:
		"""
		{params.samtools} view -H {input.bam1} | sed "/@SQ\tSN:mouse*/d" | sed "s/human//g" | \
		{params.samtools} reheader - {input.bam2} > {output.bam} 2> {log.err}
		"""

# process_bam3
rule process_bam3:
	input:
		bam = config['dirs']['outdirs']['bamhybrid'] + "{file}" + "/" + "{file}_star_hg19_header_correct.bam"
	output:
		bam = temp(config['dirs']['outdirs']['bamhybrid'] + "{file}" + "/" + "{file}_star_hg19_sorted_by_name.bam")
	log:
		out = config['dirs']['logdir'] + "{file}" + "_process_bam3.log",
		err = config['dirs']['logdir'] + "{file}" + "_process_bam3.err"
	params:
		sambamba = config['tools']['sambamba'],
		tmpdir = config['dirs']['tmpdir']
	threads: 4
	shell:
		"""
		{params.sambamba} sort -l 0 -t 8 -n --tmpdir {params.tmpdir} -o {output.bam} {input.bam} 2> {log.err} 1> {log.out}
		"""

rule process_bam4:
	input:
		bam = config['dirs']['outdirs']['bamhybrid'] + "{file}" + "/" + "{file}_star_hg19_sorted_by_name.bam"
	output:
		bam = temp(config['dirs']['outdirs']['bamhybrid'] + "{file}" + "/" + "{file}_star_hg19_fixmate.bam")
	log:
		out = config['dirs']['logdir'] + "{file}" + "_process_bam4.log",
		err = config['dirs']['logdir'] + "{file}" + "_process_bam4.err"
	params:
		samtools = config['tools']['samtools']
	threads: 2
	shell:
		"""
		{params.samtools} view -F 256 -bh {input.bam} | {params.samtools} fixmate - {output.bam} 2> {log.err} 1> {log.out}
		"""

rule process_bam5:
	input:
		bam = config['dirs']['outdirs']['bamhybrid'] + "{file}" + "/" + "{file}_star_hg19_fixmate.bam"
	output:
		bam = temp(config['dirs']['outdirs']['bamhybrid'] + "{file}" + "/" + "{file}_star_hg19_sorted.bam")
	log:
		out = config['dirs']['logdir'] + "{file}" + "_process_bam5.log",
		err = config['dirs']['logdir'] + "{file}" + "_process_bam5.err"
	params:
		sambamba = config['tools']['sambamba'],
		tmpdir = config['dirs']['tmpdir']
	threads: 4
	shell:
		"""
		{params.sambamba} sort -l 0 -t 8 --tmpdir {params.tmpdir} -o {output.bam} {input.bam} 2> {log.err} 1> {log.out}
		"""

rule process_bam6:
	input:
		bam = config['dirs']['outdirs']['bamhybrid'] + "{file}" + "/" + "{file}_star_hg19_sorted.bam"
	output:
		bam1 = temp(config['dirs']['outdirs']['bamhybrid'] + "{file}" + "/" + "{file}_star_hg19_unpaired.bam"),
		bam2 = temp(config['dirs']['outdirs']['bamhybrid'] + "{file}" + "/" + "{file}_star_hg19_sorted_paired.bam")
	log:
		# out = config['dirs']['logdir'] + "{file}" + "process_bam6.log",
		err = config['dirs']['logdir'] + "{file}" + "_process_bam6.err"
	params:
		samtools = config['tools']['samtools']
	threads: 2
	shell:
		"""
		{params.samtools} view -F 1 -bh {input.bam} > {output.bam1} 2> {log.err}
		{params.samtools} view -f 1 -bh {input.bam} > {output.bam2} 2> {log.err}
		"""

rule process_bam7:
	input:
		bam = config['dirs']['outdirs']['bamhybrid'] + "{file}" + "/" + "{file}_star_hg19_sorted_paired.bam"
	output:
		bam = config['dirs']['outdirs']['bamhybrid'] + "{file}" + "/" + "{file}_star_hg19_final.bam"
	log:
		out = config['dirs']['logdir'] + "{file}" + "_process_bam7.log",
		err = config['dirs']['logdir'] + "{file}" + "_process_bam7.err"
	params:
		sambamba = config['tools']['sambamba'],
		tmpdir = config['dirs']['tmpdir']
	threads: 4
	shell:
		"""
		{params.sambamba} markdup -t 8 -p --tmpdir={params.tmpdir} {input.bam} {output.bam} 2> {log.err} 1> {log.out}
		"""

rule process_bam8:
	input:
		bam = config['dirs']['outdirs']['bamhybrid'] + "{file}" + "/" + "{file}_star_hg19_final.bam"
	output:
		bai = config['dirs']['outdirs']['bamhybrid'] + "{file}" + "/" + "{file}_star_hg19_final.bam.bai"
	log:
		out = config['dirs']['logdir'] + "{file}" + "_process_bam8.log",
		err = config['dirs']['logdir'] + "{file}" + "_process_bam8.err"
	params:
		sambamba = config['tools']['sambamba']
	threads: 2
	shell:
		"""
		{params.sambamba} index -t 8 {input.bam} 2> {log.err} 1> {log.out}
		"""

rule process_bam9:
	input:
		bam = config['dirs']['outdirs']['bamhybrid'] + "{file}" + "/" + "{file}_star_hg19_final.bam"
	output:
		alignstatsdir = config['dirs']['outdirs']['alignstatsdir'] + "{file}" + "/" + "{file}_hybrid_hg19_alignstats.txt"
	log:
		out = config['dirs']['logdir'] + "{file}" + "_process_bam9.log",
		err = config['dirs']['logdir'] + "{file}" + "_process_bam9.err"
	params:
		alignstats = config['tools']['alignstats'],
		vcrome_bed = config['data']['bed']['vcrome_bed']
	threads: 2
	shell:
		"""
		{params.alignstats} -v -i {input.bam} -o {output.alignstatsdir} -t {params.vcrome_bed} -C -W 2> {log.err} 1> {log.out}
		"""

rule process_bam10:
	input:
		bam = config['dirs']['outdirs']['bamhybrid'] + "{file}" + "/" + "{file}_star_hybrid.bam"
	output:
		bam = temp(config['dirs']['outdirs']['bamhybrid'] + "{file}" + "/" + "{file}_star_mm10.bam")
	log:
		# out = config['dirs']['logdir'] + "{file}" + "process_bam10.log",
		err = config['dirs']['logdir'] + "{file}" + "_process_bam10.err"
	params:
		samtools = config['tools']['samtools'],
		mm10_bed = config['data']['bed']['mm10_bed'],
		mm10_fasta = config['data']['ref']['mm10_fasta']
	threads: 2
	shell:
		"""
		{params.samtools} view -L {params.mm10_bed} {input.bam} | grep mouse | grep -v human | sed "s/mouse/chr/g" | \
		{params.samtools} view -bt {params.mm10_fasta} - > {output.bam} 2> {log.err}
		"""

rule process_bam11:
	input:
		bam1 = config['dirs']['outdirs']['bamhybrid'] + "{file}" + "/" + "{file}_star_hybrid.bam",
		bam2 = config['dirs']['outdirs']['bamhybrid'] + "{file}" + "/" + "{file}_star_mm10.bam"
	output:
		bam = temp(config['dirs']['outdirs']['bamhybrid'] + "{file}" + "/" + "{file}_star_mm10_header_correct.bam")
	log:
		# out = config['dirs']['logdir'] + "{file}" + "process_bam11.log",
		err = config['dirs']['logdir'] + "{file}" + "_process_bam11.err"
	params:
		samtools = config['tools']['samtools']
	threads: 2
	shell:
		"""
		{params.samtools} view -H {input.bam1} | sed "/@SQ\tSN:human*/d" | sed "s/mouse/chr/g"  | \
		{params.samtools} reheader - {input.bam2} > {output.bam} 2> {log.err}
		"""

rule process_bam12:
	input:
		bam = config['dirs']['outdirs']['bamhybrid'] + "{file}" + "/" + "{file}_star_mm10_header_correct.bam"
	output:
		bam = temp(config['dirs']['outdirs']['bamhybrid'] + "{file}" + "/" + "{file}_star_mm10_sorted_by_name.bam")
	log:
		out = config['dirs']['logdir'] + "{file}" + "_process_bam12.log",
		err = config['dirs']['logdir'] + "{file}" + "_process_bam12.err"
	params:
		sambamba = config['tools']['sambamba'],
		tmpdir = config['dirs']['tmpdir']
	threads: 4
	shell:
		"""
		{params.sambamba} sort -l 0 -t 8 -n --tmpdir {params.tmpdir} -o {output.bam} {input.bam} 2> {log.err} 1> {log.out}
		"""

rule process_bam13:
	input:
		bam = config['dirs']['outdirs']['bamhybrid'] + "{file}" + "/" + "{file}_star_mm10_sorted_by_name.bam"
	output:
		bam = temp(config['dirs']['outdirs']['bamhybrid'] + "{file}" + "/" + "{file}_star_mm10_fixmate.bam")
	log:
		out = config['dirs']['logdir'] + "{file}" + "_process_bam13.log",
		err = config['dirs']['logdir'] + "{file}" + "_process_bam13.err"
	params:
		samtools = config['tools']['samtools']
	threads: 2
	shell:
		"""
		{params.samtools} view -F 256 -bh {input.bam} | {params.samtools} fixmate - {output.bam} 2> {log.err} 1> {log.out}
		"""	

rule process_bam14:
	input:
		bam = config['dirs']['outdirs']['bamhybrid'] + "{file}" + "/" + "{file}_star_mm10_fixmate.bam"
	output:
		bam = temp(config['dirs']['outdirs']['bamhybrid'] + "{file}" + "/" + "{file}_star_mm10_sorted.bam")
	log:
		out = config['dirs']['logdir'] + "{file}" + "_process_bam14.log",
		err = config['dirs']['logdir'] + "{file}" + "_process_bam14.err"
	params:
		sambamba = config['tools']['sambamba'],
		tmpdir = config['dirs']['tmpdir']
	threads: 4
	shell:
		"""
		{params.sambamba} sort -l 0 -t 8 --tmpdir {params.tmpdir} -o {output.bam} {input.bam} 2> {log.err} 1> {log.out}
		"""

rule process_bam15:
	input:
		bam = config['dirs']['outdirs']['bamhybrid'] + "{file}" + "/" + "{file}_star_mm10_sorted.bam"
	output:
		bam = config['dirs']['outdirs']['bamhybrid'] + "{file}" + "/" + "{file}_star_mm10_final.bam"
	log:
		out = config['dirs']['logdir'] + "{file}" + "_process_bam15.log",
		err = config['dirs']['logdir'] + "{file}" + "_process_bam15.err"
	params:
		sambamba = config['tools']['sambamba'],
		tmpdir = config['dirs']['tmpdir']
	threads: 2
	shell:
		"""
		{params.sambamba} markdup -t 8 -p --tmpdir={params.tmpdir} {input.bam} {output.bam} 2> {log.err} 1> {log.out}
		"""

rule process_bam16:
	input:
		bam = config['dirs']['outdirs']['bamhybrid'] + "{file}" + "/" + "{file}_star_mm10_final.bam"
	output:
		bai = config['dirs']['outdirs']['bamhybrid'] + "{file}" + "/" + "{file}_star_mm10_final.bam.bai"
	log:
		out = config['dirs']['logdir'] + "{file}" + "_process_bam16.log",
		err = config['dirs']['logdir'] + "{file}" + "_process_bam16.err"
	params:
		sambamba = config['tools']['sambamba']
	threads: 2
	shell:
		"""
		{params.sambamba} index -t 8 {input.bam} 2> {log.err} 1> {log.out}
		"""

rule process_bam17:
	input:
		bam = config['dirs']['outdirs']['bamhybrid'] + "{file}" + "/" + "{file}_star_mm10_final.bam"
	output:
		alignstatsdir = config['dirs']['outdirs']['alignstatsdir'] + "{file}" + "/" + "{file}_hybrid_mm10_alignstats.txt"
	log:
		out = config['dirs']['logdir'] + "{file}" + "_process_bam17.log",
		err = config['dirs']['logdir'] + "{file}" + "_process_bam17.err"
	params:
		alignstats = config['tools']['alignstats'],
		vcrome_bed = config['data']['bed']['vcrome_bed']
	threads: 2
	shell:
		"""
		{params.alignstats} -v -i {input.bam} -o {output.alignstatsdir} -t {params.vcrome_bed} -C -W 2> {log.err} 1> {log.out}
		"""

# convert hybrid bam to fastq and compress
rule bam2fastq:
	input:
		bam = config['dirs']['outdirs']['bamhybrid'] + "{file}" + "/" + "{file}_star_hg19_final.bam"
	output:
		fq1 = config['dirs']['outdirs']['bamhybrid'] + "{file}" + "/" + "{file}_hybrid_hg19_1_sequence.txt",
		fq2 = config['dirs']['outdirs']['bamhybrid'] + "{file}" + "/" + "{file}_hybrid_hg19_2_sequence.txt"
	log:
		out = config['dirs']['logdir'] + "{file}" + "_bam2fastq.log",
		err = config['dirs']['logdir'] + "{file}" + "_bam2fastq.err"
	params:
		java = config['binaries']['java'],
		picard = config['tools']['picard']
	threads: 4
	shell:
		"""
		{params.java} -Xmx4g -jar {params.picard}/picard.jar SamToFastq INPUT={input.bam} FASTQ={output.fq1} SECOND_END_FASTQ={output.fq2} 2> {log.err} 1> {log.out}
		"""

# compress fastqs
rule compressfastq:
	input:
		fq = config['dirs']['outdirs']['bamhybrid'] + "{file}" + "/" + "{file}_hybrid_hg19_{rep}_sequence.txt"
	output:
		fq = config['dirs']['outdirs']['bamhybrid'] + "{file}" + "/" + "{file}_hybrid_hg19_{rep}_sequence.txt.bz2",
	params:
		bzip2 = config['binaries']['bzip2']
	threads: 2
	shell:
		"""
		{params.bzip2} -c {input.fq} > {output.fq}
		"""



