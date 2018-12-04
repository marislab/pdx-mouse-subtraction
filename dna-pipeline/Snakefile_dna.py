from util.varsub import varsub

configfile: "config_dna.yaml"
varsub(config)

shell.prefix("source ~/.bash_profile")

with open(config['samples']) as f:
	SAMPLES = f.read().splitlines()
	print(SAMPLES)

rule all:
	input:
		config['data']['ref']['genomefile_hybrid'] + ".sa",
		expand(config['dirs']['outdirs']['realignstatsdir'] + "{file}" + "/" + "{file}_hybrid_alignstats.txt", file = SAMPLES),
		expand(config['dirs']['outdirs']['realign'] + "{file}.hybrid_hg19.bam.bai", file = SAMPLES),
		expand(config['dirs']['outdirs']['realignstatsdir'] + "{file}" + "/" + "{file}_hybrid_hg19_alignstats.txt", file = SAMPLES),
		expand(config['dirs']['outdirs']['realign'] + "{file}.hybrid_mm10.bam.bai", file = SAMPLES),
		expand(config['dirs']['outdirs']['realignstatsdir'] + "{file}" + "/" + "{file}_hybrid_mm10_alignstats.txt", file = SAMPLES),
		expand(config['dirs']['outdirs']['realign'] + "{file}.hybrid_hg19_mm10.bam.bai", file = SAMPLES),
		expand(config['dirs']['outdirs']['realignstatsdir'] + "{file}" + "/" + "{file}_hybrid_hg19_mm10_alignstats.txt", file = SAMPLES)

# convert to fastq
rule bam2fastq:
	input:
		bam = config['dirs']['outdirs']['dna_bamdir'] + "{file}" + "/" + "{file}-PDX.bam"
	output:
		fq1 = temp(config['dirs']['outdirs']['dna_fastqdir'] + "{file}" + "/" + "{file}_1.fastq"),
		fq2 = temp(config['dirs']['outdirs']['dna_fastqdir'] + "{file}" + "/" + "{file}_2.fastq")
	log:
		out = config['dirs']['logdir'] + "{file}" + "_bam2fastq.log",
		err = config['dirs']['logdir'] + "{file}" + "_bam2fastq.err"
	params:
		java = config['binaries']['java'],
		picard = config['tools']['picard']
	threads: 4
	shell:
		"""
		{params.java} -Xmx4g -jar {params.picard}/picard.jar SamToFastq \
		VALIDATION_STRINGENCY=LENIENT \
		INPUT={input.bam} FASTQ={output.fq1} SECOND_END_FASTQ={output.fq2} 2> {log.err} 1> {log.out}
		"""

rule bwa_index:
	input:
		genomefile_hybrid = config['data']['ref']['genomefile_hybrid']
	output:
		genomefile_hybrid = config['data']['ref']['genomefile_hybrid'] + ".sa"
	log:
		out = config['dirs']['logdir'] + "_bwa_index.log",
		err = config['dirs']['logdir'] + "_bwa_index.err"
	params:
		bwa = config['tools']['bwa']
	threads: 2
	shell:
		"""
		{params.bwa} index {input.genomefile_hybrid} 2> {log.err} 1> {log.out}
		"""

# align each fastq separately
rule bwa_aln:
	input:
		fq = config['dirs']['outdirs']['dna_fastqdir'] + "{file}" + "/" + "{file}_{rep}.fastq",
	output:
		sai = temp(config['dirs']['outdirs']['realign'] + "{file}" + "/" + "{file}_{rep}.fastq.sai"),
	log:
		err = config['dirs']['logdir'] + "{file}" + "_bwa_aln.err"
	params:
		bwa = config['tools']['bwa'],
		genomefile_hybrid = config['data']['ref']['genomefile_hybrid']
	threads: 4
	shell:
		"""
		{params.bwa} aln -t 8 {params.genomefile_hybrid} {input.fq} > {output.sai} 2> {log.err}
		"""

# bwa sampe
rule bwa_sampe:
	input:
		sai1 = config['dirs']['outdirs']['realign'] + "{file}" + "/" + "{file}_1.fastq.sai",
		sai2 = config['dirs']['outdirs']['realign'] + "{file}" + "/" + "{file}_2.fastq.sai",
		fq1 = config['dirs']['outdirs']['dna_fastqdir'] + "{file}" + "/" + "{file}_1.fastq",
		fq2 = config['dirs']['outdirs']['dna_fastqdir'] + "{file}" + "/" + "{file}_2.fastq",
	output:
		bam = temp(config['dirs']['outdirs']['realign'] + "{file}.tmp.bam")
	log:
		out = config['dirs']['logdir'] + "{file}" + "_bwa_sampe.log",
		err = config['dirs']['logdir'] + "{file}" + "_bwa_sampe.err"
	params:
		sample = '{file}',
		bwa = config['tools']['bwa'],
		samtools = config['tools']['samtools'],
		sambamba = config['tools']['sambamba'],
		genomefile_hybrid = config['data']['ref']['genomefile_hybrid'],
		tmpdir = config['dirs']['tmpdir']
	threads: 2
	shell:
		"""
		{params.bwa} sampe -P -r '@RG\tID:0\tSM:{params.sample}\tPL:Illumina' \
		{params.genomefile_hybrid} {input.sai1} {input.sai2} {input.fq1} {input.fq2} | \
		{params.samtools} view -F 256 -bu - | {params.samtools} fixmate - - | \
		{params.sambamba} sort -l 1 -t 8 --tmpdir {params.tmpdir} -o {output.bam} 2> {log.err} 1> {log.out}
		"""

# sambamba
rule run_sambamba:
	input:
		bam = config['dirs']['outdirs']['realign'] + "{file}.tmp.bam"
	output:
		bam = temp(config['dirs']['outdirs']['realign'] + "{file}.tmp.markdups.bam")
	log:
		out = config['dirs']['logdir'] + "{file}" + "_run_sambamba.log",
		err = config['dirs']['logdir'] + "{file}" + "_run_sambamba.err"
	params:
		tmpdir = config['dirs']['tmpdir'],
		sambamba = config['tools']['sambamba']
	threads: 2
	shell:
		"""
		{params.sambamba} markdup -t 8 -p --tmpdir={params.tmpdir} {input.bam} {output.bam} 2> {log.err} 1> {log.out}
		"""

# sambamba_index
rule sambamba_index:
	input:
		bam = config['dirs']['outdirs']['realign'] + "{file}.tmp.markdups.bam"
	output:
		bai = config['dirs']['outdirs']['realign'] + "{file}.tmp.markdups.bam.bai"
	log:
		out = config['dirs']['logdir'] + "{file}" + "_sambamba_index.log",
		err = config['dirs']['logdir'] + "{file}" + "_sambamba_index.err"
	params:
		sambamba = config['tools']['sambamba']
	threads: 2
	shell:
		"""
		{params.sambamba} index -t 8 {input.bam} 2> {log.err} 1> {log.out}
		"""

# validate bam
rule bam_validate:
	input:
		bam = config['dirs']['outdirs']['realign'] + "{file}.tmp.markdups.bam"
	output:
		txt = config['dirs']['outdirs']['realign'] + "{file}_validate.txt"
	log:
		err = config['dirs']['logdir'] + "{file}" + "_bam_validate.err"
	params:
		bamutil = config['tools']['bamutil']
	threads: 2
	shell:
		"""
		{params.bamutil} validate --in {input.bam} --verbose > {output.txt} 2> {log.err}
		"""

# align stats
rule realign_stats:
	input:
		bam = config['dirs']['outdirs']['realign'] + "{file}.tmp.markdups.bam"
	output:
		realignstats = config['dirs']['outdirs']['realignstatsdir'] + "{file}" + "/" + "{file}_hybrid_alignstats.txt"
	log:
		out = config['dirs']['logdir'] + "{file}" + "_realign_stats.log",
		err = config['dirs']['logdir'] + "{file}" + "_realign_stats.err"
	params:
		alignstats = config['tools']['alignstats'],
		realignstats = config['dirs']['outdirs']['realignstatsdir'] + "{file}" + "/",
		vcrome_bed = config['data']['bed']['vcrome_bed']
	threads: 2
	shell:
		"""
		mkdir -p {params.realignstats}

		{params.alignstats} -v -i {input.bam} -o {output.realignstats} -t {params.vcrome_bed} -C -W 2> {log.err} 1> {log.out}
		"""

rule process_bam1:
	input:
		bam = config['dirs']['outdirs']['realign'] + "{file}.tmp.markdups.bam"
	output:
		bam = temp(config['dirs']['outdirs']['realign'] + "{file}.hg19_tmp.bam")
	log:
		err = config['dirs']['logdir'] + "{file}" + "_process_bam1.err"
	params:
		samtools = config['tools']['samtools'],
		hg19_bed = config['data']['bed']['hg19_bed'],
		hg19_fasta = config['data']['ref']['hg19_fasta']
	threads: 4
	shell:
		"""
		{params.samtools} view -L {params.hg19_bed} {input.bam} | grep human | grep -v mouse | sed "s/human//g" | \
		{params.samtools} view -bt {params.hg19_fasta} - > {output.bam} 2> {log.err}
		"""

rule process_bam2:
	input:
		bam1 = config['dirs']['outdirs']['realign'] + "{file}.tmp.markdups.bam",
		bam2 = config['dirs']['outdirs']['realign'] + "{file}.hg19_tmp.bam"
	output:
		bam = temp(config['dirs']['outdirs']['realign'] + "{file}_header_correct_hg19.bam")
	log:
		err = config['dirs']['logdir'] + "{file}" + "_process_bam2.err"
	params:
		samtools = config['tools']['samtools']
	threads: 2
	shell:
		"""
		{params.samtools} view -H {input.bam1} | sed "/@SQ\tSN:mouse*/d" | sed "s/human//g" | \
		{params.samtools} reheader - {input.bam2} > {output.bam} 2> {log.err}
		"""

rule process_bam3:
	input:
		bam = config['dirs']['outdirs']['realign'] + "{file}_header_correct_hg19.bam"
	output:
		bam = temp(config['dirs']['outdirs']['realign'] + "{file}_hg19_sorted_by_name.bam")
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
		bam = config['dirs']['outdirs']['realign'] + "{file}_hg19_sorted_by_name.bam"
	output:
		bam = temp(config['dirs']['outdirs']['realign'] + "{file}_hg19_sorted.bam")
	log:
		out = config['dirs']['logdir'] + "{file}" + "_process_bam4.log",
		err = config['dirs']['logdir'] + "{file}" + "_process_bam4.err"
	params:
		samtools = config['tools']['samtools'],
		sambamba = config['tools']['sambamba'],
		tmpdir = config['dirs']['tmpdir']
	threads: 2
	shell:
		"""
		{params.samtools} view -F 256 -Sbh {input.bam} | \
		{params.samtools} fixmate - - | \
		{params.sambamba} sort -l 0 -t 8 --tmpdir {params.tmpdir} -o {output.bam} 2> {log.err} 1> {log.out}
		"""

rule process_bam5:
	input:
		bam = config['dirs']['outdirs']['realign'] + "{file}_hg19_sorted.bam"
	output:
		bam = temp(config['dirs']['outdirs']['realign'] + "{file}.hybrid_unpaired_hg19.bam")
	log:
		err = config['dirs']['logdir'] + "{file}" + "_process_bam5.err"
	params:
		samtools = config['tools']['samtools']
	threads: 2
	shell:
		"""
		{params.samtools} view -F 1 -Sbh {input.bam} > {output.bam} 2> {log.err}
		"""

rule process_bam6:
	input:
		bam = config['dirs']['outdirs']['realign'] + "{file}_hg19_sorted.bam"
	output:
		bam = temp(config['dirs']['outdirs']['realign'] + "{file}.hybrid_paired_hg19.bam")
	log:
		err = config['dirs']['logdir'] + "{file}" + "_process_bam6.err"
	params:
		samtools = config['tools']['samtools']
	threads: 2
	shell:
		"""
		{params.samtools} view -f 1 -Sbh {input.bam} > {output.bam} 2> {log.err}
		"""

rule process_bam7:
	input:
		bam = config['dirs']['outdirs']['realign'] + "{file}.hybrid_paired_hg19.bam"
	output:
		bam = temp(config['dirs']['outdirs']['realign'] + "{file}.hybrid_hg19.bam")
	log:
		out = config['dirs']['logdir'] + "{file}" + "_process_bam7.log",
		err = config['dirs']['logdir'] + "{file}" + "_process_bam7.err"
	params:
		sambamba = config['tools']['sambamba'],
		tmpdir = config['dirs']['tmpdir']
	threads: 2
	shell:
		"""
		{params.sambamba} markdup -t 8 -p --tmpdir={params.tmpdir} {input.bam} {output.bam} 2> {log.err} 1> {log.out} 
		"""

rule process_bam8:
	input:
		bam = config['dirs']['outdirs']['realign'] + "{file}.hybrid_hg19.bam"
	output:
		bam = config['dirs']['outdirs']['realign'] + "{file}.hybrid_hg19.bam.bai"
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
		bam = config['dirs']['outdirs']['realign'] + "{file}.hybrid_hg19.bam"
	output:
		realignstats = config['dirs']['outdirs']['realignstatsdir'] + "{file}" + "/" + "{file}_hybrid_hg19_alignstats.txt"
	log:
		out = config['dirs']['logdir'] + "{file}" + "_process_bam9.log",
		err = config['dirs']['logdir'] + "{file}" + "_process_bam9.err"
	params:
		alignstats = config['tools']['alignstats'],
		vcrome_bed = config['data']['bed']['vcrome_bed']
	threads: 2
	shell:
		"""
		{params.alignstats} -v -i {input.bam} -o {output.realignstats} -t {params.vcrome_bed} -W 2> {log.err} 1> {log.out}
		"""

rule process_bam10:
	input:
		bam = config['dirs']['outdirs']['realign'] + "{file}.tmp.markdups.bam"
	output:
		bam = temp(config['dirs']['outdirs']['realign'] + "{file}.mm10_tmp.bam")
	log:
		err = config['dirs']['logdir'] + "{file}" + "_process_bam10.err"
	params:
		mm10_bed = config['data']['bed']['mm10_bed'],
		mm10_fasta = config['data']['ref']['mm10_fasta'],
		samtools = config['tools']['samtools'],
	threads: 2
	shell:
		"""
		{params.samtools} view -L {params.mm10_bed} {input.bam} | grep mouse | grep -v human | sed "s/mouse/chr/g" | \
		{params.samtools} view -bt {params.mm10_fasta} - > {output.bam} 2> {log.err}
		"""

rule process_bam11:
	input:
		bam1 = config['dirs']['outdirs']['realign'] + "{file}.tmp.markdups.bam",
		bam2 = config['dirs']['outdirs']['realign'] + "{file}.mm10_tmp.bam"
	output:
		bam = temp(config['dirs']['outdirs']['realign'] + "{file}_header_correct_mm10.bam")
	log:
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
		bam = config['dirs']['outdirs']['realign'] + "{file}_header_correct_mm10.bam"
	output:
		bam = temp(config['dirs']['outdirs']['realign'] + "{file}_mm10_sorted_by_name.bam")
	log:
		out = config['dirs']['logdir'] + "{file}" + "_process_bam12.log",
		err = config['dirs']['logdir'] + "{file}" + "_process_bam12.err"
	params:
		sambamba = config['tools']['sambamba'],
		tmpdir = config['dirs']['tmpdir']
	threads: 2
	shell:
		"""
		{params.sambamba} sort -l 0 -t 8 -n --tmpdir {params.tmpdir} -o {output.bam} {input.bam} 2> {log.err} 1> {log.out}
		"""

rule process_bam13:
	input:
		bam = config['dirs']['outdirs']['realign'] + "{file}_mm10_sorted_by_name.bam"
	output:
		bam = temp(config['dirs']['outdirs']['realign'] + "{file}_mm10_final.bam")
	log:
		out = config['dirs']['logdir'] + "{file}" + "_process_bam13.log",
		err = config['dirs']['logdir'] + "{file}" + "_process_bam13.err"
	params:
		samtools = config['tools']['samtools'],
		sambamba = config['tools']['sambamba'],
		tmpdir = config['dirs']['tmpdir']
	threads: 2
	shell:
		"""
		{params.samtools} view -F 256 -Sbh {input.bam} | \
		{params.samtools} fixmate - - | {params.sambamba} sort -l 0 -t 8 --tmpdir {params.tmpdir} -o {output.bam} 2> {log.err} 1> {log.out} 
		"""

rule process_bam14:
	input:
		bam = config['dirs']['outdirs']['realign'] + "{file}_mm10_final.bam"
	output:
		bam = temp(config['dirs']['outdirs']['realign'] + "{file}.hybrid_mm10.bam")
	log:
		out = config['dirs']['logdir'] + "{file}" + "_process_bam14.log",
		err = config['dirs']['logdir'] + "{file}" + "_process_bam14.err"
	params:
		sambamba = config['tools']['sambamba'],
		tmpdir = config['dirs']['tmpdir']
	threads: 2
	shell:
		"""
		{params.sambamba} markdup -t 8 -p --tmpdir={params.tmpdir} {input.bam} {output.bam} 2> {log.err} 1> {log.out} 
		"""

rule process_bam15:
	input:
		bam = config['dirs']['outdirs']['realign'] + "{file}.hybrid_mm10.bam"
	output:
		bam = config['dirs']['outdirs']['realign'] + "{file}.hybrid_mm10.bam.bai"
	log:
		out = config['dirs']['logdir'] + "{file}" + "_process_bam15.log",
		err = config['dirs']['logdir'] + "{file}" + "_process_bam15.err"
	params:
		sambamba = config['tools']['sambamba']
	threads: 2
	shell:
		"""
		{params.sambamba} index -t 8 {input.bam} 2> {log.err} 1> {log.out} 
		"""

rule process_bam16:
	input:
		bam = config['dirs']['outdirs']['realign'] + "{file}.hybrid_mm10.bam"
	output:
		realignstats = config['dirs']['outdirs']['realignstatsdir'] + "{file}" + "/" + "{file}_hybrid_mm10_alignstats.txt"
	log:
		out = config['dirs']['logdir'] + "{file}" + "_process_bam16.log",
		err = config['dirs']['logdir'] + "{file}" + "_process_bam16.err"
	params:
		alignstats = config['tools']['alignstats'],
		vcrome_bed = config['data']['bed']['vcrome_bed']
	threads: 2
	shell:
		"""
		{params.alignstats} -v -i {input.bam} -o {output.realignstats} -t {params.vcrome_bed} -C -W 2> {log.err} 1> {log.out} 
		"""

rule process_bam17:
	input:
		bam = config['dirs']['outdirs']['realign'] + "{file}.tmp.markdups.bam"
	output:
		bam = config['dirs']['outdirs']['realign'] + "{file}.hybrid_hg19_mm10.bam"
	log:
		out = config['dirs']['logdir'] + "{file}" + "_process_bam17.log",
		err = config['dirs']['logdir'] + "{file}" + "_process_bam17.err"
	params:
		samtools = config['tools']['samtools'],
		genomefile_hybrid_index = config['data']['ref']['genomefile_hybrid_index']
	threads: 2
	shell:
		"""
		{params.samtools} view {input.bam} | grep human | grep mouse | \
		{params.samtools} view -bt {params.genomefile_hybrid_index} - > {output.bam} 2> {log.err} 1> {log.out} 
		"""

rule process_bam18:
	input:
		bam = config['dirs']['outdirs']['realign'] + "{file}.hybrid_hg19_mm10.bam"
	output:
		bam = config['dirs']['outdirs']['realign'] + "{file}.hybrid_hg19_mm10.bam.bai"
	log:
		out = config['dirs']['logdir'] + "{file}" + "_process_bam18.log",
		err = config['dirs']['logdir'] + "{file}" + "_process_bam18.err"
	params:
		samtools = config['tools']['samtools']
	threads: 2
	shell:
		"""
		{params.samtools} index -t 8 {input.bam} 2> {log.err} 1> {log.out} 
		"""

rule process_bam19:
	input:
		bam = config['dirs']['outdirs']['realign'] + "{file}.hybrid_hg19_mm10.bam"
	output:
		realignstats = config['dirs']['outdirs']['realignstatsdir'] + "{file}" + "/" + "{file}_hybrid_hg19_mm10_alignstats.txt"
	log:
		out = config['dirs']['logdir'] + "{file}" + "_process_bam19.log",
		err = config['dirs']['logdir'] + "{file}" + "_process_bam19.err"
	params:
		alignstats = config['tools']['alignstats'],
		vcrome_bed = config['data']['bed']['vcrome_bed']
	threads: 2
	shell:
		"""
		{params.alignstats} -v -i {input.bam} -o {output.realignstats} -t {params.vcrome_bed} -C -W 2> {log.err} 1> {log.out} 
		"""

rule gatk_process1:
	input:
		bam = config['dirs']['outdirs']['realign'] + "{file}.hybrid_hg19.bam"
	output:
		intervals = temp(config['dirs']['outdirs']['realign'] + "{file}.GATKrealign.intervals")
	log:
		out = config['dirs']['logdir'] + "{file}" + "_gatk_RealignerTargetCreator.log",
		err = config['dirs']['logdir'] + "{file}" + "_gatk_RealignerTargetCreator.err"
	params:
		java = config['binaries']['java'],
		hg19_fasta = config['data']['ref']['hg19_fasta'],
		gatk_jar = config['tools']['gatk_jar'],
		gatk_ref = config['data']['ref']['gatk_ref']
	threads: 4
	shell:
		"""
		{params.java} -Xmx23g -jar {params.gatk_jar}/GenomeAnalysisTK.jar \
		-T RealignerTargetCreator -R {params.hg19_fasta} -nt 4 \
		-I {input.bam} \
		--downsampling_type NONE \
		-known {params.gatk_ref}/1000G_phase1.indels.b37.vcf.gz \
		-known {params.gatk_ref}/Mills_and_1000G_gold_standard.indels.b37.vcf.gz \
		-o {output.intervals} 2> {log.err} 1> {log.out}
		"""

rule gatk_process2:
	input:
		bam = config['dirs']['outdirs']['realign'] + "{file}.hybrid_hg19.bam",
		intervals = config['dirs']['outdirs']['realign'] + "{file}.GATKrealign.intervals"
	output:
		bam = temp(config['dirs']['outdirs']['realign'] + "{file}.realigned.bam")
	log:
		out = config['dirs']['logdir'] + "{file}" + "_gatk_IndelRealigner.log",
		err = config['dirs']['logdir'] + "{file}" + "_gatk_IndelRealigner.err"
	params:
		java = config['binaries']['java'],
		hg19_fasta = config['data']['ref']['hg19_fasta'],
		gatk_jar = config['tools']['gatk_jar'],
		gatk_ref = config['data']['ref']['gatk_ref']
	threads: 4
	shell:
		"""
		{params.java} -Xmx23g -jar {params.gatk_jar}/GenomeAnalysisTK.jar \
		-T IndelRealigner -R {params.hg19_fasta} \
		-I {input.bam} \
		-targetIntervals {input.intervals} --downsampling_type NONE \
		--consensusDeterminationModel USE_READS \
		-known {params.gatk_ref}/1000G_phase1.indels.b37.vcf.gz \
		-known {params.gatk_ref}/Mills_and_1000G_gold_standard.indels.b37.vcf.gz \
		--maxReadsInMemory 250000 --bam_compression 1 \
		-o {output.bam} 2> {log.err} 1> {log.out}
		"""

rule gatk_process3:
	input:
		bam = config['dirs']['outdirs']['realign'] + "{file}.realigned.bam"
	output:
		grp = temp(config['dirs']['outdirs']['realign'] + "{file}.GATKrecal_data.grp")
	log:
		out = config['dirs']['logdir'] + "{file}" + "_gatk_BaseRecalibrator.log",
		err = config['dirs']['logdir'] + "{file}" + "_gatk_BaseRecalibrator.err"
	params:
		java = config['binaries']['java'],
		hg19_fasta = config['data']['ref']['hg19_fasta'],
		gatk_jar = config['tools']['gatk_jar'],
		gatk_ref = config['data']['ref']['gatk_ref']
	threads: 4
	shell:
		"""
		{params.java} -Xmx23g -jar {params.gatk_jar}/GenomeAnalysisTK.jar \
		-T BaseRecalibrator -I {input.bam} \
		-R {params.hg19_fasta} \
		-o {output.grp} -nct 4 --downsampling_type NONE \
		-knownSites {params.gatk_ref}/dbsnp_142_b37.vcf.gz \
		-knownSites {params.gatk_ref}/1000G_phase1.indels.b37.vcf.gz \
		-knownSites {params.gatk_ref}/Mills_and_1000G_gold_standard.indels.b37.vcf.gz \
		-cov ReadGroupCovariate -cov QualityScoreCovariate \
		-cov CycleCovariate -cov ContextCovariate 2> {log.err} 1> {log.out}
		"""

rule gatk_process4:
	input:
		bam = config['dirs']['outdirs']['realign'] + "{file}.realigned.bam",
		grp = config['dirs']['outdirs']['realign'] + "{file}.GATKrecal_data.grp"
	output:
		bam = config['dirs']['outdirs']['realign'] + "{file}.hybrid_hg19.realigned.recal.bam"
	log:
		out = config['dirs']['logdir'] + "{file}" + "_gatk_PrintReads.log",
		err = config['dirs']['logdir'] + "{file}" + "_gatk_PrintReads.err"
	params:
		tmpdir = config['dirs']['tmpdir'],
		java = config['binaries']['java'],
		hg19_fasta = config['data']['ref']['hg19_fasta'],
		gatk_jar = config['tools']['gatk_jar'],
		gatk_ref = config['data']['ref']['gatk_ref']
	threads: 4
	shell:
		"""
		{params.java} -Xmx23g -Djava.io.tmpdir={params.tmpdir} -jar {params.gatk_jar}/GenomeAnalysisTK.jar \
		-T PrintReads -I {input.bam} \
		-R {params.hg19_fasta} \
		-BQSR {input.grp} --emit_original_quals \
		--downsampling_type NONE -nct 4 --bam_compression 6 \
		-o {output.bam} 2> {log.err} 1> {log.out}
		"""
