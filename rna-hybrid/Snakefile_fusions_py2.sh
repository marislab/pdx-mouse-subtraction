from util.varsub import varsub

configfile: "config_py2.yaml"
varsub(config)

shell.prefix("source ~/.bash_profile")

with open(config['samples']) as f:
	SAMPLES = f.read().splitlines()
	print(SAMPLES)

rule all:
	input:
		expand(config['dirs']['outdirs']['bamhybrid'] + "{file}" + "/" + "{file}_hybrid_hg19_{rep}_sequence.txt.bz2", file = SAMPLES, rep = [1,2]),
		expand(config['dirs']['outdirs']['starfusiondir'] + "{file}" + "/" + "star-fusion.fusion_predictions.abridged.annotated.coding_effect.tsv", file = SAMPLES),
		expand(config['dirs']['outdirs']['defusedir'] + "{file}" + "/" + "results.filtered.tsv", file = SAMPLES)


rule decompress_files:
	input:
		fq = config['dirs']['outdirs']['bamhybrid'] + "{file}" + "/" + "{file}_hybrid_hg19_{rep}_sequence.txt.bz2",
	output:
		fq = temp(config['dirs']['outdirs']['bamhybrid'] + "{file}" + "/" + "{file}_hybrid_hg19_{rep}_sequence.fq")
	params:
		bunzip2 = config['binaries']['bunzip2']
	threads: 2
	shell:
		"""
		{params.bunzip2} -c {input.fq} > {output.fq}
		"""

rule star_fusion:
	input:
		fq1 = config['dirs']['outdirs']['bamhybrid'] + "{file}" + "/" + "{file}_hybrid_hg19_1_sequence.fq",
		fq2 = config['dirs']['outdirs']['bamhybrid'] + "{file}" + "/" + "{file}_hybrid_hg19_2_sequence.fq"
	output:
		starfusionout = config['dirs']['outdirs']['starfusiondir'] + "{file}" + "/" + "star-fusion.fusion_predictions.abridged.annotated.coding_effect.tsv"
	log:
		out = config['dirs']['logdir'] + "{file}" + "_star_fusion.log",
		err = config['dirs']['logdir'] + "{file}" + "_star_fusion.err"
	params:
		starfusion = config['tools']['starfusion'],
		stargenomedir = config['data']['ref']['stargenomedir'],
		starfusiondir = config['dirs']['outdirs']['starfusiondir'] + "{file}" + "/",
		trinity = config['tools']['trinity']
	threads: 6
	shell:
		"""
		export TRINITY_HOME={params.trinity}
		echo $TRINITY_HOME

		{params.starfusion} \
		--CPU 4 \
		--genome_lib_dir {params.stargenomedir} \
		--left_fq {input.fq1} \
		--right_fq {input.fq2} \
		--output_dir {params.starfusiondir} \
		--FusionInspector inspect \
		--examine_coding_effect \
		--denovo_reconstruct \
		--annotate 2> {log.err} 1> {log.out}
		"""

# defuse
rule run_defuse:
	input:
		fq1 = config['dirs']['outdirs']['bamhybrid'] + "{file}" + "/" + "{file}_hybrid_hg19_1_sequence.fq",
		fq2 = config['dirs']['outdirs']['bamhybrid'] + "{file}" + "/" + "{file}_hybrid_hg19_2_sequence.fq"
	output:
		results = config['dirs']['outdirs']['defusedir'] + "{file}" + "/" + "results.filtered.tsv",
		todel = temp(config['dirs']['outdirs']['defusedir'] + "{file}" + "/" + "jobs")
	log:
		out = config['dirs']['logdir'] + "{file}" + "_defuse.log",
		err = config['dirs']['logdir'] + "{file}" + "_defuse.err"
	params:
		tmpdir = config['mytmpdir'],
		defuse = config['tools']['defuse'],
		defuse_refdir = config['data']['ref']['defuse_refdir'],
		defuse_config = config['data']['ref']['defuse_config'],
		defusedir = config['dirs']['outdirs']['defusedir'] + "{file}" + "/",
		submitter_type = 'sge'
	threads: 50
	shell:
		"""
		mkdir -p {params.defusedir}

		{params.defuse} -c {params.defuse_config} \
		-d {params.defuse_refdir} \
		-o {params.defusedir} \
		-1 {input.fq1} -2 {input.fq2} \
		-l {params.tmpdir} \
		-s {params.submitter_type} \
		-p {threads} 2> {log.err} 1> {log.out}
		"""