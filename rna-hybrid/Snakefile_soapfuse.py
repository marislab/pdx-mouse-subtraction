from util.varsub import varsub

configfile: "config.yaml"
varsub(config)

shell.prefix("source ~/.bash_profile")

with open(config['samples']) as f:
	SAMPLES = f.read().splitlines()
	print(SAMPLES)

rule all:
	input:
		expand(config['dirs']['outdirs']['soapfusedir'] + "{file}" + "/" + "{file}_sampleList.txt", file = SAMPLES),
		
# soapfuse
rule run_soapfuse:
	input:
		fq1 = config['dirs']['outdirs']['bamhybrid'] + "{file}" + "/" + "{file}_hybrid_hg19_1_sequence.fq",
		fq2 = config['dirs']['outdirs']['bamhybrid'] + "{file}" + "/" + "{file}_hybrid_hg19_2_sequence.fq"
	output:
		fq1 = config['dirs']['outdirs']['soapfusedir'] + "{file}" + "/" + "fastq" + "/" + "{file}_1.fq",
		fq2 = config['dirs']['outdirs']['soapfusedir'] + "{file}" + "/" + "fastq" + "/" + "{file}_2.fq",
		slist = config['dirs']['outdirs']['soapfusedir'] + "{file}" + "/" + "{file}_sampleList.txt",
		fusions = config['dirs']['outdirs']['soapfusedir'] + "{file}" + "/" + "{file}.final.Fusion.specific.for.genes"
	log:
		out = config['dirs']['logdir'] + "{file}" + "_soapfuse.log",
		err = config['dirs']['logdir'] + "{file}" + "_soapfuse.err"
	params:
		outdir = config['dirs']['outdirs']['soapfusedir'] + "{file}" + "/" + "fastq" + "/",
		sample = "{file}",
		soapfusedir = config['dirs']['outdirs']['soapfusedir'],
		library = "fastq",
		perl = config['binaries']['perl'],
		soapfuse = config['tools']['soapfuse']
	threads: 4
	shell:
		"""
		mkdir -p {params.soapfusedir}
		mkdir -p {params.outdir}

		ln -s {input.fq1} {output.fq1}
		ln -s {input.fq2} {output.fq2}
		echo "{params.sample} {params.library} {params.sample} 101" > {output.slist}

		{params.perl} {params.soapfuse}/SOAPfuse-RUN.pl \
		-c {params.soapfuse}/config/config.txt \
		-fd {params.soapfusedir} \
		-l {output.slist} \
		-o {params.soapfusedir}  2> {log.err} 1> {log.out}
		"""

# defuse
# rule run_defuse:
# 	input:
# 		fq1 = config['dirs']['outdirs']['bamhybrid'] + "{file}" + "/" + "{file}_hybrid_hg19_1_sequence.fq",
# 		fq2 = config['dirs']['outdirs']['bamhybrid'] + "{file}" + "/" + "{file}_hybrid_hg19_2_sequence.fq"
# 	output:
# 		results = config['dirs']['outdirs']['defusedir'] + "{file}" + "/" + "results.filtered.tsv"
# 	params:
# 		defusedir = config['dirs']['outdirs']['defusedir'] + "{file}" + "/",
# 		defuse = config['tools']['defuse']
# 	threads: 4
# 	shell:
# 		"""
# 		mkdir -p {params.defusedir}

# 		cd {params.defusedir} && {params.defuse} {input.fq1} {input.fq2}
# 		"""
