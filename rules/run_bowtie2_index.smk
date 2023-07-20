import os

rule run_bowtie2_index:
	input:
		ref = os.path.join(DIRS['in'], "{reference}")
	output:
		os.path.join(DIRS['out_bowtie_index'], "{reference}.1.bt2"),
		os.path.join(DIRS['out_bowtie_index'], "{reference}.rev.1.bt2")
	params:
		index_base = os.path.join(DIRS['out_bowtie_index'], '{reference}')
	conda:
		"../envs/bowtie2.yaml"
	version:
		"0.1"
	threads: 10
	shell:
		'''
		bowtie2-build --threads {threads} {input.ref} {params.index_base}
		'''
