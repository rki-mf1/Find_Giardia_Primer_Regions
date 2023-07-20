import os

rule make_blastdb:
	input:
		reference = lambda wildcards: config["references"][wildcards.ref]
	output:
		expand(os.path.join(DIRS['blast'], "Chrom_{{ref}}_Blast_DB.{end}"), end= ["nin", "nsq", "nhr"])
	conda:
		"../envs/blast.yaml"
	params:
		name = os.path.join(DIRS['blast'], "Chrom_{ref}_Blast_DB")
	version:
		"0.1"
	threads: 10
	shell:
		'''
		makeblastdb -dbtype nucl -in {input.reference} -out {params.name}
		'''
