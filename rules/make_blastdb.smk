import os

rule make_blastdb:
	input:
		ref = config['chr_ref_fasta']
	output:
		expand(os.path.join(DIRS['blast'], "Chrom_Ref_Blast_DB.{end}"), end= ["nin", "nsq", "nhr"])
	conda:
		"../envs/blast.yaml"
	params:
		name = os.path.join(DIRS['blast'], "Chrom_Ref_Blast_DB")
	version:
		"0.1"
	threads: 10
	shell:
		'''
		makeblastdb -dbtype nucl -in {input.ref} -out {params.name}
		'''
		
