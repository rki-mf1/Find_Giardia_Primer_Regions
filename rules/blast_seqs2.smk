import os

rule blast_seqs:
	input:
		seq = os.path.join(DIRS['out_mafft_in'], "seqs_{id}_" + FILE_HASH + ".fna"),
		db = expand(os.path.join(DIRS['blast'], "Chrom_{{ref}}_Blast_DB.{end}"), end= ["nin", "nsq", "nhr"])
	output:
		blasthit = os.path.join(DIRS['blast'], "{ref}", "{ref}_blasthit_{id}_" + FILE_HASH + ".txt")
	params:
		blast_db = os.path.join(DIRS['blast'], "Chrom_{ref}_Blast_DB")
	conda:
		"../envs/blast.yaml"
	version:
		"0.1"
	shell:
		'''
		blastn -query {input.seq} -max_target_seqs 1 -dust no -db {params.blast_db} -out {output.blasthit} -outfmt "6 qaccver saccver qlen nident sstart send qstart qend"
		'''
