import os

rule blast_seqs:
	input:
		#seq = os.path.join(DIRS['out_mafft_in'], "seqs_{id}_" + FILE_HASH + ".fna"),
		seq = dynamic(os.path.join(DIRS['out_mafft_in'], "seqs_{id}_" + FILE_HASH + ".fna")),
		#seq = os.path.join(DIRS['out_mafft_in'], "{sample}.fna"),
		db = expand(os.path.join(DIRS['blast'], "Chrom_Ref_Blast_DB.{end}"), end= ["nin", "nsq", "nhr"])
	output:
		blasthit = os.path.join(DIRS['blast'], "blasthit_{id}_" + FILE_HASH + ".txt")
		#blasthit = os.path.join(DIRS['blast'], "blasthit_{sample}.txt")
	params:
		blast_db = os.path.join(DIRS['blast'], "Chrom_Ref_Blast_DB")
	conda:
		"../envs/blast.yaml"
	version:
		"0.1"
	shell:
		'''
		blastn -query {input.seq} -max_target_seqs 1 -dust no -db {params.blast_db} -out {output.blasthit} -outfmt "6 qaccver saccver qlen nident"
		'''
