import os
from collections import defaultdict

rule combine_msa_anno:
	input:
		anotation = os.path.join(DIRS['out_final'], "tmp", "annotated_alignment_{w_id}_k" + str(config['sliding_window_diff']) + "_ash" + str(config['min_ash']) + "_reads" + str(config['min_number_reads']) + "_" + FILE_HASH + ".fasta"),
		msa = os.path.join(DIRS['out_mafft_out'], "algn_{w_id}_" + FILE_HASH + ".fna")
	output:
		ano_msa = os.path.join(DIRS['out_final'], "annotated_alignment_{w_id}_k" + str(config['sliding_window_diff']) + "_ash" + str(config['min_ash']) + "_reads" + str(config['min_number_reads']) + "_" + FILE_HASH + ".fasta")
	version:
		"0.1"
	shell:
		'''
			cat {input.anotation} {input.msa} > {output.ano_msa}
		'''