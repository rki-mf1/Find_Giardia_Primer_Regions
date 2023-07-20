rule run_alistat:
	input:
		algn = os.path.join(DIRS['out_mafft_out'], "algn_{id}_"  + FILE_HASH + ".fna")
	output:
		stats = os.path.join(DIRS['out_alistat'], "alistat_algn_{id}_"  + FILE_HASH + ".txt")
	version:
		"0.1"
	conda:
		"../envs/hmmer3.yaml"
	threads: 1
	shell:
		'''
		esl-alistat {input.algn} > {output.stats}
		'''
