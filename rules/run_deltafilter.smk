rule run_deltafilter:
	input:
		delta = os.path.join(DIRS['out_mummer_delta'], "/{subj}_aligned_to_{ref}.delta")
	output:
		delta = os.path.join(DIRS['out_mummer_delta'], "/{subj}_aligned_to_{ref}.filtered.delta")
	conda:
		"../envs/mummer3.yaml"
	version:
		"0.2"
	shell:
		#-1 1-to-1 alignment
		'''
		delta-filter -1 {input.delta} > {output.delta}
		'''
