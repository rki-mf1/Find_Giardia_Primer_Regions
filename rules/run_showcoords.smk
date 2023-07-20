rule run_showcoords:
	input:
		delta = os.path.join(DIRS['out_mummer_delta'], "{subj}_aligned_to_" + REF_FILE + ".filtered.delta")
	output:
		stats = os.path.join(DIRS['out_mummer_coords'], "{subj}_aligned_to_{ref}.filtered_coords.txt")
	conda:
		"../envs/mummer3.yaml"
	version:
		"0.1"
	shell:
		'''
		show-coords -c -T {input.delta} > {output.stats}
		'''
