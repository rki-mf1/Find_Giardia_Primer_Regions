rule run_showaligns:
	input:
		delta = os.path.join(DIRS['out_mummer_delta'], "{sbj}_aligned_to_{ref}.filtered.delta")
	output:
		algn = os.path.join(DIRS['out_mummer_aligns'], "{sbj}_aligned_to_{ref}.algn.gz")
	params:
		o = os.path.join(DIRS['out_mummer_aligns'], "{sbj}_aligned_to_{ref}.algn")
	conda:
		"../envs/mummer3.yaml"
	version:
		"0.1"
	threads: 1
	shell:
		'''
		echo "" > {params.o}
		grep "^>" {input.delta} | cut -c2- | sed 's/ [0-9]\+ [0-9]\+$//' | sort | uniq | xargs -I % sh -c "show-aligns -r {input.delta} % >> {params.o}"
		gzip {params.o}
		'''