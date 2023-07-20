rule run_nucmer:
	input:
		ref = os.path.join(DIRS['in'], "{ref}.fasta"),
		fasta = os.path.join(DIRS['in'], "{subj}.fasta")
	output:
		delta = os.path.join(DIRS['out_mummer_delta'], "{subj}_aligned_to_{ref}.delta")
	params:
		o = os.path.join(DIRS['out_mummer_delta'], "{subj}_aligned_to_{ref}")
	conda:
		"../envs/mummer3.yaml"
	version:
		"0.1"
	threads: 1
	shell:
		'''
		nucmer -p {params.o} {input.ref} {input.fasta}
		'''
