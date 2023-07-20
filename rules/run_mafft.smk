import gzip

rule run_mafft:
	input:
		fna = (os.path.join(DIRS['out_mafft_in'], "seqs_{id}_" + FILE_HASH + ".fna"))
	output:
		mafft = (os.path.join(DIRS['out_mafft_out'], "algn_{id}_" + FILE_HASH + ".fna"))
	version:
		"0.1"
	conda:
		"../envs/mafft.yaml"		
	threads: 10
	shell:
		"mafft --maxiterate 1000 --globalpair --thread {threads} {input.fna} > {output.mafft}"
