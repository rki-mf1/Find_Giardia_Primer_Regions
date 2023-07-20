rule sort_align_matrix:
	input:
		matrix = os.path.join(DIRS['out_data'], "alignment_matrix.tsv.gz")
	output:
		gz = os.path.join(DIRS['out_data'], "sorted_alignment_matrix.tsv.gz"),
		tmp = temp(os.path.join(DIRS['out_data'], "sorted_alignment_matrix.tsv"))
	version:
		"0.1"
	params:
		tmpdir = DIRS['out_data']
	threads:
		1000
	shell:
		#sort by ref_strand, ref_contig, ref_start
		'''
		set +o pipefail
		zcat {input.matrix} | head -n 1 > {output.tmp}
		zcat {input.matrix} | sed 1,1d | sort -T {params.tmpdir} --parallel {threads} -k 9,9 -k 3,3 -k5,5n >> {output.tmp}
		gzip -k {output.tmp}
		'''
