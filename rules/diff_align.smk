rule diff_align:
	input:
		algn = os.path.join(DIRS['out_mafft_out'], "algn_{id}_"  + FILE_HASH + ".fna")
	output:
		diff = os.path.join(DIRS['out_aligndiff'], "diff_algn_{id}_k" + str(config['sliding_window_diff']) + "_" + FILE_HASH + ".txt")
	version:
		"0.1"
	threads: 1
	run:
		aligned_seqs=list(fna_to_dict(input.algn).values())
		k = int((config['sliding_window_diff']-1)/2)
		l = len(aligned_seqs[0])
		n = len(aligned_seqs)

		with open(output.diff, "w") as outhandle:
			outhandle.write("pos\tdiff")
			for i in range(l):
				j = i+1
				s = max(0, i - k)
				e = min(l, j + k)
				kmers = set([x[s:e] for x in aligned_seqs])
				if len(kmers) == n:
					d = "T"
				else:
					d = "F"
				outhandle.write("\n" + str(j) + "\t" + d)
