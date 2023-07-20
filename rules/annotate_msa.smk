import os
from collections import defaultdict

rule annotate_msa:
	input:
		tsv = os.path.join(DIRS['out_base'], "ranking_k" + str(config['sliding_window_diff']) + "_ash" + str(config['min_ash']) + "_reads" + str(config['min_number_reads']) + "_" + FILE_HASH + ".txt")
	output:
		dynamic(os.path.join(DIRS['out_final'], "tmp", "annotated_alignment_{m_id}_k" + str(config['sliding_window_diff']) + "_ash" + str(config['min_ash']) + "_reads" + str(config['min_number_reads']) + "_" + FILE_HASH + ".fasta"))
	version:
		"0.1"
	run:
		pois = {}
		#read input file and skip first line
		with open(input.tsv, "r") as inhandle:
			inhandle.readline()
			#loop over lines
			for line in inhandle:
				fields = line.strip().split("\t")
				if len(fields) == 0:
					continue
				#fields columns only one important	
				alignment_fname = fields[1]
				aln_data = alignment_fname.split("_")
				m_id = aln_data[1]
				pois[alignment_fname] = []########################
				#open detailed alignment statistics
				with open(os.path.join(DIRS['out_data'], alignment_fname), "r") as handle:
					handle.readline()
					for lines in handle:
						details = lines.strip().split("\t")
						if len(details) == 0:
							continue
						#last element in line is POI T or F
						if details[-1] == "T":
							pois[alignment_fname].append("Y")
						else:
							pois[alignment_fname].append("-")
					# write data
					outname = os.path.join(DIRS['out_final'], "tmp", "annotated_alignment_" + str(m_id) + "_k" + str(config['sliding_window_diff']) + "_ash" + str(config['min_ash']) + "_reads" + str(config['min_number_reads']) + "_" + FILE_HASH + ".fasta")
					with open(outname, "w") as outhandle:
						## write header (blast dabases are shown in same order as listed in REFERENCES)
						outhandle.write(">POIs\n" + "".join(pois[alignment_fname]) + "\n")

