import re
import os
import gzip
from collections import defaultdict

rule filter_align_matrix:
	input:
		matrix = os.path.join(DIRS['out_data'], "sorted_alignment_matrix.tsv.gz")
	output:
		os.path.join(DIRS['out_data'], "filtered_alignment_matrix.tsv.gz")
	version:
		"0.1"
	run:
		with gzip.open(input.matrix, "rt") as inhandle, gzip.open(output[0], "wt") as outhandle:
			outhandle.write("meta_align_id\t" + inhandle.readline().strip())

			prev_line = (None, None, None) #ref_contig, ref_strand, ref_start
			prev_group = (None, None, None) #ref_contig, ref_strand, ref_start
			prev_algn_ids = set()
			meta_align_id = 0
			sbjs = set()

			while True:
				line = inhandle.readline()
				fields = line.strip().split("\t")
				if line:
					this_line = (fields[2], fields[8], fields[4]) #ref_contig, ref_strand, ref_start
				else:
					this_line = (None, None, None)

				#group lines refereing to same reference position
				if prev_line == this_line:
					group.append(fields)
					sbjs.add(fields[1])

				else:
					#process groups of lines showing that the reference is covered by all subjects exactly once
					if len(sbjs) == SBJ_COUNT and len(group) == SBJ_COUNT:
						#check if new meta alignment
						this_group = prev_line
						algn_ids = set([x[0] for x in group])

						if this_group[0] != prev_group[0] or this_group[1] != prev_group[1] or abs(int(this_group[2])-int(prev_group[2])) > 1 or algn_ids != prev_algn_ids:
							meta_align_id += 1

						#write group
						for g in group:
							outhandle.write("\n" + str(meta_align_id)  + "\t" + "\t".join(g))

						#update group related prev data
						prev_group = this_group
						prev_algn_ids = algn_ids

					if not line:
						break

					#update group line prev data
					sbjs = {fields[1], }
					group = [fields]
					prev_line = this_line
