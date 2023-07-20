import re
import os
import gzip
from collections import defaultdict

rule create_algn_matrix:
	input:
		algns = expand(os.path.join(DIRS['out_mummer_aligns'], "{sbj}_aligned_to_" + REF_FILE + ".algn.gz"), sbj=SBJ_FILES)
	output:
		os.path.join(DIRS['out_data'], "alignment_matrix.tsv.gz")
	version:
		"0.1"
	run:
		contig_pattern = re.compile("^-- Alignments between (.+) and (.+)$")
		coords_pattern = re.compile("[[|] ([+-])|([0-9]+) - ([0-9]+)")
		algn_id = 0

		with gzip.open(output[0], "wt") as outhandle:
			outhandle.write("algn_id\tsbj\tref_contig\tsbj_contig\tref_start\tref_end\tsbj_start\tsbj_end\tref_strand\tsbj_strand\tref_seq\tsbj_seq")

			for fname in input.algns:
				with gzip.open(fname, "rt") as inhandle:
					for line in inhandle:

						if line.startswith("-- Alignments between "):
							match = contig_pattern.search(line.rstrip("\r\n"))
							ref_contig = match.group(1)
							sbj_contig = match.group(2)
							sbj = get_iso(sbj_contig)

						elif line.startswith("-- BEGIN alignment "):
							algn_id += 1
							matches = coords_pattern.findall(line)
							ref_strand = matches[0][0]
							ref_start = int(matches[1][1])-1
							ref_end = int(matches[1][2])-1

							sbj_strand = matches[2][0]
							sbj_start = int(matches[3][1])-1
							sbj_end = int(matches[3][2])-1

							ref_seq = []
							sbj_seq = []

							ref_step = int(ref_strand + "1")
							sbj_step = int(sbj_strand + "1")

						elif len(line) > 0 and line[0].isdigit():
							ref_seq.extend(line.split(" ")[-1].strip())
							line = inhandle.readline()
							sbj_seq.extend(line.split(" ")[-1].strip())

						elif line.startswith("--   END alignment "):
							rows = []
							for ref_nuc, sbj_nuc in zip(ref_seq, sbj_seq):
								if ref_nuc == ".":
									sbj_start += sbj_step
									rows[-1][11] += sbj_nuc
									rows[-1][7] = sbj_start
								elif sbj_nuc == ".":
									ref_start += ref_step
									rows[-1][10] += ref_nuc
									rows[-1][5] = ref_start
								else:
									ref_start += ref_step
									sbj_start += sbj_step
									rows.append([algn_id, sbj, ref_contig, sbj_contig, ref_start, ref_start, sbj_start, sbj_start, ref_strand, sbj_strand, ref_nuc, sbj_nuc])
							for row in rows:
								outhandle.write("\n" + "\t".join([str(x) for x in row]))
