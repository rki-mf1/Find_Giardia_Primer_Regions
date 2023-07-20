import os
import re

rule create_ASH_template:
	input:
		algn = os.path.join(DIRS['out_mafft_out'], "algn_{id}_"  + FILE_HASH + ".fna")
	output:
		template = os.path.join(DIRS['out_ash'], "alignment_{id}", "ash_template_{id}_"  + FILE_HASH + ".txt")
	version:
		"0.1"
	run:
		coords_pattern = re.compile("^>([^ ]+) (complement\()?([0-9]+):([0-9]+)")

		with open(input.algn, "r") as inhandle, open(output.template, "w") as outhandle:
			for line in inhandle:
				match = coords_pattern.search(line)
				if match:
					contig_name = match.group(1)
					strand = "-" if match.group(2) else "+"
					start = match.group(3)
					end = match.group(4)
					
					outhandle.write(contig_name + "\t" + strand + "\t" + start + "\t" + end + "\n")
