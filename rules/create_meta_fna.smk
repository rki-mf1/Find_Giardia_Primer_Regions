import re
import os
import gzip
import pandas as pd
from collections import defaultdict
from Bio.Seq import Seq
#from Bio.Alphabet import IUPAC

def convert_to_fasta_entry(datalist, contig_seqs, ref=False):
	contig_key = 3 if ref else 4
	start_key = 5 if ref else 7
	end_key = 6 if ref else 8
	strand_key = 9 if ref else 10
	seq_key = -2 if ref else -1

	contig = [x[contig_key] for x in datalist]
	if len(set(contig)) > 1:
		exit("error: contigs not consistent for meta alignment " + str(datalist[0][0]))
	contig = contig[0]

	strand = [x[strand_key] for x in datalist]
	if len(set(strand)) > 1:
		exit("error: strands not consistent for meta alignment " + str(datalist[0][0]))
	strand = strand[0]

	x = set([int(x[start_key]) for x in datalist] + [int(x[end_key]) for x in datalist])
	s = min(x)
	e = max(x)

	if strand == "+":
		seq = contig_seqs[contig][s-1:e]
		header = ">" + contig + " " + str(s) + ":" + str(e)
	else:
		header = ">" + contig + " complement(" + str(s+2) + ":" + str(e+2) + ")"
		#seq = str(Seq(contig_seqs[contig][s+1:e+2], IUPAC.unambiguous_dna).reverse_complement()) ###old biopython
		seq = str(Seq(contig_seqs[contig][s+1:e+2]).reverse_complement())

	return header + "\n" + seq

rule create_meta_fna:
	input:
		matrix=os.path.join(DIRS['out_data'], "filtered_alignment_matrix.tsv.gz")
	output:
		dynamic(os.path.join(DIRS['out_mafft_in'], "seqs_{a}_" + FILE_HASH + ".fna"))
		#directory(os.path.join(DIRS['out_mafft_in']))
	version:
		"0.1"
	run:
		contig_seqs = fna_to_dict(*FILENAMES)

		#meta_algn_id\tsub_algn_id\tsbj\tref_contig\tsbj_contig\tref_start\tref_end\tsbj_start\tsbj_end\tref_strand\tsbj_strand\tref_seq\tsbj_seq\n
		m = 0
		with gzip.open(input.matrix, "rt") as inhandle:
			inhandle.readline()
			rows = []
			prev = None #meta_align_id
			while True:
				line = inhandle.readline().strip()
				if len(line) > 0:
					fields = line.split("\t")
					meta_algn_id = fields[0]
				else:
					meta_algn_id = None

				if prev == meta_algn_id:
					rows.append(fields)

				else:

					#process prev meta alignment
					if rows:

						#unique sbj contig list contig sorted by sbj
						sbj_contigs = set([(x[2], x[4], x[10]) for x in rows])
						sbj_contigs = sorted(sbj_contigs, key = lambda y: y)

						entries = []
						min_fragment_length = []
						for sbj_data in sbj_contigs:
							sbj, sbj_contig, sbj_strand = sbj_data
							data = [x for x in rows if x[4] == sbj_contig and x[10] == sbj_strand]

							#add reference entry
							if not entries:
								entries.append(convert_to_fasta_entry(data, contig_seqs, True))
								min_fragment_length.append(len(entries[-1].split("\n")[-1]))

							#add subject entry
							entries.append(convert_to_fasta_entry(data, contig_seqs))
							min_fragment_length.append(len(entries[-1].split("\n")[-1]))

						if min(min_fragment_length) >= config['minlen']:
							m += 1
							fname = os.path.join(DIRS['out_mafft_in'], "seqs_" + str(m) + "_" + FILE_HASH + ".fna")
							with open(fname, "w+") as outhandle:
								outhandle.write("\n".join(entries))

					#process current line
					if not line:
						break
					prev = meta_algn_id
					rows = [fields]
