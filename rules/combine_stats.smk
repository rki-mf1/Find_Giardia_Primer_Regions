import os
import re

rule combine_stats:
	input:
		mafft = os.path.join(DIRS['out_mafft_out'], "algn_{id}_"  + FILE_HASH + ".fna"),
		ash_stats = expand(os.path.join(DIRS['out_ash'], "alignment_{{id}}", "stats_{{id}}_{name}_"  + FILE_HASH + ".txt"), name=[x for x in ALL_NAMES if x not in NAMES_WITHOUT_FQ]),
		diff_stat = os.path.join(DIRS['out_aligndiff'], "diff_algn_{id}_k" + str(config['sliding_window_diff']) + "_" + FILE_HASH + ".txt"),
		ash_template = os.path.join(DIRS['out_ash'], "alignment_{id}", "ash_template_{id}_" + FILE_HASH + ".txt")
	output:
		tsv = os.path.join(DIRS['out_data'], "alignment_{id}_k" + str(config['sliding_window_diff']) + "_ash" + str(config['min_ash']) + "_reads" + str(config['min_number_reads']) + "_" + FILE_HASH + ".txt")
	version:
		"0.1"
	threads: 1
	run:
		#processing template to store coding strands
		strands = {}
		seq_coords = {}
		with open(input.ash_template, "r") as handle:
			for line in handle:
				fields = line.strip().split("\t")
				if len(fields) == 0:
					continue
				sample = get_iso(fields[0])
				strands[sample] = fields[1]		
				seq_coords[sample] = (int(fields[2]), int(fields[3]))
				
		#MSA differentiating regions
		diffs = []
		with open(input.diff_stat, "r") as handle:
			handle.readline()
			for line in handle:
				line = line.strip()
				if len(line) == 0:
					continue
				diffs.append(line.split("\t")[-1])
				
		if "T" not in diffs:
			open(output.tsv, "w").close()
		
		else:				
			#process msa
			msa = {}
			msa_len = None
			msa_seqs = {}
			for id, seq in fna_to_dict(input.mafft).items():
				sample = get_iso(id)
				msa[sample] = seq if strands[sample] == "+" else seq[::-1]
				msa_seqs[sample] = seq
				if not msa_len:
					msa_len = len(seq)
			
			
			#register gaps in MSA
			gaps = {}
			for sample, seq in msa.items():
				gaps[sample] = set([x[0] for x in enumerate(seq) if x[1] == "-"])
				
			#processing ash_stats
			ashs = {}
			reads = {}
			for fname in input.ash_stats:
				sample = None
				msa_pos = -1
				n = 0
				with open(fname, "r") as handle:
					handle.readline()
					for line in handle:
						fields = line.strip().split("\t")
						if len(fields) == 0:
							continue
						n += 1
						msa_pos += 1
						if not sample:
							sample = get_iso(fields[0])
							ashs[sample] = []
							reads[sample] = []
						while msa_pos in gaps[sample]:
							ashs[sample].append("-")
							reads[sample].append("-")
							msa_pos += 1
						ash = 0 if int(fields[3]) == 0 else int(max(fields[5:8]))/int(fields[3])
						ashs[sample].append(ash)
						reads[sample].append(fields[3])
						
				##consider gaps at sequence
				while msa_pos+1 in gaps[sample]:
					ashs[sample].append("-")
					reads[sample].append("-")
					msa_pos += 1
				
				if msa_len != len(ashs[sample]):
					exit("error: msa length differs from ash values + gaps")
					
				if strands[sample] == "-":
					ashs[sample] = ashs[sample][::-1]
					reads[sample] = reads[sample][::-1]	
					
			#combining
			rng = range(msa_len)
					
			
			##msa position
			rows = [[x+1] for x in rng]

			##coords
			for sample in ALL_NAMES:
				if strands[sample] == "+":
					seq_pos = seq_coords[sample][0]
					step = 1
				else:
					seq_pos = seq_coords[sample][1]
					step = -1	
					
				this_gaps = set([x[0] for x in enumerate(msa_seqs[sample]) if x[1] == "-"]) 
				for msa_pos in rng:
					if msa_pos in this_gaps:
						rows[msa_pos].append("-")
					else:
						rows[msa_pos].append(seq_pos)
					seq_pos += step
			
			##MSA bases
			for sample in ALL_NAMES:
				
				for msa_pos in rng:		
					rows[msa_pos].append(msa_seqs[sample][msa_pos])		
					
			##MSA number of reads
			for sample in [x for x in ALL_NAMES if x not in NAMES_WITHOUT_FQ]:
				for msa_pos in rng:		
					rows[msa_pos].append(reads[sample][msa_pos])	

			##ASH
			for sample in [x for x in ALL_NAMES if x not in NAMES_WITHOUT_FQ]:
				for msa_pos in rng:
					rows[msa_pos].append(ashs[sample][msa_pos])		
					
			##diff region
			for msa_pos in rng:		
				rows[msa_pos].append(diffs[msa_pos])
				
			##poi
			for msa_pos in rng:	
				this_ashs = [ashs[x][msa_pos] for x in FQ_SAMPLES if x in ashs and ashs[x][msa_pos] != "-"]
				this_reads = [reads[x][msa_pos] for x in FQ_SAMPLES if x in reads and reads[x][msa_pos] != "-"]
				if len(this_ashs) == 0:
					rows[msa_pos].append("-")	
				elif diffs[msa_pos] == "T" and min(this_ashs) >= config['min_ash'] and (int(min(this_reads))) >= int(config['min_number_reads']):
					rows[msa_pos].append("T")				
				else:
					rows[msa_pos].append("F")				

			##writing output
			with open(output.tsv, "w") as handle:
				#headline
				handle.write("position - MSA\t")
				handle.write("\t".join(["position - " + x for x in ALL_NAMES]) + "\t")
				handle.write("\t".join(["base - " + x for x in ALL_NAMES]) + "\t")
				handle.write("\t".join(["#reads - " + x for x in ALL_NAMES if x not in NAMES_WITHOUT_FQ]) + "\t")
				handle.write("\t".join(["ASH - " + x for x in ALL_NAMES if x not in NAMES_WITHOUT_FQ]) + "\t")
				handle.write("differentiating kmer" + "\t")
				handle.write("POI" + "\n")
				#data rows
				for row in rows:
					row = [str(x) for x in row]
					handle.write("\t".join(row) + "\n")
					
			

		