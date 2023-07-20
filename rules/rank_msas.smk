import os
from collections import defaultdict

rule rank_msas:
	input:
		ash_stats = dynamic(os.path.join(DIRS['out_data'], "alignment_{m}_k" + str(config['sliding_window_diff']) + "_ash" + str(config['min_ash']) + "_reads" + str(config['min_number_reads']) + "_" + FILE_HASH + ".txt")),
		blasthits = dynamic(expand(os.path.join(DIRS['blast'], "{ref}", "{ref}_blasthit_{{m}}_" + FILE_HASH + ".txt"), ref=REFERENCES))
	output:
		tsv = os.path.join(DIRS['out_base'], "ranking_k" + str(config['sliding_window_diff']) + "_ash" + str(config['min_ash']) + "_reads" + str(config['min_number_reads']) + "_" + FILE_HASH + ".txt")
	version:
		"0.1"
	threads: 1
	run:
		pois = {}
		lens = {}
		chromosomes = {}
		
		# process each alignment
		## check ash
		### loop over ash files
		for fname in input.ash_stats:
			alignment_id = os.path.basename(fname).split("_")[1]
			pois[fname] = 0
			lens[fname] = 0
			#ash data
			with open(fname, "r") as handle:
				handle.readline()
				for line in handle:
					line = line.strip()
					if len(line) == 0:
						continue
					lens[fname] += 1
					fields = line.split("\t")
					if fields[-1] == "T":
						pois[fname] += 1
						
			## check blast hits for the respective aligment considering all blast dbs
			### loop over blast dbs
			chromosomes[fname] = []
			for blast_db in REFERENCES:
				#check %identity (store in list), number of hit occurence (store in variable), coords of best hit (highest %identity)
				blasthit_file = os.path.join(DIRS['blast'], blast_db, blast_db + "_blasthit_" + alignment_id + "_" + FILE_HASH + ".txt")
				identities = defaultdict(set)
				n = defaultdict(int) 
				ident = -1
				coords = None
				with open(blasthit_file, "r") as handle:
					for line in handle:
						if len(line.strip()) == 0:
							continue
						fields = line.strip().split("\t")
						qcov = abs(int(fields[6]) - int(fields[7]))/int(fields[2])*100
						if qcov >= 50:
							this_ident = round(int(fields[3])/int(fields[2])*100, 2)
							identities[fields[1]].add(this_ident)
							n[fields[1]] += 1
							if this_ident > ident:
								coords = (fields[4], fields[5])
								ident = this_ident
				
				# store data in list in same order as blast database are listed in REFRENCE
				if len(identities) > 0:
					chrom_data = { x:  (n[x], min(identities[x]), max(identities[x])) for x in identities }
					data = []
					for chrom in sorted(chrom_data.keys()):
						data.append(chrom + "(N: " + str(chrom_data[chrom][0]) + "; min_ident: " + str(chrom_data[chrom][1]) + "%; max_ident: " + str(chrom_data[chrom][2]) + "%; max_ident_coords: " + coords[0] + "-" + coords[1] + ")" )
					chromosomes[fname].append(", ".join(data))
				else:
					chromosomes[fname].append("no blast hit")
		
		# write data
		with open(output.tsv, "w") as handle:
			
			## write header (blast dabases are shown in same order as listed in REFERENCES)
			handle.write("#\tmsa_file\tmsa_length\tpoi_count\t" + "\t".join(REFERENCES) + "\n")
			i = 0
			
			#write alignments in descending order of number of pois (blast results are shown in the same order as the repective database is listed in REFERENCES)
			for key, value in sorted(pois.items(), key=lambda item: item[1], reverse=True):
				if lens[key] > 0 and pois[key] > 0:
					i += 1
					handle.write(str(i) + "\t" + os.path.basename(key) + "\t" + str(lens[key]) + "\t" + str(pois[key]) + "\t" + "\t".join(chromosomes[key]) + "\n")
