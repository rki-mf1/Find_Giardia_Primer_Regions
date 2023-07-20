import os
import pprint
import hashlib
import sys
from Bio import SeqIO

pp = pprint.PrettyPrinter(indent=4)
DEBUG = pp.pprint

if (config["sliding_window_diff"]%2) == 0:
	sys.exit("error: sliding_window_diff must be an odd integer.")

#DIRECTORIES FOR INPUT/OUTPUT
DIRS = {
	"in": config['fasta_dir'],
	"in_fastq": config['fastq_dir'],
	"out_base": config['outdir'],
	"out_algn": os.path.join(config['outdir'], "algn"),
	"out_seq": os.path.join(config['outdir'], "seq"),
	"out_mummer": os.path.join(config['outdir'], "mummer"),
	"out_mummer_delta": os.path.join(config['outdir'], "mummer", "delta"),
	"out_mummer_coords": os.path.join(config['outdir'], "mummer", "coords"),
	"out_mummer_aligns": os.path.join(config['outdir'], "mummer", "aligns"),
	"out_mafft_in": os.path.join(config['outdir'], "mafft", "seqs"),
	"out_mafft_out": os.path.join(config['outdir'], "mafft", "aligns"),
	"out_data": os.path.join(config['outdir'], "data"),
	"out_alistat": os.path.join(config['outdir'], "align_features", "alistat"),
	"out_aligndiff": os.path.join(config['outdir'], "align_features", "diff"),
	"out_bowtie_index": os.path.join(config['outdir'], "bowtie2", "index"),
	"out_bowtie_mapping": os.path.join(config['outdir'], "bowtie2", "mapping"),
	"out_ash": os.path.join(config['outdir'], "ash"),
	"out_final": os.path.join(config['outdir'], "final"),
	"log": os.path.join(config['outdir'], "log"),
	"blast":os.path.join(config['outdir'], "blast")
}


#FUNCTIONS
def read_samplesheet():
	samples = {}
	with open(config['samplesheet']) as handle:
		handle.readline()
		for line in handle:
			if len(line.strip()) == 0:
				continue
			r1, r2, ref, name = line.strip().split("\t")
			
			if name in samples:
				sys.exit("error: " + name + " not unique in sample sheet.")

			samples[name] = (r1, r2, ref)
	return samples

def fna_to_dict(*fnames):
	contigs = {} #key: header, val: sequence
	for fname in fnames:
		for record in SeqIO.parse(fname, "fasta"):
			if record.id in contigs and contigs[record.id] != str(record.seq):
				sys.exit("error: " + record.id + " refers to different sequences")
			contigs[record.id] = str(record.seq)
	return contigs

def get_iso(contig_id):
	return ISO_PATTERN.search(contig_id).group(1)

#CONSTANTS
FILES = [x for x in os.listdir(DIRS['in']) if not x.endswith(".fai")]
FILENAMES = [ os.path.join(DIRS['in'], x) for x in FILES ]
FILE_HASH = hashlib.md5(" ".join(FILES).encode('utf-8')).hexdigest()
REF_FNAME = os.path.join(DIRS['in'], FILES[0])
REF_FILE = os.path.splitext(os.path.basename(FILES[0]))[0]
SBJ_FILES = [ os.path.splitext(os.path.basename(x))[0] for x in FILES[1:] ]
SBJ_COUNT = len(SBJ_FILES)
ISO_PATTERN = re.compile(config['iso_regex'])
FQ_SAMPLES = read_samplesheet()
ALL_NAMES = sorted(FQ_SAMPLES.keys())
NAMES_WITHOUT_FQ = sorted([x for x in FQ_SAMPLES.keys() if FQ_SAMPLES[x][0] == "-"])
REFERENCES = sorted(config["references"].keys())

rule all:
	input:
		dynamic(os.path.join(DIRS['out_final'], "annotated_alignment_{w_id}_k" + str(config['sliding_window_diff']) + "_ash" + str(config['min_ash']) + "_reads" + str(config['min_number_reads']) + "_" + FILE_HASH + ".fasta"))

		

#mummer actions
include: "rules/run_nucmer.smk"
include: "rules/run_deltafilter.smk"
include: "rules/run_showaligns.smk"

#matrix actions
include: "rules/create_algn_matrix.smk"
include: "rules/sort_align_matrix.smk"
include: "rules/filter_align_matrix.smk"

#mafft actions
include: "rules/create_meta_fna.smk"
include: "rules/run_mafft.smk"
include: "rules/run_alistat.smk"
include: "rules/diff_align.smk"

#bowtie2 actions
include: "rules/run_bowtie2_index.smk"
include: "rules/run_bowtie2.smk"

#ash actions
include: "rules/create_ASH_template.smk"
include: "rules/run_bam_index.smk"
include: "rules/extract_short_bam_stats.smk"

#analysis
include: "rules/combine_stats.smk"
include: "rules/rank_msas.smk"

#blast
include: "rules/make_blastdb2.smk"
include: "rules/blast_seqs2.smk"

#annotate alignment
include: "rules/annotate_msa.smk"
include: "rules/combine_msa_anno.smk"