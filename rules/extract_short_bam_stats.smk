import os

def input_extract_short_bam_stats(wildcards):
	outfiles = {}
	outfiles['template'] = os.path.join(DIRS['out_ash'], "alignment_" + wildcards.id, "ash_template_" + wildcards.id + "_"  + FILE_HASH + ".txt")
	outfiles['bam'] = os.path.join(DIRS['out_bowtie_mapping'], wildcards.name + ".bam")
	outfiles['fasta'] = os.path.join(DIRS['in'], FQ_SAMPLES[wildcards.name][2])
	return outfiles

rule extract_short_bam_stats:
	input:
		unpack(input_extract_short_bam_stats)
	output:
		stat_file = os.path.join(DIRS['out_ash'], "alignment_{id}", "stats_{id}_{name}_"  + FILE_HASH + ".txt")
	conda:
		"../envs/pysamstats.yaml"
	version:
		"0.1"
	params:
		script = os.path.join(workflow.basedir, "scripts", "pysamstats.sh")
	shell:
		'''
			{params.script} {input.template} {wildcards.name} {input.fasta} {input.bam} {output.stat_file}
		'''
