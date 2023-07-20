import os

rule run_bam_index:
	input:
		bam = os.path.join(DIRS['out_bowtie_mapping'], "{name}.bam")
	output:
		bai = os.path.join(DIRS['out_bowtie_mapping'], "{name}.bam.bai"),
	version:
		"0.1"
	conda:
		"../envs/bowtie2.yaml"
	shell:
		'''
		samtools index -b {input.bam}
		'''