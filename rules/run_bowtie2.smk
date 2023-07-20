import os

def input_run_bowtie2(wildcards):
	r1, r2, ref = FQ_SAMPLES[wildcards.name]
	files = {}
	files['ref'] = os.path.join(DIRS['in'], ref)
	files['index'] = os.path.join(DIRS['out_bowtie_index'], ref + ".1.bt2")
	files['r1'] = os.path.join(DIRS['in_fastq'], r1)
	files['r2'] = os.path.join(DIRS['in_fastq'], r2)
	return files

rule run_bowtie2:
	input:
		unpack(input_run_bowtie2)
	output:
		bam = os.path.join(DIRS['out_bowtie_mapping'], "{name}.bam"),
		stats = os.path.join(DIRS['out_bowtie_mapping'], "{name}.stats.txt")
	log:
		os.path.join(DIRS['log'], "bowtie_samtools.{name}.log")
	conda:
		"../envs/bowtie2.yaml"
	version:
		"0.1"
	threads: 20
	shell:
		'''
		index_base=$(echo {input.index} | sed 's/.1.bt2$//')
		bowtie2 --very-sensitive --gbar 20 -X 15000 -p {threads} --qc-filter -x $index_base -1 {input.r1} -2 {input.r2} 2> {output.stats} | samtools view -Sbh -@ {threads} 2> {log} | samtools sort -@ {threads} -o {output.bam}
		samtools index {output.bam}
		'''
