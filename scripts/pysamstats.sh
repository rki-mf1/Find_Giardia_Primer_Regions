#!/bin/bash
templatefile="$1"
name="$2"
fasta="$3"
bam="$4"
out="$5"

template_line=$(grep $name $templatefile)
contig=$(echo $template_line | cut -f 1 -d " ")
start=$(echo $template_line | cut -f 3 -d " ")
end=$(echo $template_line | cut -f 4 -d " ")
end=$((end+1))

pysamstats --type variation -D 10000000 -c "$contig" -s "$start" -e "$end" -u -d --fields "chrom,pos,ref,reads_all,matches,A,T,G,C,N" -f "$fasta" "$bam" --output "$out"