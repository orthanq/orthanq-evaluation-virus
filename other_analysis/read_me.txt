Here, the labmix sample was quantified using virus-vg.
The resulting fasta containing reconstructed contigs was aligned against the reference fasta containing all five lineages using minimap2.

used commands for virus-vg:
1) pear -f fastp-trimmed/pe/SRR961514.1.fastq -r fastp-trimmed/pe/SRR961514.2.fastq -o SRR961514_PEAR_output
2) savage.py -p1 SRR961514_PEAR_output.unassembled.forward.fastq -p2 SRR961514_PEAR_output.unassembled.reverse.fastq --ref hiv_reference_sequence.fasta -o SRR961514_savage_out --split 15 --min_overlap_len 150
3) python virus-vg/scripts/build_graph_msga.py -f fastp-trimmed/pe/SRR961514.1.fastq -r fastp-trimmed/pe/SRR961514.2.fastq -c SRR961514_savage_out/contigs_stage_c.fasta -vg vg
4) python optimize_strains.py -m 70 -c 140 node_abundance.txt contig_graph.final.gfa