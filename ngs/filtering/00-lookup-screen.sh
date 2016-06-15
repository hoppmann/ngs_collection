#!/bin/bash
#General definitions
gemini="/data/ngs/bin/gemini/bin/gemini"
dbDir=$1

# list of genes for lookup
genes=($2)

#columns="gene, chrom, start, ref, alt, depth, qual, is_conserved, gerp_bp_score, aaf_1kg_eur, impact_severity, impact, exon, is_lof, codon_change, aa_change, pfam_domain, polyphen_pred, polyphen_score, sift_pred, sift_score,  (gt_depths).(*), (gts).(*), (gt_types).(*)"

#set column variable via column bash script
. /h/hoppmann/scripts/ngs/filtering/00-columns.sh


for gene in ${genes[@]}
do
	echo "#### $gene ####" 
	echo "" 
	
	$gemini query -q "select $columns from variants where gene ='$gene' " \
	--header \
	$dbDir  
	echo ""
	echo ""
done

