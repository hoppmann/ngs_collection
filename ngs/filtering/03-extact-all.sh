#columns="gene, chrom, start, ref, alt, depth, is_conserved, gerp_bp_score, aaf_1kg_eur, impact_severity, impact, exon, is_lof, codon_change, aa_change, pfam_domain, polyphen_pred, polyphen_score, sift_pred, sift_score, (gt_depths).(*), (gts).(*), (gt_types).(*)"

#get columns from column script
. /h/hoppmann/scripts/ngs/filtering/00-columns.sh

gemini query -q "select $columns from variants" --header $1 > $2
