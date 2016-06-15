if [ ! -d out ]
then 
mkdir out
fi

columns="gene, chrom, start, ref, alt, depth, qual, transcript, is_conserved, gerp_bp_score, aaf_1kg_eur, impact_severity, impact, exon, is_lof, codon_change, aa_change, pfam_domain, polyphen_pred, polyphen_score, sift_pred, sift_score, (gt_depths).(*), (gts).(*), (gt_types).(*)"

. /h/hoppmann/scripts/ngs/filtering/00-columns.sh

#define database
dbDir="/data/ngs/rohdaten/kinderklinik/2015-03-05_run7/out/11-databases/213-VEP.db"
gemini="/data/ngs/bin/gemini/bin/gemini"


#define out dir
OUT="out/213-filter-HIGH.xls"

#write header inforamtion
echo "# autosomal recessive inheritance (konsaguin) " > $OUT
echo "# minor allele frequency of max 0.01 or null in 1kgp " >> $OUT
echo "# homozygos in patient, heterozygos in parents" >> $OUT
echo "# severity of HIGH or null" >> $OUT
echo "" >> $OUT

gemini query -q "select $columns from variants \
where (impact_severity == 'HIGH' or impact_severity is null) \
and (aaf_1kg_eur < 0.01 or aaf_1kg_eur is null) " \
--gt-filter "gt_types.213_1 == HOM_ALT and gt_types.213_2 == HET and gt_types.213_3 == HET" \
--header \
$dbDir >> $OUT

#add empty line
echo "" >> $OUT



#define out dir
OUT="out/213-filter-MED.xls"

#write header inforamtion
echo "# autosomal recessive inheritance (konsaguin) " > $OUT
echo "# minor allele frequency of max 0.01 or null in 1kgp " >> $OUT
echo "# homozygos in patient, heterozygos in Parents" >> $OUT
echo "# severity of MED or HIGH or null" >> $OUT
echo "" >> $OUT

gemini query -q "select $columns from variants \
where (impact_severity == 'MED' or impact_severity == 'HIGH' or impact_severity is null) \
and (aaf_1kg_eur < 0.01 or aaf_1kg_eur is null) " \
--gt-filter "gt_types.213_1 == HOM_ALT and gt_types.213_2 == HET and gt_types.213_3 == HET" \
--header \
$dbDir >> $OUT

#add empty line
echo "" >> $OUT



#define out dir
OUT="out/213-filter-dominant-HIGH.xls"

#write header inforamtion
echo "# autosomal dominant inheritance" > $OUT
echo "# minor allele frequency of max 0.01 or null in 1kgp " >> $OUT
echo "# HET in patient, HOM_REF in Parents" >> $OUT
echo "# severity HIGH or null" >> $OUT
echo "" >> $OUT

$gemini query -q "select $columns from variants \
where (impact_severity == 'HIGH' or impact_severity is null) \
and (aaf_1kg_eur < 0.01 or aaf_1kg_eur is null) " \
--gt-filter "gt_types.213_1 == 1 and (gt_types.213_2 == 0 or gt_types.213_2 == 2) and (gt_types.213_3 == 0 or gt_types.213_3 == 2)" \
--header \
$dbDir >> $OUT



#define out dir
OUT="out/213-filter-dominant-MED.xls"

#write header inforamtion
echo "# autosomal dominant inheritance" > $OUT
echo "# minor allele frequency of max 0.01 or null in 1kgp " >> $OUT
echo "# HET in patient, HOM_REF in Parents" >> $OUT
echo "# severity MED, HIGH or null" >> $OUT
echo "" >> $OUT

$gemini query -q "select $columns from variants \
where (impact_severity == 'MED' or impact_severity == 'HIGH' or impact_severity is null) \
and (aaf_1kg_eur < 0.01 or aaf_1kg_eur is null) " \
--gt-filter "gt_types.213_1 == 1 and (gt_types.213_2 == 0 or gt_types.213_2 == 2) and (gt_types.213_3 == 0 or gt_types.213_3 == 2)" \
--header \
$dbDir >> $OUT




