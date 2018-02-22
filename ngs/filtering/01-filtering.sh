
#prepare out folder
outFolder="trioAnalysis"
if [ ! -d $outFolder ]
then 
mkdir $outFolder
fi

. /h/hoppmann/scripts/ngs/filtering/00-columns.sh

## gemini path
#gemini="/data/ngs/bin/gemini/anaconda/bin/gemini"

#define database and extract name
dbDir="/data/ngs/rohdaten/kinderklinik/solid/2015-02-16_run6/out/laterAnalyis/02-databases/510.db"
bn=$(basename $dbDir)
name=${bn%.*}

# define relatives
patient="923_1"
relative_1="923_2"
relative_2="923_3"


##################################
######## Recessive filter ########
##################################


#define out dir
OUT="$outFolder/$name-REC-HIGH.xls"
echo "$OUT"

# query database
gemini query -q "select $columns from variants \
where (impact_severity == 'HIGH' or impact_severity is null) \
and (aaf_gnomad_afr < 0.01 or aaf_gnomad_afr > 0.99 or aaf_gnomad_afr is null) \
and (aaf_gnomad_eas < 0.01 or aaf_gnomad_eas > 0.99 or aaf_gnomad_eas is null) \
and (aaf_gnomad_nfe < 0.01 or aaf_gnomad_nfe > 0.99 or aaf_gnomad_nfe is null) \
and (aaf_gnomad_sas < 0.01 or aaf_gnomad_sas > 0.99 or aaf_gnomad_sas is null) \
" \
--gt-filter "gt_types.$patient == 3 and gt_types.$relative_1 == 1 and gt_types.$relative_2 == 1" \
--header \
$dbDir > $OUT




# define out dir
OUT="$outFolder/$name-REC-MED.xls"
echo "$OUT"

# query database
gemini query -q "select $columns from variants \
where (impact_severity == 'MED' or impact_severity == 'HIGH' or impact_severity is null) \
and (aaf_gnomad_afr < 0.01 or aaf_gnomad_afr > 0.99 or aaf_gnomad_afr is null) \
and (aaf_gnomad_eas < 0.01 or aaf_gnomad_eas > 0.99 or aaf_gnomad_eas is null) \
and (aaf_gnomad_nfe < 0.01 or aaf_gnomad_nfe > 0.99 or aaf_gnomad_nfe is null) \
and (aaf_gnomad_sas < 0.01 or aaf_gnomad_sas > 0.99 or aaf_gnomad_sas is null) \
" \
--gt-filter "gt_types.$patient == 3 and gt_types.$relative_1 == 1 and gt_types.$relative_2 == 1" \
--header \
$dbDir > $OUT




#################################
######## Dominant filter ########
#################################


#define out dir
OUT="$outFolder/$name-DOM-HIGH.xls"
echo "$OUT"

# query database
gemini query -q "select $columns from variants \
where (impact_severity == 'HIGH' or impact_severity is null) \
and (aaf_gnomad_afr < 0.01 or aaf_gnomad_afr > 0.99 or aaf_gnomad_afr is null) \
and (aaf_gnomad_eas < 0.01 or aaf_gnomad_eas > 0.99 or aaf_gnomad_eas is null) \
and (aaf_gnomad_nfe < 0.01 or aaf_gnomad_nfe > 0.99 or aaf_gnomad_nfe is null) \
and (aaf_gnomad_sas < 0.01 or aaf_gnomad_sas > 0.99 or aaf_gnomad_sas is null) \
" \
--gt-filter "gt_types.$patient == 1 and (gt_types.$relative_1 == 0 or gt_types.$relative_1 == 2) and (gt_types.$relative_2 == 0 or gt_types.$relative_2 == 2)" \
--header \
$dbDir > $OUT



##define out dir
OUT="$outFolder/$name-DOM-MED.xls"
echo "$OUT"

# query database
gemini query -q "select $columns from variants \
where (impact_severity == 'MED' or impact_severity == 'HIGH' or impact_severity is null) \
and (aaf_gnomad_afr < 0.01 or aaf_gnomad_afr > 0.99 or aaf_gnomad_afr is null) \
and (aaf_gnomad_eas < 0.01 or aaf_gnomad_eas > 0.99 or aaf_gnomad_eas is null) \
and (aaf_gnomad_nfe < 0.01 or aaf_gnomad_nfe > 0.99 or aaf_gnomad_nfe is null) \
and (aaf_gnomad_sas < 0.01 or aaf_gnomad_sas > 0.99 or aaf_gnomad_sas is null) \
" \
--gt-filter "gt_types.$patient == 1 and (gt_types.$relative_1 == 0 or gt_types.$relative_1 == 2) and (gt_types.$relative_2 == 0 or gt_types.$relative_2 == 2)" \
--header \
$dbDir > $OUT




