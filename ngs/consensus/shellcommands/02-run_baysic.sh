putDir="11-baysic"
mkdir $outDir

inDir="06a-different_caller/"

perl /media/anselm/Daten/programs/baysic/baysic.pl \
--statsOutFile combined.stats \
--pvalCutoff 0.8 \
--vcf 01-vcf/gatk.vcf \
--vcf 01-vcf/samtools.vcf \
--vcf 01-vcf/platypus.vcf \
--countsOutFile combined.cts \
--vcfOutFile combined_80perc.vcf
