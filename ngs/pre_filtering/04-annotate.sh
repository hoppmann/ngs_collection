OUT=04-annotated


# progs
vep="/data/ngs/bin/variant_effect_predictor/variant_effect_predictor.pl"
vt="/data/ngs/bin/vt-master/vt"

#extract number of possible cpu
CPU=$(cat /proc/cpuinfo | grep processor | wc -l)

mkdir -p $OUT

for i in $(ls $1/*.vcf)
do
        echo $i
	filename=$(basename "$i")
        filename="${filename%.*}"        

	# tabix files
	bgzip -c $i > $i.gz
	tabix $i.gz

	
        
        #prepare for anotation-c 
	echo $i
        zless $i \
        | sed 's/ID=AD,Number=./ID=AD,Number=R/' \
        | $vt decompose -s - \
        | $vt normalize -r /data/ngs/resources/bundle_2.8/ucsc.hg19.fasta - > $OUT/$filename.vt.vcf
        
        # annotate using VEP
	$vep \
        -i $OUT/$filename.vt.vcf \
        -o $OUT/$filename-vep.vcf \
        --fork $CPU \
        --cache \
	--cache_version 77 \
        --assembly GRCh37 \
        --force_overwrite \
        --sift b \
        --polyphen b \
        --symbol \
        --numbers \
        --biotype \
        --total_length \
        --vcf \
        --offline \
        --fields Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE
        
        
        
        # annotate using snpEff
	echo "Running SNPEff!"
	java -jar /data/ngs/bin/snpEff/4.2/snpEff.jar \
	GRCh37.75 \
	-classic \
	-formatEff \
	$OUT/$filename-vep.vcf \
	> $OUT/$filename-vep-snpeff.vcf

	# indexing and ziping out file
	echo "Indexing files!"
	bgzip -c $OUT/$filename-vep-snpeff.vcf > $OUT/$filename-vep-snpeff.vcf.gz
	tabix -p vcf $OUT/$filename-vep-snpeff.vcf.gz        
        
        
        
        
        
        
        
        
done

