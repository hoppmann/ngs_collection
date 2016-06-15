OUT=04-annotated

#extract number of possible cpu
CPU=$(cat /proc/cpuinfo | grep processor | wc -l)

mkdir -p $OUT

for i in $(ls $1/*.vcf)
do
        echo $i
	filename=$(basename "$i")
        filename="${filename%.*}"        
        
        #prepare for anotation

	echo $i
        zless $i | sed 's/ID=AD,Number=./ID=AD,Number=R/' | /data/ngs/bin/vt-master/vt decompose -s - | /data/ngs/bin/vt-master/vt normalize -r /data/ngs/programs/bundle_2.8/ucsc.hg19.fasta - > $OUT/$filename.vt.vcf
        
        
	/data/ngs/bin/variant_effect_predictor/variant_effect_predictor.pl \
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
done

