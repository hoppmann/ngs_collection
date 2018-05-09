#fileIn=$1
#failed=$(cat $fileIn | awk '{ if ($5 < 50 ) print}' | wc -l)
#total=$(cat $fileIn | wc -l)
#echo "$failed/$total" | bc -l


# get file in name
fileIn=$1
fileIn="out.xls"
outDir="04-getCoverage"


# get name for intermediate file
fullFile=$(basename "$fileIn")
fileName=${fullFile%.*}
intermedFile="$outDir/$fileName.mean.xls"


# for each line get mean coverage and add to file
# depth / ( end - start)
cat $fileIn | awk '{ meanDepth=($5 / ( $3-$2 ))} {roundMean=sprintf("%.0f", meanDepth)} { print $0 "\t" roundMean}' > $intermedFile

# get percent of coverage of panel exones with decent (>= 20 coverage)
failed=$(cat $intermedFile | awk '{ if ($6 < 20 ) print}' | wc -l)
total=$(cat $intermedFile | wc -l)
coverage=$(echo "1- $failed/$total" | bc -l)

echo "Coverage = $coverage"
