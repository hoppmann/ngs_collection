
fileIn=$1

file=$(basename "$fileIn")
file="${file%.*}"

fileOut=$file-all.out

echo $fileOut

#get columns from column script
. /h/hoppmann/scripts/ngs/filtering/00-columns.sh




gemini query -q \
"select $columns from variants" \
--header $fileIn > $fileOut
