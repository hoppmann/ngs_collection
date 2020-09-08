
SCRIPT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"
pipelinePath=`dirname $SCRIPT`



patID=$1
panel=$2
resultDir=$3
# bamFolder=$4
bamFile="04-baseRecalibration/$patID.bam"

$pipelinePath/panelCoverage/00-single_bam.sh $bamFile $panel $resultDir
