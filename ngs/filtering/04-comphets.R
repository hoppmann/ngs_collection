library(combinat)


###################################
######## to adapt each run ########
###################################

#get command line options
options = commandArgs(trailingOnly = TRUE)

optionsText=paste("Options have to be in following order:","name(familyID), idPatient, idFather, idMother, input file.", "For idFather and idMother only the last number after the familyID is needed", sep="\n")

if (length(options) != 5) {
	print (options)
	stop(optionsText)
}


#check if all options set
name = options[1] # familiy id name
idPatient = options[2]
idFather = options[3]
idMother = options[4]
fileIN = options[5]

#print out input options
print (options)


######## cut out due to commandline options ########

#to adapt (has to be family ID)
# name = "30_13"
# 
#  # can be found in last loop or in original file needed is the column number of gt_type.father gt_type.mother
# 
# idPatient = 1
# idMother = 2
# idFather = 3
# idBrother = 4


# fileIN = "consensus-vep-snpEff-all.out"


######## beyond in script again ########

#create directory for output
outDir = "09-filter"
dir.create(file.path(outDir), showWarnings = FALSE)


#read in raw data from gemini
data = read.table (fileIN, h = T, sep = "\t")

#################################
######## comp het search ########
#################################

#define variables
impacts = c("HIGH", "MED")
#get coulumn number of father and mother
gtPatient = match (paste("gt_types.", name, "_", idPatient, sep = ""), colnames (data))
gtOne = match (paste("gt_types.", name, "_", idFather, sep = ""), colnames (data))
gtTwo = match (paste("gt_types.", name, "_", idMother, sep = ""), names(data))
# gtBrother = match (paste("gt_types.", name, "_", idBrother, sep = ""), names(data))

for (impact in impacts) {
	outName = paste(name, "_", impact, sep = "")
	print (outName)
	
	#transform colum from factor to numeric
	data$aaf_1kg_eur = as.numeric(as.character(data$aaf_1kg_eur))
	
	#create columnName for later subsetting dependent on that column
	columnName = paste("gt_types.", name, "_", idPatient, sep = "")
	
	#depending on iteration step extract priliminary set of potential genes
	if (impact == "MED" ){
		
		#create subsets with impactseverity HIGH or MED & aaf_1kgp < 0.01 or NONE & HET
		potentialGenes = data[(data$impact_severity == "HIGH" | data$impact_severity == "MED") & data[[columnName]] == 1  & (is.na(data$aaf_1kg_eur) | data$aaf_1kg_eur < 0.01 | data$aaf_1kg_eur > 0.99), ]
		
	} else if (impact == "HIGH") {
		
		#create subsets with impactseverity HIGH & aaf_1kgp < 0.01 or NONE & HET
		potentialGenes = data[(data$impact_severity == "HIGH" ) & data[[columnName]] == 1  & (is.na(data$aaf_1kg_eur) | data$aaf_1kg_eur < 0.01 | data$aaf_1kg_eur > 0.99), ]
	}
	
	#get list of unique genes
	genes = c(as.character(unique(potentialGenes$gene)))
	
	#extract from list genes with more then one mutation
	tmp = data.frame()
	for (i in 1:length(genes)) {
		# for (i in 1:1) {
		currentList = subset (potentialGenes, potentialGenes$gene == genes[i])
		if (nrow(currentList) > 1) {
			tmp = rbind(tmp, currentList)
		}
	}
	
	#set new set of potential genes and get list of individual genes
	potentialGenes = tmp
	genes = unique(potentialGenes$gene)
	
	
	#### extract all mutation pairs not inherited from one single parent
	
	#init variables
	tmp = data.frame()
	cont = 1
	
	for (i in 1:length(genes)) {
		currentList = subset (potentialGenes, potentialGenes$gene == genes[i])
		
		if (length(currentList$end) < 2){
			next
		}
		
		###get all possible pairs of compHets
		combinations = data.frame(combn(currentList$end, 2))
		
		#foreach combination check if patient is identical with both parents
		for (j in 1:length(combinations[1,])) {
			currentPair = t(currentList[currentList$end %in% combinations[,j],])
			
			# if Patient and a single parents share both positions exclude from list
			if (currentPair[gtOne,1] == currentPair[gtPatient,1] && currentPair[gtOne,2] == currentPair[gtPatient,2]){
			} else if (currentPair[gtTwo,1] == currentPair[gtPatient,1] && currentPair[gtTwo,2] == currentPair[gtPatient,2]) {
				
				# 	# check that brother has same mutation as patient
				# } else if (currentPair[gtBrother,1] != currentPair[gtPatient, 1] | currentPair[gtBrother,2] != currentPair[gtPatient, 2]){
			} else {
				
				#check that Parents are not all 0 => at leaste one must be inherited
				if (currentPair[gtOne,1] == 0 & currentPair[gtOne,2] == 0 & currentPair[gtTwo,1] == 0 &currentPair[gtTwo,2] == 0){
					
				} else {
					
					#if not inherited from one parent, give combination a uniqe number and save combination 
					currentPair = data.frame(t(currentPair))
					
					#recover column names for later replacement in tmp
					columnNames = colnames(currentPair)
					columnNames = append(columnNames, "PairNumber")
					
					currentPair$PairNumber=cont
					cont = cont + 1 
					tmp = rbind(tmp, currentPair)
					colnames(tmp) = c(columnNames)
				}
			}
		}
	}
	
	
	#save in corresponding variable
	potentialGenePairs = tmp
	
	
	#remove unneccesary variables
	rm (currentPair, combinations, currentList, tmp, i, j, genes, cont, columnNames)
	
	#get uniqe genes and unique positions within these genes for second output file
	genes = unique (potentialGenePairs$gene)
	uniqueStart = data.frame()
	
	if (length(genes ) > 0) {
		for (i in 1:length(genes)){
			current = subset(potentialGenePairs, potentialGenePairs$gene == genes[i])
			uniqueStart = rbind (uniqueStart, current[!duplicated(current$end),])
		}
	}
	
	
	#remove unused variables 
	rm (current, potentialGenes, genes, i)
	
	#find out number of different genes found an print it to screen
	print (paste(length(unique(potentialGenePairs$gene)), " different compound heterozygous genes found."))
	
	#save name of different genes found
	# write.table(unique(outList$gene), "/media/anselm/Daten/tmp/filtern/out/30_11_1_unique_genes.txt", row.names = F, col.names = T, quote = F)
	
	#save output
	write.table (potentialGenePairs, paste(outDir, "/", outName,"_comp_het_pairs.out", sep=""), row.names=F, col.names = T, quote = F, sep = "\t")
	write.table (uniqueStart, paste(outDir, "/", outName, "_comp_het_positions.out", sep=""), row.names=F, col.names = T, quote = F, sep = "\t")
}

rm(data, gtOne, gtTwo, columnName, impact, impacts, name, outName, potentialGenePairs, uniqueStart, idFather, idMother)
