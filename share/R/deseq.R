library("DESeq")

# load arguments
ARGV <- commandArgs(TRUE)

countTable <- read.table(pipe(sprintf("uniq %s", ARGV[1])), header=TRUE, row.names=1)

if(ARGV[2]==0) {
	conds <- factor(unlist(strsplit(ARGV[3], ",", fixed=TRUE)))
	cds <- newCountDataSet( countTable, conds )
	cds <- estimateSizeFactors( cds )
	cds <- estimateDispersions(cds, fitType="local")
	res <- nbinomTest(cds, "nDP", "DP")
	res[grep(ARGV[4], res$id),]
} else if(ARGV[2]==1) {
	countTable <- round(countTable[,grep(ARGV[3], colnames(countTable))], digits=0)
	conds <- unlist(colnames(countTable))
	cds <- newCountDataSet(countTable, conds)
	cds <- estimateSizeFactors(cds)
	sizeFactors(cds)
}
