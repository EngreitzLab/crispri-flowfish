## Plot correlations in alelle counts among experimental replicates

suppressPackageStartupMessages(library(optparse))

option.list <- list(
  make_option("--samplesheet", type="character", help="Input samplesheet with all sample info"),
  make_option("--groupcol", type="character", default="ExperimentID", help="Column representing an 'experiment'; replicates of this experiment will be compared"),
  make_option("--filecol", type="character", default="ExperimentID", help="Column with paths to count files"),
  make_option("--outsummary", type="character", help="Output summary PDF with boxplots of correlation values")
)
opt <- parse_args(OptionParser(option_list=option.list))

suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

groupcol <- opt$groupcol
filecol <- opt$filecol
samplesheet <- read.delim(opt$samplesheet, check.names=F, stringsAsFactors=F)
groups <- unique(samplesheet[,groupcol])


allcor <- data.frame(Bin=c(),R=c(),stringsAsFactors=F)

for (group in groups) {
	print(group)
	currsheet <- subset(samplesheet, get(groupcol) == group)
	bins <- unique(currsheet$Bin)
	groupdata <- lapply(unique(currsheet[,filecol]), function(x) {
		if (file.exists(x)) {
			res <- read.delim(x, check.names=F, stringsAsFactors=F)
			return(res)
		}
	}) %>% setNames(unique(currsheet[,filecol]))
	groupdata[sapply(groupdata, is.null)] <- NULL

	if (length(groupdata) > 1) {
		for (bin in bins) {
			bindata <- lapply(names(groupdata), function(g) {
				if (bin %in% colnames(groupdata[[g]])) {
					res <- groupdata[[g]][c('OligoID',bin)]
					colnames(res)[2] <- basename(g)
					return(res)
				}
			})
			bindata[sapply(bindata, is.null)] <- NULL
			merged <- bindata[[1]]
			if (length(bindata) > 1) {
				for (i in 2:length(bindata))
					merged <- merged %>% merge(bindata[[i]], by='OligoID', fill=0)

				if (!is.null(opt$excludeMostFrequent)) {
					freqs <- apply(merged %>% select(-OligoID), 1, mean)
					merged <- merged[order(freqs, decreasing=T),]
				}
				cors <- cor(merged %>% select(-OligoID), use='pairwise.complete.obs')
				allcor <- rbind(allcor, data.frame(Bin=bin, R=cors[lower.tri(cors)]),stringsAsFactors=F)
			}
		}
	}
}


write.table(allcor, file=paste0(opt$outsummary,".tsv"), row.names=F, col.names=T, sep='\t', quote=F)

pdf(file=opt$outsummary, width=4, height=2)
p <- ggplot(allcor, aes(x=Bin, y=R)) + geom_boxplot() + ylim(0,1) + ylab("Replicate correlations")
print(p)
dev.off()
