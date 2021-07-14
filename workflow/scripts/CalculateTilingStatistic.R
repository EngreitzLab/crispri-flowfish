##################################################################
## Jesse Engreitz
## April 28, 2016
##
## Example command:
## use .r-3.1.1-bioconductor-3.0; R_LIBS_SITE=; reuse Python-2.7
## Rscript CalculateTilingStatistic.R --input /seq/lincRNA/cfulco/MYCTiling/160303.EssentialityTiling/160328.HiSeq.CRISPRiScreensMoreMaterial/160329.ScreenSummaries/160501.CiTile.database.GATA1test.bed --output /seq/lincRNA/RAP/Paper/MYCTiling/Figures/WindowSize/test.bed --scoreColumn log2FC --T0columns RAW_CiTile.1.T0,RAW_CiTile.8.T0


suppressPackageStartupMessages(library("optparse"))

option.list <- list(
  make_option(c("-i", "--input"), type="character", help="Input database file containing at least the following columns: chr, start, end, GuideSequence, an arbitrarily named score column, and a 'target' column that includes negative_control labeled guides"),
  make_option(c("-o", "--output"), type="character", help="Output tab-delimited text file"),
  make_option(c("--T0columns"), type="character", help="Comma-separated list of columns to apply the minT0count filter"),
  make_option(c("-c", "--scoreColumn"), type="character", default="log2FC", help="Column name calculating scores and statistics"),
  make_option(c("-w", "--grnaWindowSize"), type="numeric", default=20, help="Number of guides to tile in windows"),
  make_option(c("-m", "--maxSpan"), type="numeric", default=500, help="Max span of the guides in a given window"),
  make_option(c("-s", "--minOffTargetScore"), type="numeric", default=50, help="Minimum off-target score"),
  make_option(c("-t", "--minT0count"), type="numeric", default=50, help="Minimum T0 counts"),
  make_option(c("--minGuideLength"), type="numeric", default=20, help="Minimum length of guideRNAs to use"),
  make_option(c("--maxGuideLength"), type="numeric", default=21, help="Maximum length of guideRNAs to use"),
  make_option(c("--maxGsInGuideSequence"), type="numeric", default=10, help="Maximum number of G's allows in guide sequence"),
  make_option("--codeDir", type="character", default="/seq/lincRNA/RAP/Promoters/LanderLab-EP-Prediction/", help="Absolute path to code repository")
)
opt <- parse_args(OptionParser(option_list=option.list))

suppressPackageStartupMessages(source(paste0(opt$codeDir, "/src/libs/JuicerUtilities.R")))
suppressPackageStartupMessages(source(paste0(opt$codeDir, "/src/CRISPRScreen/JuicerCRISPRi.R")))



## PARAMETERS
t0columns <- strsplit(gsub("-",".",opt$T0columns), ",")[[1]]
window.size <- opt$grnaWindowSize
max.span <- opt$maxSpan
data.col <- gsub("-",".",opt$scoreColumn)

x <- read.delim(gzfile(opt$input))
stopifnot(all(c("chr","start","end",data.col,"target","GuideSequence") %in% colnames(x)))

x$nG <- sapply(as.character(as.matrix(x$GuideSequence)), function(x) sum(strsplit(x, "")[[1]] == "G"))
x$GuideLength <- sapply(as.character(as.matrix(x$GuideSequence)), nchar)

to.use <- subset(x, OffTargetScore >= opt$minOffTargetScore & 
                    nG <= opt$maxGsInGuideSequence & 
                    GuideLength >= opt$minGuideLength & 
                    GuideLength <= opt$maxGuideLength)
for (col.name in t0columns) {
  to.use <- subset(to.use, get(col.name) >= opt$minT0count)
}
cat(nrow(to.use),"of",nrow(x),"total guideRNAs passed the filters\n")
stopifnot(nrow(to.use) > 0)

control.vals <- subset(to.use, target == "negative_control")[,data.col]
cat("Using",length(control.vals),"as negative control guides\n")
stopifnot(length(control.vals) > 0)
# write.table(to.use, file=paste0(opt$input,".filtered.guides.txt"), sep='\t', quote=F, row.names=F, col.names=T)


result <- calculateSlidingWindowStatistic(subset(to.use, target != "negative_control"), control.vals, window.size, max.span, data.col=data.col)

write.table(result, file=opt$output, sep='\t', quote=F, row.names=F, col.names=T)

to.write.mean <- with(result, data.frame(chr=chr, start=start, end=end, score=-mean))
write.table(to.write.mean, file=paste0(opt$output,".mean.bedgraph"), sep='\t', quote=F, row.names=F, col.names=F)

# to.write.u <- with(result, data.frame(chr=chr,start=start,end=end,score=log10(fdr.utest)*sign(mean)))
# write.table(to.write.u, file=paste0(opt$output,".UFDR.bedgraph"), sep='\t', quote=F, row.names=F, col.names=F)

# to.write.gRNA <- with(subset(to.use,!is.na(start)), data.frame(chr=chr, start=start, end=end, score=-get(data.col)))
# write.table(to.write.gRNA, file=paste0(opt$output,".IndividualGuides.bedgraph", sep='\t', quote=F, row.names=F, col.names=F)


