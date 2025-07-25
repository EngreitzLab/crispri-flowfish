#!/usr/bin/env Rscript
##################################################################
## EstimateEffectsFromBins.R
## Jesse Engreitz and Ben Doughty
##
## Based on previous script from Tejal Patwardhan, with a critical fix to the log likelihood function
## July 7, 2017
## 
## Input: count data, sort params, design file
## Black box: computes weighted average, mle mean, standard deviation for each guides
## Output: .bed and .bedgraph for each guide that passes filters
##

rm(list=ls())
suppressPackageStartupMessages(library("optparse"))
option.list <- list(
  make_option(c("-dlo", "--designDocLocation"), type="character", help="Design document location, ex: ~/Documents/Harvard/LanderResearch/flo1/raw/160705.Ess.Design.txt"),
  make_option(c("-clo", "--countsLocation"), type="character", help="Counts document location, ex: ~/Documents/Harvard/LanderResearch/flo1/raw/4.tsv"),
  make_option(c("-sp", "--sortParamsloc"), type="character", help="MUST HAVE A MEAN, BOUNDS, AND BARCODE COLUMN!!, ex: ~/Documents/Harvard/LanderResearch/flo1/raw/sortParams/gata1ff3.txt"),
  make_option(c("-om", "--outputmle"), type="character", help="File to write the MLE mean into for each guide"),
  make_option(c("-l", "--log"), type="character", help="Log results of MLE")
)
opt <- parse_args(OptionParser(option_list=option.list))

## Set up variables
designDocLocation <- opt$designDocLocation
countsLocation <- opt$countsLocation
sortParamsloc <- opt$sortParamsloc
outputmle <- opt$outputmle
log <- opt$log

MAXMEAN <- 1000000 # in FACS space, max value a guide can take
MINMEAN <- 10  # in FACS space, min value a guide can take

## load packages
# suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stats4))
suppressPackageStartupMessages(library(methods))
set.seed(100)

## FUNCTIONS

## One method of estimating the fraction of cells in the seventh bin
# takes as input the counts for allele i across the bins, the fraction of allele i in the input, and the total number of cells
# returns the number of "missing" cells (i.e. total * fraction - sum(observed))
estimateSeventhBinInput <- function(bin.counts, input.fraction, total.count) {
  return(total.count * input.fraction - sum(bin.counts))
}

## A different method of estimating the fraction of cells in the seventh bin
# takes as input a guess of starting mean+sd, bin counts and bin boundaries
# performs one step of MLE procedure, and uses this to estimate the fraction of cells that fall in 7th bin
# essentially performing one step of EM procedure
### TODO ###
# iterate EM until convergence
estimateSeventhBinEM <- function(mu.guess, sd.guess, bin.counts, bins, minmean=MINMEAN, maxmean=MAXMEAN, minvar=0, maxvar=1) {
    o <- c(bin.counts, 0) # have to artificially pad with a 0 so the likelihood function doesn't complain

    est <- mle(minuslog=ll, 
             start=list(mu=mu.guess,std=sd.guess), 
             fixed=list(observations=o, bins=bins), 
             method="L-BFGS-B",
             lower=c(log10(minmean), minvar),
             # lower=c(0, 0), since we are in log space, I think there are no constraints on the mean 
             upper=c(log10(maxmean), maxvar))#, 
             #control=list(ndeps = c(0.05,0.05)))

    MU <- est@coef[[1]]
    SI <- est@coef[[2]]
    
    ps <- vector()
    
    for (i in c(1:dim(bins)[1])) {
      ps <- c(ps, pnorm(bins[i,2], mean=MU, sd=SI, log.p=FALSE) - pnorm(bins[i,1], mean=MU, sd=SI, log.p=FALSE))
    }
    
    # use this to guess how many cells are in last bin
    seventh.bin.count <- sum(bin.counts) / sum(ps) * (1 - sum(ps))

    return(seventh.bin.count)
}

## NLL function
# takes as input a mean and standard deviation in log space, as well as a set of bin boundaries and a set of observations
# per bin (which is of length nbins+1, to account for missed cells)
# returns the negative of the log likelihood of observing those counts with the given parameters
ll <- function(mu, std, observations, bins) {
  # write("Invoked NLL",file=log,append=TRUE)
  # write(dim(bins),file=log,append=TRUE)

  #writeLines("Invoked NLL", log)
  if (std < 0) {
    return(10^10) ## make the sum really high if it picks a negative standard deviation
  }
  
  # initialize vector of probabilities of falling in bins 1..N
  pe <- vector()
  
  for (i in c(1:dim(bins)[1])) {
    pe <- c(pe, pnorm(bins[i,2], mean=mu, sd=std, log.p=FALSE) - pnorm(bins[i,1], mean=mu, sd=std, log.p=FALSE))
  }  

  # add n+1th bin = p(fall outside bin)
  pe <- c(pe, 1-sum(pe)) # add a "bin" for the remaining cells
  pe[pe <= 0] <- 10^-10 # remove any 0s from the probabilities (shouldn't happen but might)

  # assert that the lengths match
  if (length(pe) != length (observations)) {
    stop("Bin counts and bin boundaries don't have matching dimensions")
  }

  sum <- 0  # -log likelihood function
  for (i in c(1:length(observations))) {
    sum <- sum - (log(pe[i]) * observations[i]) # log(p^k)
  }

  return(sum)
}

# Wrap the maximum likelihood estimator
# getNormalMLE <- function(mu.i, sd.i, bin.counts, total.count, bins, minmean=MINMEAN, maxmean=MAXMEAN, minvar=0, maxvar=1) {
# getNormalMLE <- function(mu.i, sd.i, bin.counts, seventh.bin.count, bins, minmean=MINMEAN, maxmean=MAXMEAN, minvar=0, maxvar=1) {
getNormalMLE <- function(mu.i, sd.i, bin.counts, bins, input.present, idx, total.count, mS, minmean=MINMEAN, maxmean=MAXMEAN, minvar=0.1, maxvar=1) {
  ## mu.i = initial guess of the mean
  ## si.i = initial guess of the standard deviation
  ## bin.counts = vector of counts per bin (observed)
  ## bins = matrix of bin boundaries (nrows = nbins, ncols = 2)
  ## returns mean and standard deviation in log-space

  # cast to numeric vector
  bin.counts <- as.numeric(as.matrix(bin.counts))


  # first, need to figure out how to estimate the seventh bin counts
  seventh.bin.count <- 0
  if (input.present) {
    seventh.bin.count <- estimateSeventhBinInput(bin.counts, mS$input.fraction[idx], total.count)
    method <- "input"
  }

  seventh.bin.count <- estimateSeventhBinEM(mu.i, sd.i, bin.counts, bins)
  print(idx)
  method <- "EM"
  
  
  # add on "seventh" bin counts
  # o <- c(bin.counts, total.count - sum(bin.counts))
  o <- c(bin.counts, seventh.bin.count)

  ## Observation vector now has one entry for each bin, plus one entry for the estimated number of counts falling outside any of the bins (n+1 = "7")
  ## bins still has only n (=6) entries, for the observed bins, but the log-likelihood function `ll` adds in the n+1th bin (7th bin = 1-sum(p_i)) 

  est <- mle(minuslog=ll, 
           start=list(mu=mu.i,std=sd.i),
           # start=list(mu=MU,std=SI), 
           fixed=list(observations=o, bins=bins), 
           method="L-BFGS-B",
           lower=c(log10(minmean), minvar),
           # lower=c(0, 0), since we are in log space, I think there are no constraints on the mean and variance 
           upper=c(log10(maxmean), maxvar))#, 
           #control=list(ndeps = c(0.05,0.05)))
  MU <- est@coef[[1]]
  SI <- est@coef[[2]]

  result <- c(MU, SI, method)
  names(result) <- c("mean", "sd", "method")
  # print(result)
  return(result)
}

## LOAD DATA 
# loadReadCounts <- function(designfile, countsfile) {
loadReadCounts <- function(countsfile) {
  # designDoc <- read.delim(designfile)
  counts <- read.delim(countsfile)
  
  # create a merged sheet
  # mS <- merge(designDoc, counts)

  if (nrow(counts)==0) {
    writeLines("Unable to merge design doc and counts file", log)
    stop()
  }

  # if (! all(colnames(mS)[colnames(mS) %in% sort.params$bins$name] == sort.params$bins$name) ) {
  #   stop("Did not find columns matching bin names in the merged count table or did not find them in the same order")
  # }

  return(counts)
}

getBinNames <- function(counts) {
  col.names <- colnames(counts)
  bin.names <- col.names [! col.names %in% c('OligoID','All')]
  return(bin.names)
}

# loadSortParams_Astrios <- function(filename, mS, total.binname="Total") {
loadSortParams_Astrios <- function(filename, bin.names, total.binname="Total") {
  ## Returns a list with the following items:
  ## bins: data.frame with the following columns:
  ##   name
  ##   mean (log10)
  ##   lowerBound (log10)
  ##   upperBound (log10)
  ##   count
  ## totalCount: Total count of cells sorted

  print(bin.names)

  sort.params <- read.delim(filename)

  ## Check the sort params
  required.cols <- c("Mean","Bounds","Barcode","Count")
  if (! all(required.cols %in% colnames(sort.params)) ) {
    stop("Sort parameters file did not have the needed :", required.cols)
  }
  if (!(total.binname %in% sort.params$Barcode)) {
    stop(paste0("Barcode column in sort parameters file needs an entry called '",total.binname,"' to represent the total cell count from FACS"))
  }

  ## Extract info
  bin.indices <- which(sort.params$Barcode != total.binname)
  bin.names <- gsub("-",".",as.character(sort.params$Barcode[bin.indices]))
  bin.means <- log10(sort.params$Mean[bin.indices])
  bin.bounds <- do.call(rbind, lapply(strsplit(gsub("\\(| |\\)","", sort.params$Bounds[bin.indices]),","), as.numeric))
  total.count <- sort.params$Count[sort.params$Barcode == total.binname]
  seed <- log10(1+(sort.params[1,'Std.Dev.'] / sort.params[1,'Mean'])^2)

  # check that it has all the bins that the countsFile has
  if (sum(sort.params$Barcode %in% bin.names) != length(bin.names)) {
    stop("Sort parameters file did not have all the bins")
  }

  # if (! all(colnames(mS)[colnames(mS) %in% sort.params$bins$name] == sort.params$bins$name) ) {
  #   stop("Did not find columns matching bin names in the merged count table or did not find them in the same order")
  # }
  
  # filt.names <- vector()

  # for (i in c(1:length(bin.names))) {
  #   name <- bin.names[i]
  #   if (sum(mS[,name]) > 0) {
  #     filt.names <- c(filt.names, bin.names[i])
  #   }
  # } 

  bins <- data.frame(name=bin.names, mean=bin.means, lowerBound=log10(bin.bounds[,1]), upperBound=log10(bin.bounds[,2]), count=sort.params$Count[bin.indices], stringsAsFactors=F)
  rownames(bins) <- bin.names
  museed <- log10(bin.bounds[4,1])
  return(list(bins=bins[bin.names,], totalCount=total.count, sd.seed=seed, mu.seed=museed))
}

loadSortParams_Astrios_nototal <- function(filename, bin.names, total.binname="Total") {
  ## Returns a list with the following items:
  ## bins: data.frame with the following columns:
  ##   name
  ##   mean (log10)
  ##   lowerBound (log10)
  ##   upperBound (log10)
  ##   count
  ## totalCount: Total count of cells sorted

  print(bin.names)

  sort.params <- read.delim(filename)

  ## Check the sort params
  required.cols <- c("Mean","Bounds","Barcode","Count")
  if (! all(required.cols %in% colnames(sort.params)) ) {
    stop("Sort parameters file did not have the needed :", required.cols)
  }
  if (!(total.binname %in% sort.params$Barcode)) {
    stop(paste0("Barcode column in sort parameters file needs an entry called '",total.binname,"' to represent the total cell count from FACS"))
  }

  ## Extract info
  bin.indices <- which(sort.params$Barcode != total.binname)
  bin.names <- gsub("-",".",as.character(sort.params$Barcode[bin.indices]))
  bin.means <- log10(sort.params$Mean[bin.indices])
  bin.bounds <- do.call(rbind, lapply(strsplit(gsub("\\(| |\\)","", sort.params$Bounds[bin.indices]),","), as.numeric))
  total.count <- 0 #sort.params$Count[sort.params$Barcode == total.binname]
  seed <- log10(1+(sort.params[1,'Std.Dev.'] / sort.params[1,'Mean'])^2)

  # check that it has all the bins that the countsFile has
  if (sum(sort.params$Barcode %in% bin.names) != length(bin.names)) {
    stop("Sort parameters file did not have all the bins")
  }

  # if (! all(colnames(mS)[colnames(mS) %in% sort.params$bins$name] == sort.params$bins$name) ) {
  #   stop("Did not find columns matching bin names in the merged count table or did not find them in the same order")
  # }
  
  # filt.names <- vector()

  # for (i in c(1:length(bin.names))) {
  #   name <- bin.names[i]
  #   if (sum(mS[,name]) > 0) {
  #     filt.names <- c(filt.names, bin.names[i])
  #   }
  # } 

  bins <- data.frame(name=bin.names, mean=bin.means, lowerBound=log10(bin.bounds[,1]), upperBound=log10(bin.bounds[,2]), count=sort.params$Count[bin.indices], stringsAsFactors=F)
  rownames(bins) <- bin.names
  museed <- log10(bin.bounds[4,1])
  return(list(bins=bins[bin.names,], totalCount=total.count, sd.seed=seed, mu.seed=museed))
} 

loadSortParams_BigFoot <- function(filename, bin.names, total.binname="Total", full.file=TRUE) {
  ## Returns a list with the following items:
  ## bins: data.frame with the following columns:
  ##   name
  ##   mean (log10)
  ##   lowerBound (log10)
  ##   upperBound (log10)
  ##   count
  ## totalCount: Total count of cells sorted

  sort.params <- read.csv(filename)

  ## Check the sort params
  # check that it has all the required columns
  required.cols <- c("Mean","Min","Max","Barcode","Count")
  if (! all(required.cols %in% colnames(sort.params)) ) {
    stop("Sort parameters file did not have the needed :", required.cols)
  }
  # check that it has a total count
  if (!(total.binname %in% sort.params$Barcode)) {
    stop(paste0("Barcode column in sort parameters file needs an entry called '",total.binname,"' to represent the total cell count from FACS"))
  }
  # check that it has all the bins that the countsFile has
  if (sum(sort.params$Barcode %in% bin.names) != length(bin.names)) {
    stop("Sort parameters file did not have all the bins")
  }

  ## Extract info from sortParams for each bin
  bin.indices <- which(sort.params$Barcode %in% bin.names)
  bin.names.inorder <- sort.params$Barcode[bin.indices]
  # bin.names.inorder <- gsub("-",".",as.character(sort.params$Barcode[bin.indices]))

  if (full.file) {
    # need to treat the rows differently and cast them away from "Factors"
    bin.means <- log10(as.numeric(as.character(sort.params$Mean)[bin.indices]))
    bin.mins <- log10(as.numeric(as.character(sort.params$Min)[bin.indices]))
    bin.maxs <- log10(as.numeric(as.character(sort.params$Max)[bin.indices]))
    seed <- log10(1+(as.numeric(as.character(sort.params$StdDev[sort.params$Barcode == total.binname])) / as.numeric(as.character(sort.params$Mean[sort.params$Barcode == total.binname])))^2)
  } else {
    bin.means <- log10(sort.params$Mean[bin.indices])
    bin.mins <- log10(sort.params$Min[bin.indices])
    bin.maxs <- log10(sort.params$Max[bin.indices])
    seed <- log10(1+(sort.params$StdDev[sort.params$Barcode == total.binname] / sort.params$Mean[sort.params$Barcode == total.binname])^2)
  }
  
  bin.counts <- sort.params$Count[bin.indices]

  # compute some extra numbers which will help our estimation
  total.count <- sort.params$Count[sort.params$Barcode == total.binname]

  bins <- data.frame(name=bin.names.inorder, mean=bin.means, lowerBound=bin.mins, upperBound=bin.maxs, count=bin.counts, stringsAsFactors=F)
  rownames(bins) <- bin.names.inorder
  return(list(bins=bins, totalCount=total.count, sd.seed=seed))
}

loadSortParams_BigFoot_noTotal <- function(filename, bin.names, full.file=FALSE) {
  ## Returns a list with the following items:
  ## bins: data.frame with the following columns:
  ##   name
  ##   mean (log10)
  ##   lowerBound (log10)
  ##   upperBound (log10)
  ##   count
  ## totalCount: Total count of cells sorted

  sort.params <- read.delim(filename)

  # print(sort.params)

  ## Check the sort params
  # check that it has all the required columns
  required.cols <- c("Mean","Min","Max","Bin","Count")
  if (! all(required.cols %in% colnames(sort.params)) ) {
    stop("Sort parameters file did not have the needed :", required.cols)
  }
  # check that it has a total count
  # if (!(total.binname %in% sort.params$Barcode)) {
  #   stop(paste0("Barcode column in sort parameters file needs an entry called '",total.binname,"' to represent the total cell count from FACS"))
  # }
  # check that it has all the bins that the countsFile has
  # if (sum(sort.params$Bin %in% bin.names) != length(bin.names)) {
    # stop("Sort parameters file did not have all the bins")
  # }

  ## Extract info from sortParams for each bin
  bin.indices <- which(sort.params$Bin %in% bin.names)
  bin.names.inorder <- sort.params$Bin[bin.indices]
  # bin.names.inorder <- gsub("-",".",as.character(sort.params$Barcode[bin.indices]))

  if (full.file) {
    # need to treat the rows differently and cast them away from "Factors"
    bin.means <- log10(as.numeric(as.character(sort.params$Mean)[bin.indices]))
    bin.mins <- log10(as.numeric(as.character(sort.params$Min)[bin.indices]))
    bin.maxs <- log10(as.numeric(as.character(sort.params$Max)[bin.indices]))
    seed <- 1 # log10(1+(as.numeric(as.character(sort.params$StdDev[sort.params$Barcode == total.binname])) / as.numeric(as.character(sort.params$Mean[sort.params$Barcode == total.binname])))^2)
  } else {
    bin.means <- log10(sort.params$Mean[bin.indices])
    bin.mins <- log10(sort.params$Min[bin.indices])
    bin.maxs <- log10(sort.params$Max[bin.indices])
    seed <- 0.5 ## HOW TO CHOOSE SEED??? # log10(1+(sort.params$StdDev[sort.params$Barcode == total.binname] / sort.params$Mean[sort.params$Barcode == total.binname])^2)
  }
  
  bin.counts <- sort.params$Count[bin.names.inorder]

  # compute some extra numbers which will help our estimation
  total.count <- 0 # sort.params$Count[sort.params$Barcode == total.binname]

  mu.seed <- bin.mins[3]

  bins <- data.frame(name=bin.names.inorder, mean=bin.means, lowerBound=bin.mins, upperBound=bin.maxs, count=bin.counts, stringsAsFactors=F)
  rownames(bins) <- bin.names.inorder
  return(list(bins=bins, totalCount=total.count, sd.seed=seed, mu.seed=mu.seed))
}

addMetaData <- function(mS, bin.names) {
  mS$sumreads <- rowSums(mS[,bin.names])
  mS$numobsbins <- rowSums(mS[,bin.names]>0)
  return(mS)
}



rescaleReadCounts <- function(mS, sort.params, bin.names, input="InputCount") {
  binNames <- bin.names
  # write.table(mS[,binNames]/colSums(mS[,binNames]), file='/seq/lincRNA/Ben/base_editing/190415_BFF_spike_test/IL2RA-2/rs61839660-TATCTATTTTGGTCCCAAAC.NoSpike.BEFORE.txt', sep="\t", quote=F, row.names=F, col.names=T)
  mS[,binNames] <- sapply(binNames, function(binName) mS[,binName] / sum(mS[,binName]) * sort.params$bins[binName,"count"])
  # write.table(mS[,binNames], file='/seq/lincRNA/Ben/base_editing/190415_BFF_spike_test/IL2RA-2/rs61839660-TATCTATTTTGGTCCCAAAC.NoSpike.AFTER.txt', sep="\t", quote=F, row.names=F, col.names=T)
  mS[is.na(mS)] <- 0
  return(mS)
}


addWeightedAverage <- function(mS, sort.params, bin.names) {
  mS$sumcells <- rowSums(mS[,bin.names])
  sort.params$bins$mean <- 0.5 * (sort.params$bins$lowerBound + sort.params$bins$upperBound)
  sort.params$bins['A','mean'] <- sort.params$bins['A','upperBound']
  sort.params$bins['F','mean'] <- sort.params$bins['F','lowerBound']

  mS$WeightedAvg <- ifelse(mS$sumcells == 0, 0, rowSums(t(sort.params$bins[bin.names,'mean'] * t(as.matrix(mS[,bin.names])))) / mS$sumcells)
  return(mS)
}

## Set up read count table and sort params
# mS <- loadReadCounts(designDocLocation, countsLocation)
counts <- loadReadCounts(countsLocation)
bin.names <- getBinNames(counts)
# sort.params <- loadSortParams_Astrios_nototal(sortParamsloc, bin.names)
sort.params <- loadSortParams_BigFoot_noTotal(sortParamsloc, bin.names)
counts <- addMetaData(counts, bin.names)
counts <- rescaleReadCounts(counts, sort.params, bin.names)
counts <- addWeightedAverage(counts, sort.params, bin.names)

# write.table(mS, file='/seq/lincRNA/Ben/base_editing/190415_BFF_spike_test/IL2RA-2/rs61839660-TATCTATTTTGGTCCCAAAC.NoSpike.AFTER.txt', sep="\t", quote=F, row.names=F, col.names=T)

# run MLE
runMLE <- function(mS, sort.params) {
  ## seed mean and standard deviation
  mu.seed <- sort.params$mu.seed #.5 
  print(mu.seed)
  sd.seed <- sort.params$sd.seed #.5 
  # print(sd.seed)
  bin.bounds <- as.matrix(sort.params$bins[,c("lowerBound","upperBound")])
  total.count <- sort.params$totalCount

  input.bin.name <- 'All'
  input.present <- FALSE

  if (input.bin.name %in% colnames(mS)) {
    mS$input.fraction <- mS[,input.bin.name] / sum(mS[,input.bin.name])
    input.present <- TRUE
  }

  ############# MAKE SURE BIN.BOUNDS AND mS[i,bin.names] ARE SORTED IN THE SAME ORDER!!! #############
  mleOuts <- data.frame(do.call(rbind, lapply(1:nrow(mS), function(i) 
    #tryCatch({ getNormalMLE(mS$WeightedAvg[i], sd.seed, mS[i,sort.params$bins$name], mS[i,inputCountCol], bin.bounds) }, error = function(err) { writeLines(paste0("MLE errored out at given initialization for guide: ", i, ", using weighted average instead"), log); result <- c(mS$WeightedAvg[i], sd.seed); names(result) <- c("mean", "sd"); return(result) }))))
    # tryCatch({ getNormalMLE(mS$WeightedAvg[i], sd.seed, mS[i,sort.params$bins$name], mS[i,inputCountCol], bin.bounds) }, error = function(err) { write(paste0("MLE errored out at given initialization for guide: ", i, ", using weighted average instead: ", err),file=log,append=TRUE); result <- c(mS$WeightedAvg[i], sd.seed); names(result) <- c("mean", "sd"); return(result) }))))
    # estimate the number of cells falling in n+1th bin
    # estimated.seventh.bin <- estimateSeventhBinInput(mS[i,bin.names], mS[i,input.bin.name] / sum(mS[,input.bin.name]), total.count)    
    tryCatch({ getNormalMLE(mu.seed, sd.seed, mS[i,bin.names], bin.bounds, input.present, i, total.count, mS) }, error = function(err) { write(paste0("MLE errored out at given initialization for guide: ", i, ", using weighted average instead: ", err),file=log,append=TRUE); result <- c(mS$WeightedAvg[i], sd.seed, "weighted_avg"); names(result) <- c("mean", "sd", "method"); return(result) }))))
  
  mS$logMean <- mleOuts$mean
  mS$logSD <- mleOuts$sd
  mS$method <- mleOuts$method
  
  return(mS)
}


# write(pe[7],file=log,append=TRUE)
counts <- runMLE(counts, sort.params)

write.table(counts, file=outputmle, sep="\t", quote=F, row.names=F, col.names=T)
write("Test",file=log,append=TRUE)
#writeLines("Test", log)
