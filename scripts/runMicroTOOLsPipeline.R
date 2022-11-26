#!/usr/bin/env Rscript

##
## Process the Stn. ALOHA microarrays used in
##    Open ocean and coastal strains of the N2-fixing cyanobacterium UCYN-A have
##    distinct transcriptomes
##       María del Carmen Muñoz-Marín, Jonathan D. Magasin, and Jonathan P. Zehr
##
## This script uses the MicroTOOLs microarray analysis pipeline, available at:
##     https://www.jzehrlab.com/microtools
## However, the UCYN-A microarray design (NCBI GEO platform accession GPL32341)
## is used rather than the MicroTOOLs design.  Quality control and robust
## multi-array averaging are used performed as well as the identification of
## genes with significant 24 h periodic expression.
##

options(width = 120)  # Better for printing tables

## Arguments used for our publication are indicated.
args <- commandArgs(trailingOnly=T)
snrCut <- as.numeric(args[1])           # detect genes if signal to noise ratio is >= 5
bottomPct <- as.numeric(args[2])        # bottom 10% of genes in each sample are noise
stopifnot(snrCut > 0 && bottomPct > 0)
cat("Will run the pipeline with stage 5 gene detection using snrCut=",snrCut," and bottomPct=",bottomPct,"\n\n")

suppressMessages(library(MicroTOOLs))
cfgFile <- 'Metadata/config_for_pipeline.tsv'
approachesToCarryForward <- 'ALOHA_Mari.quantilesNorm.medianpolishGene'
ooi <- c('UCYN-A1','UCYN-A2','UCYN-A3')  # Organisms of interest

##
## Simple checks on the metadata.
##
##  -- Check that SampleID's and array file names match up.
x <- read.table('Metadata/aloha_ucyna_arrays.tsv', sep="\t", header=T)
stopifnot(gsub('\\.txt\\.gz','',basename(as.character(x$FileName))) == gsub('_.*$','',as.character(x$SampleID)))
##  -- Check that descriptive part of SampleID matches the LabelName (e.g. "L3", "L6", etc.).
x <- read.table('Metadata/aloha_ucyna_samples.tsv', sep="\t", header=T)
stopifnot(gsub('^[0-9]+_','',basename(as.character(x$SampleID))) == x$LabelName)
rm(x)

PipelineStage0.Initialize(cfgFile)
PipelineStage1.ProbeQuality(cfgFile)
PipelineStage2.NormalizeArrays(cfgFile, c('quantiles'))
## Conversion happens for each of the approaches used in stage 2.
PipelineStage3.ConvertProbesToGenes(cfgFile, c('medianpolish'))
PipelineStage5.DetectGenes(cfgFile, approachesToCarryForward, snrCut=snrCut, bottomPct=bottomPct)

## Run script that evaluates ERCC levels and detected genes.  Find the Output directory in '.'
cat('\n\nComparing genes detected by pipeline to the ERCC.\n')
logFile <- 'log.howToDetect.txt'
system2('bash', input='Rscript ./scripts/howToDetect.R .', ## look for ./Output
        stdout=logFile, stderr=logFile)
file.copy(logFile, 'DetectionAnalysis')
file.remove(logFile)

dgTab <- read.csv(file=file.path('DetectionAnalysis','detectedGenes.csv'), row.names=1)
minSamp <- 4
cat('\nWill now use the ERCC analysis to filter the', nrow(dgTab), 'genes detected by the pipeline.\n',
    'Will keep genes that were detected by the pipeline and that also satisfy one or both of the following\n',
    'two tests in at least',minSamp,'samples:\n',
    '  1. The gene\'s observed level was above the observed level for the least concentrated ERCC.,\n',
    '     Mari estimated 43 transcripts before amplification for this ERCC.\n',
    '  2. The gene\'s predicted intensity is > the predicted intensity for the least concentrated ERCC.\n')
## Note that #1 compares, within each sample, the gene's level to the TGroup A ERCC "level". (The "level"
## is the minimum of the TGroup A probes because even the min's are high.  See howToDetect.R comments.)
## Test #2 uses a single predicted ERCC intensity (not per sample) because there is one linear model.
cat('The following table shows how many genes will pass the tests.  SampsAbove43mRNA counts the number\n',
    'of genes that pass test #2 in n=0,1,..16 samples. SampsAboveLowestERCC counts the genes that pass\n',
    'test #1.\n')
print(table(dgTab[,c('SampsAboveLowestERCC','SampsAbove43mRNA')]))
## If want to see how different minSamps would change the total detected genes.
## sapply(0:8, function(ms) nrow(subset(dgTab, SampsAbove43mRNA >= ms | SampsAboveLowestERCC >= ms)))
dgTab <- subset(dgTab, SampsAbove43mRNA >= minSamp | SampsAboveLowestERCC >= minSamp)
cat('After filtering there are',nrow(dgTab),'genes detected.\n')
## There are three least concentrated ERCC (each with ~43 transcripts spiked in):
##     00048:   78 spots on the array
##     00057:   79 spots
##     00075:   35 spots
## UCYN-A genes usually have 24 spots, but some have 20, 16, or fewer.  So I would
## expect the array to more easily detect ERCC if the hybridization time does not
## permit the sample to explore the array.


## Before DE stage 6, remove the genes we just filtered.  Do not bother to
## remove them from detectedGenesByEset.  (Stage 6 loads detectedGenesByEset but
## does not use them.)
cat('\n\nReplacing geneEsetList.RData so that it has the filtered genes. (Will save a "beforeErccFilter" version.\n')
geslR <- 'Output/Gene_Results/geneEsetList.RData'
load(geslR)
file.rename(geslR, 'Output/Gene_Results/geneEsetList.beforeErccFilter.RData')
stopifnot(rownames(dgTab) %in% rownames(geneEsetList[[approachesToCarryForward]]))
geneEsetList[[approachesToCarryForward]] <- geneEsetList[[approachesToCarryForward]][rownames(dgTab),]
save(geneEsetList, file=geslR)

PipelineStage6.DifferentialExpression(cfgFile, approachesToCarryForward,
                                      logFC=log2(1.2), organismsOfInterest=ooi)
cat("\nDone with with the usual pipeline stages!\n")


## Fourier analysis number one!
cat('\n\nLooking for periodic genes (Fourier analysis) using findCyclicGenes.R.',
    'The output will be in CyclicAnalysis, including a separate log file.\n')
logFile <- 'log.findCyclicGenes.txt'
system2('bash', input='Rscript ./scripts/findCyclicGenes.R .', ## look for ./Output
        stdout=logFile, stderr=logFile)
file.copy(logFile, 'CyclicAnalysis')
file.remove(logFile)


## Fourier analysis number two!
cat('\n\nWill not do an additional Fourier analysis of ALOHA but will use only the genes',
    'and timepoints (mainly) from mBio 2019.  This helps us understand why much more',
    'diel genes are detected at SIO than at ALOHA (after controlling for differences in',
    'the number of detected genes and timepoints used in the Fourier analysis).',
    'The output is in CyclicAnalysis.mbio2019ified, including a separate log file.\n')
logFile <- 'log.findCyclicGenes.mbio2019ified.txt'
system2('bash', input='Rscript ./scripts/findCyclicGenes.R . TRUE', ## look for ./Output
        stdout=logFile, stderr=logFile)
file.copy(logFile, 'CyclicAnalysis.mbio2019ified')
file.remove(logFile)

## 2020 Nov. 11: Make a table of detected genes (since the pipeline only makes
## one for DE genes).  Mari requested such a table on ~Nov. 6.  The table
## should have expression levels averaged among replicates, and feature data.
## Also tack on some results from the cyclic analysis.
stopifnot(length(approachesToCarryForward)==1)
eset <- geneEsetList[[approachesToCarryForward]]
stopifnot(colnames(exprs(eset))==rownames(pData(eset)))
MeanOfReplicates <- function(tp) {
    ## Get the replicates for time point 'tp'. Then for each gene take the mean
    ## over the replicates.
    timePoints <- as.character(pData(eset)$LabelName)
    tidx <- which(timePoints == tp)
    rowMeans(exprs(eset)[,tidx,drop=F])
}
dat <- mapply(MeanOfReplicates, unique(pData(eset)$LabelName))
stopifnot(rownames(dat) == rownames(eset))
stopifnot(colnames(dat) == unique(pData(eset)$LabelName))  # Order happens to be L3-->2L6, good.
dat <- cbind(dat, fData(eset))
## Add some cyclic analysis results.
dat$Cyclic <- 'no'
dat$Cyclic.cluster <- 'none'
cycTab <- read.csv(file.path('CyclicAnalysis','cyclicGenes.csv'), row.names=1)
stopifnot(rownames(cycTab) %in% rownames(dat))
gids <- rownames(cycTab)[rownames(cycTab) %in% rownames(dat)]
dat[gids,'Cyclic'] <- 'yes'
dat[gids,'Cyclic.cluster'] <- as.character(cycTab[gids,'Cluster'])
cat('Saving table of detected genes with cyclic status and annotation to detectedGenes.csv\n')
cat('Includes the mean (over replicates) log2 transcript levels.\n')
write.csv(dat, 'detectedGenes.csv')
rm(eset,dat,cycTab,gids)

cat("\n\nTotally DONE!  For safekeeping, you should create a directory and move all results to it:\n",
    "   -- MicroTOOLs pipeline output:  report.html and Output directory\n",
    "   -- DetectionAnalysis directory from the additional ERCC based gene detection\n",
    "   -- CyclicAnalysis from the Fourier analysis\n",
    "   -- CyclicAnalysis.mbio2019ified from the Fourier analysis restricted to mBio 2019\n",
    "      genes and timepoints\n\n")
