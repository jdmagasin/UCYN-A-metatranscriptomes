#!/usr/bin/env Rscript

##
## Detect cyclically expressed genes using the Cycle package, similar to the
## vignette for the cycle package.
##

args <- commandArgs(trailingOnly=T)
pipeDir <- args[1]  # can pass in '.'
stopifnot(dir.exists(pipeDir))

## Make the background models reproducible, Since backgroundData() seems to just
## call rnorm(), set.seed() should do the job.
set.seed(8675309)

mbio2019ify <- FALSE  # Set true if you want the analysis done with respect to the
                      # genes detected in the mBio 2019 paper (and at ALOHA of course)
                      # and also to drop two of the timepoints so that there are 8 total
                      # as in mBio2019, rather than 10 (at ALOHA).
if (length(args) > 1 && (args[2] %in% c('T','TRUE'))) {
    mbio2019ify <- TRUE
}
## More info on mbio2019ify:
##   I'm trying to figure out why despite much more detected genes at ALOHA I find a
##   much lower % of diel genes. Set this flag to make the Fourier analysis be with
##   respect to the ~same genes detected as in mBio 2019 and the ~same timepoints.
##   This should control for effects that more timpepoints or detected genes have on
##   the background distribution of Fourier scores.
## Note that if this flag is T, then the script will look in the working directory
## for:  OtherDielExprStudies/MunozMartin_mBio2019/detectedGenes.txt

if (!mbio2019ify) {
    outDir <- file.path(pipeDir,'CyclicAnalysis')
} else {
    outDir <- file.path(pipeDir,'CyclicAnalysis.mbio2019ified')
}
dir.create(outDir)

## The pipeline records which genes were detected in each sample in detGenes.
## Should I use detGenes to set to NA the intensities of genes that were not
## detected in each sample?  If not, then low (normalized) values will be
## used. If so, then the code below will impute values (as in the vignette).
## For PipeRun.snrCut5.bottomPct10, I tried both with and without and it does
## not make much difference.  191 or 192 cyclic genes are identified, and 7 of
## them are specific to the former or latter approach.
setUndetectedGenesToNA <- TRUE

suppressMessages(library(MicroTOOLs))
load(file.path(pipeDir,'Output/Gene_Results/geneEsetList.RData'))
eset <- geneEsetList[[1]]
rm(geneEsetList)
load(file.path(pipeDir,'Output/Gene_Results/detectedGenesByEset.RData'))
detGenes <- detectedGenesByEset[[1]]
rm(detectedGenesByEset)


## Before July 7, 2020, I used knn.
Impute.knn <- function(eset)
{
    cat('Imputing missing values using knn.\n')
    ##eset <- fill.NA(eset, mode="wknn")  # Leaves some NA's (?)
    eset <- fill.NA(eset, mode="knn")
    ## Vignette suggests wknn or knn. There is also mean, but would taking the
    ## mean across light and dark samples flatten the gene's expression profile?
    eset
}

Impute <- Impute.knn


if (setUndetectedGenesToNA) {
    ## Set to NA genes that were not detected in the sample. See note above.
    ## NOTE: detGenes is as output by the pipeline and so it includes genes that
    ## did not pass the ERCC-based filtering.  That's okay because genes in the
    ## eset are only those that passed ERCC-based filtering.
    cat('Genes not detected by the pipeline in a specific sample will now be set to NA.\n')
    stopifnot(names(detGenes) == colnames(eset))
    for (sn in colnames(eset)) {
        nidx <- which(!rownames(eset) %in% names(detGenes[[sn]]))
        cat('Setting',length(nidx),'undetected genes to NA in sample',sn,'\n')
        exprs(eset)[nidx,sn] <-NA
    }
    ## Remove genes with > 25% undetected.  Then impute values for retained genes.
    cat('Dropping genes that have NA in > 25% of the samples.\n')
    eset <- filter.NA(eset, thres=0.25)  # (5 genes dropped for PipeRun.snrCut5.bottomPct10)
    eset <- Impute(eset)
}
stopifnot(!is.na(exprs(eset)))

## Define the hours since time 0 (65308_L3).  I did this manually by looking at pData
sampHours <- c(0,3,3,6,6,9,9,12,15,18,18,21,21,24,27,27)
stopifnot(length(sampHours) == ncol(eset))

## Now combine replicate times. If you do not, then they will be treated as
## independent (see vignette).  (Could transpose exprs(eset) and then
## aggregate(), but I like the progress messages.)
for (sh in unique(sampHours)) {
    sidx <- which(sampHours == sh)
    if (length(sidx) == 1) next
    ## The first sample at the time gets the means of all samps at that time
    cat('Taking gene means samples from hour=',sh,'. Samples are',paste0(colnames(eset)[sidx],collapse=','),'\n')
    exprs(eset)[,sidx[1]] <- rowMeans(exprs(eset)[,sidx])
}
## Drop the non-first-at-hour samples
eset <- eset[,match(unique(sampHours), sampHours)]
colnames(eset) <- gsub('[0-9]+_', '', colnames(eset))  # Or just replace with the LabelName:)
## Fixup sampHours so that it is in sync with eset's colnames
sampHours <- unique(sampHours)
stopifnot(ncol(eset) == length(sampHours))
## Critical: sampHours must be ordered as the samples (cols) in the eset)

## Here is where I drop 2 timepoints from the eset, going from 10 to 8. The SIO
## (mBio 2019) data have just 8 timepoints. Maybe that is why they have lower
## Fourier scores than ALOHA (means of 4.3 vs. 6.0 for all diel genes).  If
## fewer timepoints --> lower scores, then that should also be true for the
## background scores.  *If* somehow the background scores decrease FASTER with
## fewer timepoints, that might explain why proportionally (to detected) more
## diel genes were found in mBio 2019 than at ALOHA.
if (mbio2019ify) {
    ## SIO:          L6   L9          D3    D6       2D12   2L3   2L9  2L12
    ## ALOHA:   L3   L6   L9   L12    D3    D6   D9   D12   2L3   2L6
    ## drop      x                               x
    eset <- eset[,setdiff(colnames(eset), c('L3','D9'))]
    cat('Threw out L3 and D9 so that ALOHA would have 8 timepoints, similar to mBio 2019.\n')
    sampHours <- sampHours[-c(1,7)]
}

## Here is where I optionally limit the eset to just the UCYN-A genes detected
## in mBio 2019.
if (mbio2019ify) {
    mbioLocusTags <- readLines('OtherDielExprStudies/MunozMartin_mBio2019/detectedGenes.txt')
    cat('Eset has',nrow(eset),'genes before restricting to the',length(mbioLocusTags),
        'genes detected in mBio 2019 paper.\n')
    eset <- eset[fData(eset)$ORF %in% mbioLocusTags,]
    cat('After restricting, eset now has',nrow(eset),'genes.\n')
}

## Standardize expression levels as in the vignette.
cat('Standarizing each gene\'s intensities across the samples.\n')
## Actually, fdrfourier() does this anyway -- peek at the source.  No harm though.
eset.s <- standardise(eset)  # Note that I'm standardizing the log2 intensities.

## From vignette, when explaining that levels are autocorrelated so Gaussian and
## randomized null models lead to overestimates of periodic genes.
cat('Later we will use the AR1 background model to get the same autocorrelation anyway.',
    'Let\'s check whether autocorrelation is high enough to worry about...\n')
auto.corr <- 0
for (i in 2:dim(exprs(eset.s))[[2]]){
    auto.corr[i] <- cor(exprs(eset.s)[,i-1],exprs(eset.s)[,i])
}
auto.corr  # Ours are low (max of 0.32 for PipeRun.snrCut5.bottomPct10).

## Calculate FDR based by comparison to the specified model.
## Cycle periodHrs and sample timpepoints are both in hours
## Critical: sampHours must be ordered as the samples (cols) in the eset)
DoFourierAnalysis <- function(eset, timepoints, periodHrs=24, fmod='ar1', NN=1000)
{
    fdrFourRes <- fdrfourier(eset, T=periodHrs, times=timepoints,
                             background.model=fmod, N=NN, progress=TRUE)
    lowFdrGenes <- which(fdrFourRes$fdr < 0.25)
    cat('Fourier analysis with',fmod,'background identified',length(lowFdrGenes),
        'periodic genes (FDR<0.25)\n')
    fdrFourRes
}

cat('Now the real work.  For each gene calculate its Fourier score and FDR using the',
    'AR1 background model.  Make 1000 background data sets.  Get some coffee!\n')
fdrFourRes.ar1 <- DoFourierAnalysis(eset.s, timepoints=sampHours)
lowFdrGenes <- which(fdrFourRes.ar1$fdr < 0.25)
length(lowFdrGenes)
genes.cycle <- fdrFourRes.ar1$fdr[lowFdrGenes]  # genes with diel pattern
if (TRUE) {
    ## I'm using AR1, but also handy to do the other backgrounds for comparison.
    fdrFourRes.rr <- DoFourierAnalysis(eset.s, timepoints=sampHours, fmod='rr')
    fdrFourRes.gauss<- DoFourierAnalysis(eset.s, timepoints=sampHours, fmod='gauss')
    save(eset,fdrFourRes.ar1, fdrFourRes.rr, fdrFourRes.gauss,
         file=file.path(outDir,paste0('find_Cyclic_Genes.Rdata')))
} else {
    save(eset,fdrFourRes.ar1, file=file.path(outDir,paste0('find_Cyclic_Genes.Rdata')))
}
cat('That took a long time. Saved the Fourier analysis objects to',outDir,'\n')


## Check F vs FDR
png(file.path(outDir,'fourierScoreVsFDR.png'))
plot(fdrFourRes.ar1$fdr, fdrFourRes.ar1$F, xlab="FDR", ylab="Fourier score", main="FDR (ar1 background)")
dev.off()
## FDR distribution
png(file.path(outDir,'fdrDensity.png'))
plot(density(fdrFourRes.ar1$fdr), main="FDR (ar1 background)")
dev.off()

##
## Clustering of just the cyclic genes
##  -- If you run pvclust() multiple times, you will get different numbers of significant (95%)
##     clusters.  (centroid, correlation).  One big one, and 11 to 15 smaller ones.
##     If these smaller clusters look similar, might be better to plot their cycles together
##     and note that merged the clusters.
##  -- Would it be interesting to cluster cyclic with non-cyclic (non-signif F score) genes?
##
library(pvclust)
gids.cycle <- names(genes.cycle)
eset.s.cycle <- eset.s[gids.cycle,]
##hc <- hclust(dist(1 - cor(t(exprs(eset.s.cycle)), method='pearson'), method = "euclidean"))

clustMeth <- 'centroid'  # ward.D2: most signif clusters have very few genes.
                         # centroid: results in some large clusters (>70 genes) and smaller ones
                         # average:  some of the signif clusters have lots of genes.
nboot <- 1000  # 1K yields very tiny standard errors (see below)
cat('Using pvclust to cluster the cyclic genes (',clustMeth,')',
    'and then multiscale bootstrap (distance based on correlation,',nboot,'bootstraps)\n')
pvc <- pvclust(t(exprs(eset.s.cycle)),
               method.hclust=clustMeth,
               method.dist='correlation',  # also 'abscor'
               nboot=nboot,
               parallel=T)
png(file.path(outDir,'cyclicGeneDendrogram.png'), width=18, height=7, units='in', res=72)
plot(pvc, main='Cyclic genes dendrogram (p-values as%)', cex=0.5)
alpha <- 0.95
maxOnly <- TRUE  # Do not report clusters that are subsets of larger clusters
pvrect(pvc, alpha=alpha, pv="au", type="geq", max.only=maxOnly)
dev.off()

## Get the significant clusters (same as pvrect).  Rename them by their edges.
cat('Picking significant clusters (approx. unbiased p-value % >=',round(100*alpha),')\n')
clusts <- pvpick(pvc, alpha=alpha, pv='au', type='geq', max.only=maxOnly)
stopifnot(length(clusts$clusters) == length(clusts$edges))
names(clusts$clusters) <- paste0('edge.',clusts$edges)

## Evaluate the accuracy of the approximately unbiased (au) p-values.  (See paper,
## and ?print.pvclust and ?seplot and ?msplot).
cat('Here are approx unbiased p-values (au) and associated standard errors (au.se)',
    'for the significant edges (alpha >',alpha,')\n')
print(pvc, which=clusts$edges)

png(file.path(outDir,'pvclustCurveFitsForSignifEdges.png'), width=16, height=10, units='in', res=72)
msplot(pvc, edges=clusts$edges)  # would-be-title='Curve fitting for significant clusters'
dev.off()

## Summarize all the significant clusters.
cat('\nThese are the number of genes in each cluster:\n')
sort(sapply(clusts$clusters, length), decreasing=T)
cat('\n\nThese are the kinds of genes in each cluster:\n')
sapply(clusts$clusters, function(gids) { sort(unique(as.character(fData(eset.s)[gids,'GENE']))) })

cat('Keeping clusters with >= 5 genes.\n')
clusts <- clusts$clusters[names(which(sapply(clusts$clusters, length) >= 5))]


## Plot the clusters!
library(ggplot2)
library(reshape2)
x <- c('L3', 'L6', 'L9', 'L12', 'D3', 'D6', 'D9', 'D12' , '2L3', '2L6')
st2hr <- c(0, cumsum(rep(3,length(x)-1)))  # hours since L3
names(st2hr) <- x
stopifnot(names(st2hr) == colnames(eset.s))
dusk <- st2hr['D3']-3   # D3 is 3 hrs into the night, so dusk is 3 hrs before the D3 hour
dawn <- st2hr['2L3']-3  # 2L3 is 3 hrs into day
gList <- list()
for (cn in names(clusts)) {
    ## Plot standardized levels.
    dat <- melt(exprs(eset.s.cycle)[clusts[[cn]],])
    colnames(dat) <- c('GID','SampTime','Z')
    dat$Hour <- st2hr[dat$SampTime]
    ## GENE: as.char is critical!
    x <- fData(eset.s.cycle)[as.character(dat$GID), 'GENE']; x[x==''] <- 'unknown'
    x <- substr(x, 1, 10)  # <= 10 letters of gene name so legend width reasonable
    ## Genes with >1 appearance get [num] in their name. Divide to ~unmelt.
    gt <- sort(table(x)/length(levels(dat$SampTime)), decreasing=T)
    stopifnot(sum(gt) == length(clusts[[cn]]))
    x <- paste0(x, ' [',gt[x],']')
    x <- gsub(' \\[1\\]$', '', x)
    dat[['Gene']] <- factor(x, levels= names(sort(table(x),decreasing=T)))
    ## ORGANISM: as.char is critical!
    x <- fData(eset.s.cycle)[as.character(dat$GID), 'ORGANISM']; x[x=='UCYNA2'] <- 'UCYNA3'
    dat[['Organism']] <- factor(x)
    ## Each GID is a group so that we get 1 line per GID, but color by property
    g <- ggplot(dat, aes(x=Hour, y=Z, group=GID)) +
        annotate('rect', xmin=dusk, xmax=dawn, ymin=-Inf, ymax=Inf, fill='navyblue', alpha=0.20) +
        geom_line(aes(color=Gene)) + theme_bw() +
        scale_x_continuous(breaks=st2hr, labels=(st2hr + 9) %% 24) +
        labs(title=paste('Gene cluster',cn), y='Intensity (Z score)', x='Hour of day') +
        facet_wrap(~Organism, nrow=2)
    gList[[cn]] <- g
}

for (cn in names(gList)) {
    png(file.path(outDir,paste0('cluster_',cn,'.png')),
        width=8, height=4*length(levels(gList[[cn]]$data$Organism)), res=72)
    plot(gList[[cn]])
    dev.off()
}
                
## Make a table of cyclic genes.  If a gene is in one of the four main clusters,
## note that.
g2c <- rep('other',length(gids.cycle))
names(g2c) <- gids.cycle
for (x in names(clusts)) { g2c[clusts[[x]]] <- x }
dat <- data.frame(exprs(eset.s)[gids.cycle,],
                  fData(eset.s)[gids.cycle,],
                  Cluster = g2c[gids.cycle],
                  FScore  = fdrFourRes.ar1$F[gids.cycle],
                  FDR     = fdrFourRes.ar1$fdr[gids.cycle])
write.csv(dat, file=file.path(outDir,'cyclicGenes.csv'))


## Quick way to show that genes (names) are specific to each cluster (signif
## clusters with >5 genes). 'unknown' is the 1 gene name shared by most clusters.
library(venn)  # can handle more sets than gplots::venn
png(file.path(outDir,'venn_geneNameCounts_signifClustersWithAtLeast5genes.png'))
x <- lapply(clusts, function (gids) unique(fData(eset)[gids,'GENE']))
x <- x[order(sapply(x,length), decreasing=T)]  # b/c at most 1st 7 clusters drawn
venn(x)
dev.off()

quit(save='no')
