##
## Interactive script for main statistical tests in
##    Open ocean and coastal strains of the N2-fixing cyanobacterium UCYN-A have
##    distinct transcriptomes
##       María del Carmen Muñoz-Marín, Jonathan D. Magasin, and Jonathan P. Zehr
##

## Manifest:
##   Objects
##     aloha         - Has three slots:
##                      eset     -- ExpressionSet of Stn. ALOHA log2 transcript
##                                  levels.
##                      cycGenes -- Data frame of diel genes.  FDR and Fourier
##                                  score from cycle::fdrfourier(), Cluster from
##                                  pvclust().
##     sio           - Similar to 'aloha' but for the Scripps Pier data from
##                     Muñoz-Marín et al. (2019, mBio).  Cluster from mBio 2019.
##
##     platTab       - Data frame describing the microarray.
##
##     gids.A{1,2,3} - Gene identifiers (GIDs) for UCYN-A1, A2, and A3.
##
##     taraTab         - Data for South Atlantic TARA Oceans Station 78, adapted
##                       from Cornejo-Castillo et al. (2016).
##
##   Functions:
##     GetExprs               - Get a matrix of transcript levels for 'aloha' or
##                              'sio' for just the specified GIDs.
##     GetPeakAndTroughTimes  - Returns a data frame describing the time points
##                              at which maximum and minimum transcript levels
##                              were seen.
##

library(Biobase)
load('Munoz-Marin_Magasin_Zehr_PLOS_ONE.Rdata')

set.seed(19200921)

stopifnot( nrow(aloha$eset) == 1939 ) # detected genes at Stn. ALOHA
stopifnot( nrow(sio$eset)   ==  762 ) # detected genes at Scripps Pier
stopifnot( nrow(platTab)    == 2753 ) # Num genes on the Stn. ALOHA microarray
stopifnot( length(gids.A1)  == 1195 ) #  A1 genes on the microarray
stopifnot( length(gids.A2)  == 1244 ) #  A2 genes on the microarray
stopifnot( length(gids.A3)  ==  314 ) #  A3 genes on the microarray


## Genes detected at both locations are used in several tests.
crossSiteDetectedGids <- intersect(rownames(aloha$eset),rownames(sio$eset))
stopifnot(length(crossSiteDetectedGids) == 760)


################################################################################
##
## Results section 1: UCYN-A2 had higher transcripts per cell than other
## sublineages at Stn. ALOHA
##

##
## Paired t-test of 272 genes detected for both UCYN-A2 and A1. Compares
## orthologs based on gene symbols Restrict to the >0.8um fraction to pair genes
## from the same sample.  (There was no A2 in the 0.2-3um sample and no A1 in
## the 5-20um sample).
##
agg <- aggregate(cbind(metagenome.normalized,gyrB.normalized) ~ GENE + ORGANISM,
                 subset(taraTab, Filter=='>0.8um'), median)
agg <- subset(agg, !GENE%in%c('','0'))
agg <- agg[-grep('^E[0-9]+\\.', agg$GENE),]  # Drop vague enzymes
nrow(agg)  # 893
## Keep genes that we have for both A1 and A2
x <- intersect(subset(agg, ORGANISM=='UCYNA1')$GENE,
               subset(agg, ORGANISM=='UCYNA2')$GENE)
agg <- subset(agg, GENE %in% x)
nrow(agg) # 812 (length(x)==406 genes)
df <- data.frame(UCYNA1 = subset(agg, ORGANISM=='UCYNA1')[,"metagenome.normalized"],
                 UCYNA2 = subset(agg, ORGANISM=='UCYNA2')[,"metagenome.normalized"],
                 GENE1  = subset(agg, ORGANISM=='UCYNA1')[,"GENE"],
                 GENE2  = subset(agg, ORGANISM=='UCYNA1')[,"GENE"])
stopifnot(as.character(df$GENE1) == as.character(df$GENE2))
nrow(df) # 406
df <- subset(df, UCYNA1>0 & UCYNA2>0)  # As above, only examine exrpessed genes.
nrow(df) # 272 if just the >0.8um fraction
a2 <- log2(df$UCYNA2)
a1 <- log2(df$UCYNA1)
qqnorm(a2); qqline(a2)  # For t-test, check normality. Seems so for most a2 and a1
qqnorm(a1); qqline(a1)  # observations.
t.test(a2, a1, paired=T, alternative='greater')        # PAIRED, ONE-SIDED
wilcox.test(a2, a1, paired=T, alternative='greater')   # Alternative if disagree that a2 and a1 are
                                                       # approximately normal.
## Either way, p ~ 0 and the alternative hypothesis seems more likely.
## Gene-for-gene, UCYN-A2 transcript levels are higher than UCYN-A1 levels.
rm(agg,df,x,a1,a2)


##
## Unpaired t-test of all genes detected for both UCYN-A2 and A1
##
## Get normalized counts for detected genes.
a2 <- subset(taraTab, ORGANISM=='UCYNA2' & metagenome.normalized > 0)
a1 <- subset(taraTab, ORGANISM=='UCYNA1' & metagenome.normalized > 0)
## Count detected gene targets. This is less than the number of observations.
length(unique(a2$GID))  # 1020 (but 1529 observations)
length(unique(a1$GID))  # 1083 (but 2008 observations)
a2 <- log2(a2$metagenome.normalized)
length(a1) # 2008 observations
a1 <- log2(a1$metagenome.normalized)
length(a2) # 1529 observations
qqnorm(a2); qqline(a2)  # Approximately normal for most observations.
qqnorm(a1); qqline(a1)  # Similar though tails are less normal.
t.test(a2, a1, alternative="greater")                   # UNPAIRED, ONE-SIDED. Welch b/c unequal variance
wilcox.test(a2, a1, alternative='greater')
## For both tests, p ~ 0 so reject the null hypothesis.  Overall transcripts levels from
## UCYN-A2 are greater than from UCYN-A1.
rm(a1,a2)


################################################################################
##
## Results section 2: Highly transcribed genes were similar across sublineages
## and environments
##

## Check if UCYN-A3 genes with orthologs to the A1 and A2 cross-site detected
## genes tended to be more highly transcribed. This comparison is at ALOHA only
## because we only have A3 data at ALOHA.

gids <- crossSiteDetectedGids
nrow(subset(platTab[gids,], GENE!=''))      # 462 of the 760 ALOHA/SIO genes have gene symbol annotation.
nrow(subset(platTab[gids,], PATHWAY!=''))   # 495 have pathway annotation
## Now divide the A3 detected genes into those represented in the 760 vs. not.
gsyms <- names(table(platTab[gids,'GENE']))
gsyms <- gsyms[gsyms!='']
length(gsyms)  # 292
a3withHomologs <- subset(fData(aloha$eset), ORGANISM=='UCYNA3' & GENE %in% gsyms)
sort(table(a3withHomologs$PATHWAY), decreasing=T) 
a3withHomologs <- rownames(a3withHomologs)
length(a3withHomologs)  # 103
a3withoutHomologs <- rownames(subset(fData(aloha$eset), ORGANISM=='UCYNA3' & !(GENE %in% c('',gsyms))))
length(a3withoutHomologs)  # 44

## Note that exprs has log2 levels.  The TARA counts above were just counts.
a3withHomologs     <- apply(GetExprs(aloha$eset, a3withHomologs),    1, median)
a3withoutHomologs  <- apply(GetExprs(aloha$eset, a3withoutHomologs), 1, median)
qqnorm(a3withHomologs);    qqline(a3withHomologs)     # Check for approximate normality in order
qqnorm(a3withoutHomologs); qqline(a3withoutHomologs)  # to use the t-test.
t.test(a3withHomologs, a3withoutHomologs, alternative='greater')       # p = 3.1E-7; means 6.8 and 5.8; reject null
wilcox.test(a3withHomologs, a3withoutHomologs, alternative='greater')  # p = 2.4E-5; reject null
## Reject the null hypothesis.  At ALOHA, UCYN-A3 genes that had orthologs to
## A1/2 cross-site genes tended to be highly transcribed (-- higher than A3
## genes without orthologs to A1/2 cross-site genes).
rm(gids, gsyms, a3withHomologs, a3withoutHomologs)


## Test whether genes were transcribed with similar intensity at ALOHA and at
## Scripps.  Correlate the ranks of the intensities at both sites.  (Rather than
## the intensities because the data sets were normalized separately.)

CorrTest <- function(gids, pathway=NULL, corMeth='spearman')
{
    if (!is.null(pathway)) {
        gids <- rownames(subset(platTab[gids,], PATHWAY==pathway))
    } else { pathway <- 'any' }
    ## Limit to the 'gids' that were detected in both sites.
    both <- intersect(gids, intersect(rownames(aloha$eset),rownames(sio$eset)))
    alohaMeds <- apply(GetExprs(aloha$eset, both)[both,], 1, function(v) median(v))
    sioMeds   <- apply(GetExprs(sio$eset, both)[both,],   1, function(v) median(v))
    stopifnot(names(alohaMeds)==names(sioMeds))
    cat('Looking at',pathway,'pathway: ', length(both),'gids (of',length(gids),
        'gids) were detected at both sites.\n')
    cor.test(x=alohaMeds, y=sioMeds, method=corMeth)
}
gids <- rownames(subset(platTab, ORGANISM=='UCYNA1'))
CorrTest(gids)  # A1:  n = 368 genes; rho = 0.57;  p < 2.2e-16
gids <- rownames(subset(platTab, ORGANISM=='UCYNA2'))
CorrTest(gids)  # A2:  n = 392 genes; rho = 0.70;  p < 2.2e-16
rm(gids)


################################################################################
##
## Results section 3: Transcripts peaked at different times for Stn. ALOHA
## sublineages, except for genes involved in N2 fixation
##

## No statistical tests.  However, let's show how many genes changed peak times
## across habitats and how many did not.

GetCrossSitePeaksAndTroughs <- function(strain)
{
    x <- GetPeakAndTroughTimes(aloha, strain, F)  # All detected genes at ALOHA
    x$GID <- rownames(x)
    y <- GetPeakAndTroughTimes(sio,   strain, F)  # All detected genes at Scripps
    y$GID <- rownames(y)
    df <- merge(x, y, by='GID', all.x=F, suffixes = c('.aloha','.sio'))
    df <- as.data.frame(lapply(df, as.character))
    df
}

a1 <- GetCrossSitePeaksAndTroughs('UCYNA1')
a2 <- GetCrossSitePeaksAndTroughs('UCYNA2')

round(100 * sum(a1$peak.aloha != a1$peak.sio)  / nrow(a1), 1)   # 85.9% of a1 genes changed peaks
round(100 * sum(a1$peak.aloha == a1$peak.sio)  / nrow(a1), 1)   # 14.1% of a1 genes did not change peaks

round(100 * sum(a2$peak.aloha != a2$peak.sio)  / nrow(a2), 1)   # 82.1% of a1 genes changed peaks
round(100 * sum(a2$peak.aloha == a2$peak.sio)  / nrow(a2), 1)   # 17.9% of a1 genes did not change peaks

## Tables can be used to see which kinds of peak changes occurred (and
## non-changes are on the diagonal).
cat("UCYN-A1 genes counted by the time point of peak transcript level",
    "at Stn. ALOHA vs. Scripps Pier:")
table( data.frame(aloha = a1$peak.aloha, sio = a1$peak.sio) )

cat("UCYN-A2 genes counted by the time point of peak transcript level",
    "at Stn. ALOHA vs. Scripps Pier:")
table( data.frame(aloha = a2$peak.aloha, sio = a2$peak.sio) )

## For both UCYN-A1 and UCYN-A2, most genes peak at L9 (3pm) at Scripps Pier.
## Much more variable at Stn. ALOHA.  This is shown in Fig. 3.

rm(a1,a2)


################################################################################
##
## Results section 4 supplementary: Cross-site detected genes were less likely
## to have significant diel expression at Stn. ALOHA
##

## Contingency table of cross-site detected genes.
m <- matrix(c(
    aloha.diel    = length(intersect(rownames(aloha$cycGenes), crossSiteDetectedGids)),
    aloha.nonDiel = length(setdiff(crossSiteDetectedGids, rownames(aloha$cycGenes))),
    sio.diel      = length(intersect(rownames(sio$cycGenes), crossSiteDetectedGids)),
    sio.nonDiel   = length(setdiff(crossSiteDetectedGids, rownames(sio$cycGenes)))),
    nrow = 2, ncol=2,
    dimnames = list(c('diel','nonDiel'),c('ALOHA','Scripps')))

m
##         ALOHA Scripps
## diel       81     649
## nonDiel   679     111
##
## Strongly suggests that diel-ness differed across habitats.

fisher.test(m, alternative="two.sided")  # p ~ 0
## Reject the null hypothesis that there is no association between diel-ness and habitat.

rm(m)

## --- END ---
