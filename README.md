## UCYN-A Metatranscriptomes at Stn. ALOHA and Scripps Pier

Key scripts and data files for
**Open ocean and coastal strains of the N<sub>2</sub>-fixing cyanobacterium UCYN-A have distinct transcriptomes**
by María del Carmen Muñoz-Marín, Jonathan D. Magasin, and Jonathan P. Zehr
under review at _PLOS ONE_.

Biologically usable forms of nitrogen are scarce over much of the world ocean. Thus, marine microbes that can fix atmospheric N<sub>2</sub> gas into usable ammonia (NH<sub>3</sub>) are critical to marine food webs. One such microbe, the N<sub>2</sub>-fixing, the endosymbiotic cyanobacterium _Candidatus_ Atelocyanobacterium thalassa ("UCYN-A"), can contribute substantial amounts of fixed nitrogen to marine microbial communities.  UCYN-A genomes are highly strealined and lack the full complement of genes required for obtaining carbon through photosynthesis.  Instead UCYN-A obtains carbon from its haptophyte (photosynthetic) host in exchange for fixed nitrogen.

Our study compared metatranscriptomes from natural populations of UCYN-A that were collected from two very different habitats, the nutrient-poor North Pacific Subtropical Gyre (at Stn. ALOHA) and the nutrient-rich waters off the California coast (at Scripps Pier).  For our study, metatranscritomes were obtained for samples collected over ~2-days with 3 to 6 hours between samples in order to study day-night gene expression patterns.  Differences in expression patterns across habitats can improve our understanding of UCYN-A ecotypes and their hosts, which differ among sublineages.

The Scripps Pier metatranscriptomes were previously reported by [Muñoz-Marín et al. (2019)](https://journals.asm.org/doi/10.1128/mBio.02495-18).  The Stn. ALOHA data are new.  Both studies used custom UCYN-A microarrays that targeted major UCYN-A ecotypes UCYN-A1 (open ocean) and UCYN-A2 (coastal).  The Stn. ALOHA data used an updated microarray design that included known genes for UCYN-A3 (open ocean).


### Repository contents:
* scripts/
    - runMicroTOOLsPipeine.R:  Runs the [MicroTOOLs microarray pipeline](https://www.jzehrlab.com/microtools) on the Stn. ALOHA arrays.  This is provided to illustrate the major steps for processing the arrays from quality control through robust multi-array averaging and the identification of genes with significant 24 h periodic expression.  Main outputs from running the pipeline are provided in Data/workspace.Rdata.
    - findCyclicGenes.R:  This is run by the previous script to identify genes with significant 24 h periodic expression ("diel genes").
    - statistical_tests.R:  Interactive script that does statistical tests that appear in our publication.  Relies on Data/workspace.Rdata.
    
* Metadata/:  Includes the main configuration file that was used to run runMicroTOOLsPipeline.R and two other files that describe the Stn. ALOHA data. Note that the Stn. ALOHA data reside at the [NCBI Gene Expression Omnibus](https://www.ncbi.nlm.nih.gov/geo/) under accession GSE206403.

* Data/workspace.Rdata:  R image that has the Stn. ALOHA and Scripps Pier ExpressionSets and other key objects used during the analysis.

* Logs/:  Includes the log files from runMicroTOOLsPipeline.R and findCyclicGenes.R that are associated with the final data analyzed for the publication.
