# Carotid.Art.Dis.scRNA_Seq
## Code Files For Analyzing 3 Patient Carotid Artery Single Cell RNA-Seq Experiment

The code files in this repository were created and used to analyze three pairs of patient tissue samples collected during carotid endocardectomy procedures. For more details of the experiment, please see our bioRxiv paper:

### *Decoding the transcriptome of atherosclerotic plaque at single-cell resolution*
Tom Alsaigh, Doug Evans, David Frankel, Ali Torkamani
[bioRxiv 2020.03.03.968123; doi: ](https://doi.org/10.1101/2020.03.03.968123) with [GEO Data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159677)

It is expected that this archive version of the paper will be updated and formally published in a recognized journal (after 11/01/2021). Interested parties should look for it's formal publication for the authoritive presentation of the experiment.

## Analysis Process
The data used in this analysis was taken from 3 pairs of samples that were collected at different times and processed as independant pairs. Within each pair are samples of the carotid artery where an athersclorotic lession had formed (Atheroscloric Core or AC) and nearby carotid tissue relatively devoid of plaque (Proximal Adjact or PA). AC and PA sample pairs were processed using a 10x Genomics single cell platform. The sample pairs were NGS processed as multiplexed pairs.

The resulting NGS data from these experiments were processed using the 10x Genomics Cellranger suite of software. Those files (NGS and Cellranger) serve as the primary input to the analytic programs found in this repository. These source data files are available via [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159677). When running the 10x software, default parameters and the latest reference transcriptome were always used, so no documentation is included here showing the Cellranger commands. The NGS input files for Cellranger and found on GEO should be sufficient for an interested researcher to reproduce the Cellranger results, also found on GEO.

## Analysis Steps


