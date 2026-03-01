\name{MeLoELF}
\alias{MeLoELF}
\title{Analysis of MeLoELF-seq Data}
\usage{
MeLoELF(parent,
        target,
        FWD.sites,
        REV.sites,
        mdir=getwd(),
        sam.file='auto',
        sam.indices=c(10,33,34,0),
        crunch.too=T,
        process=T,
        pre.ligated=F,
        methyl.type='B',
        thresh.meth='BM',
        BM.strand.data='r',
        var.thresh=0,
        T1.err=0.01,
        T2.err=0.05,
        target_fwd_auc=0.9,
        read.length=NULL,
        padding=2,
        completeness=0.3,
        matching=0.9,
        exact.search=F,
        align.file='align.RData',
        processed.file='processed.RData',
        seq.file='sequences.txt',
        melo.file='methlocations.txt',
        meca.file='methcalls.txt',
        plot_title='',
        plot.nom=c('FWD','REV'))
}
\arguments{

Required:

\item{parent}{A single character string indicating the substrate's FWD/parent strand sequence (5'->3'); DEFAULT -- none.}
\item{target}{A single character string indicating the substrate's REV/target strand sequence (5'->3'); DEFAULT -- none.}
\item{FWD.sites}{A numeric vector indicating indices for CpG methylation sites in the parent strand; DEFAULT -- none.}
\item{REV.sites}{A numeric vector indicating indices for CpG methylation sites in the target strand; DEFAULT -- none.}

Defaults:

\item{mdir}{A single character string indicating the working directory for files; DEFAULT -- working directory.}
\item{sam.file}{A single character string indicating the filename for the input sam file containing the Nanopore sequencing data; DEFAULT -- searches mdir directory for sam file.}
\item{sam.indices}{A numeric vector of length four indicating the SAM file's column numbers for sequences, methyl locations, and methyl scores, respectively, as well as the number of sam file header lines to skip; DEFAULT -- c(10,33,34,0)}
\item{crunch.too}{A logical indicating whether alignment data should be generated instead of loaded; DEFAULT -- TRUE.}
\item{process}{A logical indicating whether alignment data should be processed and analyzed; DEFAULT -- TRUE.}
\item{pre.ligated}{A logical indicating whether the substrate was ligated before enzyme activity (see Details before using); DEFAULT -- FALSE.}
\item{methyl.type}{A single character string indicating which methyl scores should be used for analysis -- 'M' = methyl, 'H' = hydroxy-methyl, 'B' = both/combined; DEFAULT -- 'B'.}
\item{thresh.meth}{A single character string indicating which methyl score thresholding method should be used, or a single numeric value indicating the desired value for manual thresholding -- 'BM' = Mixed Beta Distributions, 'AUC' = Survival AUC Tuning, 0-1 = manual thresholding; DEFAULT -- 'BM'.}
\item{BM.strand.data}{A single character string indicating what subset of the data's CpG methyl scores should be used for thresholding if thresh.meth='BM' and var.thresh<2 (ignored for thresholding purposes if var.thresh=2) -- 'f' = FWD/parent strand scores, 'r' = REV/target strand scores, 'b' = all scores; DEFAULT -- 'r'.}
\item{var.thresh}{A numeric value indicating whether CpG sites should be thresholded individually (ignored unless pre.ligated=F and thresh.meth='BM') -- 0 = no individual site thresholding, 1 = individual thresholding for rev/target strand sites only, 2 = individual thresholding for all sites; DEFAULT -- 0.}
\item{T1.err}{A numeric value indicating what type-I error rate should be targeted during mixed beta thresholding if thresh.meth='BM' -- 0-1; DEFAULT -- 0.01.}
\item{T2.err}{A numeric value indicating what type-II error rate should be targeted during mixed beta thresholding if type-I error rate cannot be estimated and if thresh.meth='BM' -- 0-1; DEFAULT -- 0.05.}
\item{target_fwd_auc}{A numeric value indicating the percentile confidence score cutoff below which methylation data is ignored if thresh.meth='AUC' -- 0-1; DEFAULT -- 0.9.}
\item{read.length}{A numeric vector of length 2 indicating minimum then maximum read length to process (>3000 base reads are likely to be data processing artifacts & <10 are likely to have little alignment accuracy) -- 0-Inf; DEFAULT -- scales lower bound to half the REV/target sequence length, and upper bound to the 98th percentile of SAM file read lengths.}
\item{padding}{A numeric value indicating the number of additional fragments per read to attempt extracting -- 0-Inf; DEFAULT -- 2.}
\item{completeness}{A numeric value indicating the minimum length of mapped fragments relative to the reference sequence for them to be used in analysis -- 0-1; DEFAULT -- 0.3.}
\item{matching}{A numeric value indicating the minimum proportion of matching basecalls for a mapped fragment to be used in analysis -- 0-1; DEFAULT -- 0.9.}
\item{exact.search}{A logical indicating whether the default fragment mapping/alignment should be replaced with a faster, simpler motif search (see Details before using); DEFAULT -- F.}
\item{align.file}{A single character string indicating the name of the RData file to export mapping data to if crunch.too==T, or to import mapping data from if crunch.too==F and process==T; DEFAULT -- 'align.RData'.}
\item{processed.file}{A single character string indicating the name of the Rdata file to export processed data to if process==T; DEFAULT -- 'processed.RData'.}
\item{seq.file}{A single character string indicating the filename for the exported txt file containing the Nanopore read sequences; DEFAULT -- 'sequences.txt'.}
\item{melo.file}{A single character string indicating the filename for the exported txt file containing the Nanopore read CpG indices; DEFAULT -- 'methlocations.txt'.}
\item{meca.file}{A single character string indicating the filename for the exported txt file containing the Nanopore read methyl and hydroxy-methyl scores; DEFAULT -- 'methcalls.txt'.}
\item{plot_title}{A single character string indicating a common ID to append to the titles of outputted analysis graphs; DEFAULT -- none.}
\item{plot.nom}{A character string vector of length 2 indicating the naming convention to use for the FWD/parent and REV/target strands, respectively, in outputted analysis graphs; DEFAULT -- c('FWD','REV').}
}

\details{

{sam.file} -- {If sam.file is set to NULL, or if no sam file is found by auto-search, then the working directory will be searched to directly load seq.file, melo.file, and meca.file.}

{thresh.meth} -- {The BM option fits the collective methyl scores across FWD and REV strands to a mixed beta distribution, then sets the threshold based on a 0.95 (by default) quantile in the low/unmethylated distribution. The AUC option scales the threshold based on what produces an FWD/parent strand AUC~0.9 (by default) in the survival analysis; it is only suited for hemimethylated substrate reactions, where the FWD strand is available as a positive control with known methylation status.}

{var.thresh} -- {If var.thresh=1, the FWD/parent strand will be thresholded based on the composite distribution of CpG sites implicated by BM.strand.data.}

{exact.search} -- {This parameter alters the mapping algorithm to a simple motif search, which drastically shortens fragment extraction time by only scanning each read for exact matches to the FWD/parent and REV/target sequences. The format of almost all data structures and output are conserved, with the exception of fragment quality filtering-related metrics. Relatedly, unless used in conjunction with pre.ligated=T, using this setting deprecates the completeness and matching parameters.}

{pre.ligated} -- {This parameter is intended for niche applications with enzymes acting on hemi-methylated substrate. If pre.ligated=T, mapped fragments are computationally reconstructed into their original read polymers, sorted into FWD and REV homopolymers, filtered according to their read mapping and coverage quality, then collapsed into vectors of only their indexed CpG sites for analysis. This option ignores the var.thresh parameter entirely, and it repurposes the completeness and matching parameters during polymer filtering for read coverage and mapped section basematch accuracy metrics, respectively. In addition, only the methylation_survival plot and a rudimentary version of the QCgraphs plot are exported.}

---------
---------

Summary of MeLoELF-seq Analysis Pipeline

======

STEP 1:	Basecalling

POD5 files containing raw signal data are basecalled via ONTs Dorado v1.4.0 software.

user$ ~/dorado-1.4.0-linux-x64/bin/dorado basecaller sup,5mCG_5hmCG pod5/ --kit-name SQK-NBD114-24 > bam/combined.bam


STEP 2:	Subsetting

The combined BAM output file is sorted into individual sequencing sample BAM files via barcodes with Dorado.

user$ ~/dorado-1.4.0-linux-x64/bin/dorado demux -o bam/ --no-classify bam/combined.bam


STEP 3:	BAM-to-SAM File Conversion

Individual BAM files are converted to SAM files with Samtools.

user$ samtools view file.bam > file.sam


STEP 4:	MeLoELF Analysis

GitHub Repository:	www.github.com/whemphil/MeLoELF

MeLoELF is an R package for parsing the sequence, 5mC/5hmC score, and 5mC/5hmC location data provided within Dorado-basecalled SAM files, and analyzing the methylation states of the individual substrate molecules that compose the sequenced polymers.
MeLoELF was coded in R v4.5.2 on MacOS, with Unix machines in mind. Function parameters and data structures were designed for optimal analysis of enzymes acting on hemimethylated substrates, like DNMT1, but are generally applicable for de novo methylases as well.
MeLoELF analysis can be divided into five general phases.

Phase 1:  SAM File Parsing

The SAM file (identified by the mdir and sam.file parameters) is parsed with “awk” via system/bash invocation to extract the data columns corresponding to read sequences, MM:Z methyl locations/indices, and ML:B:C 5mC/5mhC scores (identified by the sam.indices parameter). For each read, the sequence information is converted to a character vector, and methyl and hydroxy-methyl score vectors of equal lengths are prepared. Default scores are set to zero for “C” basecalls and NA otherwise, then base indices implicated by the MM:Z values are updated with the corresponding ML:B:C scores. Finally, the ML:B:C scores, which are reported on a 0-256 integer scale in SAM files, are normalized to a 0-1 scale.

These three aligned vectors for each read are used in phase 2.

Phase 2:  Fragment Mapping and Extraction

The original, complementary substrate sequences (identified by parent and target parameters) are each slid across the entire read sequence, and the fraction of matching bases calculated for every possible alignment position with each reference sequence. In order of decreasing reference-read congruency, fragments are “extracted” from the read and recorded into data frames, along with their corresponding methyl score information, position within the read, and read number. This continues with recurring sliding alignments until the entire read sequence has been “extracted”, or the number of extracted fragments reaches a padded (set by padding parameter) estimate of how many substrate repeats the oligomeric read should contain.

For cases of partial read-reference congruency, fragments are “trimmed” to the most 5’ and 3’ base-matches that have a directly adjacent base match, and only information within those indices is “extracted” from the parent read. Consequently, fragments with sequence insertions/deletions are initially recorded as separate, partial fragments in the data frame, though their read position data allows for recombination later.

For each extracted fragment, two quality parameters are calculated. Completeness, which indicates the length of each extracted fragment relative to the reference sequence length, and matching, which indicates the proportion of bases within the extracted fragment that match their positional counterparts in the reference sequence.

In our MeLoELF experiments, product ligations produce oligomers with an average of ~6 fragments (~240 bases). However, in a fraction of a percent of SAM file entries, erroneous reads appear that can be thousands of times longer than the median read length. To avoid unnecessarily inflating computational time, and because these reads are a minority that likely represent sequencing artifacts rather than physical polymer sizes, reads outside a certain size range (set by the read.length parameter) are skipped entirely.

Phase 3:  Data Consolidation and Fragment Quality Filtering

Data are re-organized into consolidated data frames containing extracted fragment information across all reads. Initially, data frames contain information of all “extracted” fragments, no matter how poor their congruency to the reference sequence may be. Before subsequent analyses, extracted fragments may be filtered out according to their completeness and matching quality metrics from phase 2 (set by parameters of the same name). In addition, their ML:B:C data are pruned to only include 5mC scores, 5hmC scores, or the sum of those scores (set by the methyl.type parameter).

Methylation score data frames are column-pruned to isolate only the CpG sites of interest (set by FWD.sites and REV.sites parameters), and REV/target fragment data frames are column-reversed to complement their FWD counterparts.

Phase 4:  Methyl Score Thresholding

MeLoELF offers three approaches for binarizing methyl score data (set by the thresh.meth parameter) into C vs. 5(h)mC predictions: Area Under Curve (AUC) optimization, Beta unmixing (BM), and manual selection. BM is the MeLoELF default.

AUC thresholding operates under the assumption that the reaction substrate pre-sequencing was hemi-methylated, with the FWD/parent fragments having a known, pre-existing methylation pattern. It scales a single threshold value down, iteratively checking the subsequent survival analysis’s AUC metric (see phase 5) until it reaches a certain value on the FWD fragments (set by the target_fwd_auc parameter; 0.9 by default). That threshold is then applied to all sites.

BM thresholding works by separating (or pooling entirely) methyl scores into comparable groups based on their substrate strand and/or CpG site (set by the BM.strand.data and var.thresh parameters), then fitting the methyl score distributions with mixed beta probability density functions. It operates under the assumption that the lower and upper distributions correspond to negative/unmethylated and positive/methylated sites, respectively.

Once the data have been fit to the mixed beta distribution, the lower/unmethylated distribution is isolated and a threshold set that aims to control for type-I error (set by T1.err parameter).

In the event that a methyl score distribution is found to be best fit by a single beta distribution, its identity is inferred from the distribution’s mode if it’s sufficiently extreme, or from the relative position of unexplained distribution density in the case of more ambiguous distributions. For single, negative distributions, thresholding is set by controlling for type-I error as described above. For single, positive distributions, thresholding is instead set by controlling for type-II error rate (set by the T2.err parameter).

In the case of a manual thresholding override, the methyl scores of all sites and strands are binarized from a single, provided value.

Phase 5:  Methylation State Analyses

The first analysis output file, QCgraphs, contains the plots “Methyl Score Distributions” and “Methyl Location Distributions”. The former provides violin plots of the normalized ML:B:C scores (set by methyl.type parameter) at each CpG sites along the aligned fragments, as well as counts of the of scored bases mapped to each site. Data are from only fragments that passed quality filtering. The latter plot provides bar graphs of fractions of FWD and REV fragments with methylation at each CpG site (according to selected thresholding), as well as the fraction of fragments with those CpG sites as their most 5’ or 3’ methylations. In addition, there are exact values for the percent of FWD and REV CpG sites methylated (“CpG”), FWD and REV fragments methylated (“Sub.”), FWD and REV fragments passing quality filtering (“Qual.”), and the total coverage of mapped reads before and after fragment quality filtering ("Map").

The second analysis output file, QCgraphsB, contains distribution plots of mapped fragment positions (that pass quality filtering). The metrics for fragment position are “Relative to Read”, which plots location of the fragment’s center relative to its parent read length, “From Nearest Edge”, which plots how many multiples of the reference sequence length the mapped fragment’s center is from its parent read’s closest 5’ or 3’ end, “From 5’/3’ End”, which plot the same metric but from a specific parent read end, and “From Read Center”, which plots how many multiples of the reference sequence length the mapped fragment’s center is from the center of its parent read. Probability density functions are provided for both the combined and separate FWD and REV fragment positional data.

The third analysis output file, methylation_survival, provides survival plots across methylated fragment data for a metric of enzyme processivity. In brief, the 5’-most methylation on each methylated fragment is identified, then the number of consecutive methylations 3’ of the site divided by the number of CpG sites 3’ of the site to give a survival score. Survival scores are subsetted by FWD and REV fragments, and survival curves of those scores integrated to give the AUC metric. The “Sim. Dist.” plot corresponds to a randomization of the binarized REV fragment data frame subjected to the same analysis, meant to provide a baseline for the AUC of a distributive reaction with similar overall methylation levels.

The fourth analysis output file, methyl_distribution, provides histograms of the number of methylated CpG sites across FWD and REV fragments. The “Sim. Dist.” corresponds to the same control as above.

}

\value{

Files:

\item{sequences.txt}{A txt file containing the Nanopore read sequences parsed from the input sam file.}
\item{methlocations.txt}{A txt file containing the Nanopore read C-methylation indices parsed from the input sam file.}
\item{methcalls.txt}{A txt file containing the Nanopore read C-methyl scores parsed from the input sam file.}
\item{align.RData}{A compressed R data file containing all relevant variables created as part of the pre-analysis fragment mapping and methylation score assignment.}
\item{processed.RData}{A compressed R data file containing the processed matrices of FWD and REV mapped fragments' methyl scores that passed quality filtering.}
\item{methylation_matrices.RDS}{A compressed R data file containing the selected threshold value for methyl score binarization, and containing the data frames with survival data for the FWD, REV, and simulated control fragments.}
\item{survival_summary.csv}{A csv file providing brief summary statistics for the survival data for the FWD, REV, and simulated control fragment data, such as their respective survival plot AUCs.}
\item{revbinary.csv}{A csv file containing the binarized matrix of mapped REV fragments' methylation states.}

Graphs:

\item{methylation_survival.png}{Survival plots for the FWD, REV, and simulated control methylation state data.}
\item{methyl_distribution.png}{Grouped histograms of the numbers of methyl marks on each mapped fragment of the FWD, REV, and simulated control data sets.}
\item{QCgraphs.png}{A collection of plots with quality control information about the methylation properties of FWD and REV mapped fragments. Includes violin plots of methyl score distributions for each methylation site relative to the selected threshold, histograms of fraction methylation at individual CpG sites and of initial REV site methylation assuming directional processivity, and calculations of percent CpG site methylation, percent substrate molecule methylation, percent of mapped fragments passing quality control filtering, and percent coverage of mapped reads before and after fragment quality filtering.}
\item{QCgraphsB.png}{A collection of plots with quality control information about the mapping biases of FWD and REV fragments. Includes probability density distributions for several metrics of mapped fragment positions on their parent reads.}
\item{R Graphical Output}{If the BM thresholding method is used, a plot comparing the probability densities of the observed methyl score data and the fitted mixed beta distribution is produced. If var.thresh>0, then multiple plots will be outputted; first for all the combined scores of strands indicate by BM.strand.data, then for each CpG site starting from the 5' end of the fwd/parent strand and alternating between the rev/target and fwd/parent strands if var.thresh=2.}

}

\description{
Function to map and computationally deconstruct end-ligated fragment sequences, parse their methylation states, and analyze enzyme activity and processivity.

Designed and optimized for Unix machines -- may return errors on Windows PCs.
}



