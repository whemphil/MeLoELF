\name{MeLoELF}
\alias{MeLoELF}
\title{Analysis of MeLoELF-seq Data}
\description{
Function to map and computationally deconstruct end-ligated frament sequences, parse their methylation states, and analyze enzyme activity and processivity.
}
\usage{
MeLoELF(parent,
        target,
        FWD.sites,
        REV.sites,
        mdir=getwd(),
        sam.file=list.files(path = mdir,pattern = '*.sam'),
        crunch.too=T,
        process=T,
        pre.ligated=F,
        methyl.type='B',
        thresh.meth='BM',
        BM.strand.data='rev',
        target_fwd_auc=0.9,
        read.length=c(10,3000),
        completeness=0.9,
        matching=0.9,
        align.file='align.RData',
        processed.file='processed.RData',
        seq.file='sequences.txt',
        melo.file='methlocations.txt',
        meca.file='methcalls.txt',
        plot_title='')
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
\item{crunch.too}{A logical indicating whether alignment data should be generated instead of loaded; DEFAULT -- TRUE.}
\item{process}{A logical indicating whether alignment data should be processed and analyzed; DEFAULT -- TRUE.}
\item{pre.ligated}{A logical indicating whether the substrate was ligated before enzyme activity; DEFAULT -- FALSE.}
\item{methyl.type}{A single character string indicating which methyl scores should be used for analysis -- 'M' = methyl, 'H' = hydroxy-methyl, 'B' = both/combined; DEFAULT -- 'B'.}
\item{thresh.meth}{A single character string indicating which methyl score thresholding method should be used -- 'BM' = Mixed Beta Distributions, 'AUC' = Survival AUC Tuning; DEFAULT -- 'BM'.}
\item{BM.strand.data}{A single character string indicating what subset of the data's CpG methyl scores should be used for thresholding if thresh.meth='BM' -- 'f' = FWD/parent strand scores, 'r' = REV/target strand scores, 'b' = all scores; DEFAULT -- 'r'.}
\item{target_fwd_auc}{A numeric value indicating the percentile confidence score cutoff below which methylation data is ignored if thresh.meth='AUC' -- 0-1; DEFAULT -- 0.9.}
\item{read.length}{A numeric vector of length 2 indicating minimum then maximum read length to process (>3000 base reads are likely to be data processing artifacts & <10 are likely to have little alignment accuracy) -- 0-Inf; DEFAULT -- c(10,3000).}
\item{completeness}{A numeric value indicating the minimum proportion of reference sequence length for a mapped fragment to be used in analysis -- 0-1; DEFAULT -- 0.9.}
\item{matching}{A numeric value indicating the minimum proportion of matching basecalls for a mapped fragment to be used in analysis -- 0-1; DEFAULT -- 0.9.}
\item{align.file}{A single character string indicating the name of the RData file to export mapping data to if crunch.too==T, or to import mapping data from if crunch.too==F and process==T; DEFAULT -- 'align.RData'.}
\item{processed.file}{A single character string indicating the name of the Rdata file to export processed data to if process==T; DEFAULT -- 'processed.RData'.}
\item{seq.file}{A single character string indicating the filename for the exported txt file containing the Nanopore read sequences; DEFAULT -- 'sequences.txt'.}
\item{melo.file}{A single character string indicating the filename for the exported txt file containing the Nanopore read CpG indices; DEFAULT -- 'methlocations.txt'.}
\item{meca.file}{A single character string indicating the filename for the exported txt file containing the Nanopore read methyl and hydroxy-methyl scores; DEFAULT -- 'methcalls.txt'.}
\item{plot_title}{A single character string indicating the common ID for outputted analysis graphs; DEFAULT -- none.}
}

\details{

{thresh.meth} -- {The BM option fits the collective methyl scores across FWD and REV strands to a mixed beta distribution, then sets the threshold based on a 0.95 (by default) quantile in the low/unmethylated distribution. The AUC option scales the threshold based on what produces an FWD/parent strand AUC~0.9 (by default) in the survival analysis; it is only suited for hemimethylated substrate reactions, where the FWD strand is available as a positive control with known methylation status.}

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
\item{QCgraphs.png}{A collection of plots with quality control information about the methylation properties of FWD and REV mapped fragments. Includes violin plots of methyl score distributions for each methylation site relative to the selected threshold, histograms of fraction methylation at individual CpG sites and of initial REV site methylation assuming directional processivity, and calculations of percent CpG site methylation, percent REV substrate molecule methylation, and percent of mapped REV fragments passing quality control filtering.}
\item{R Graphical Output}{If the BM thresholding method is used, a plot comparing the probability densities of the observed methyl score data and the fitted mixed beta distribution is produced.}

}

