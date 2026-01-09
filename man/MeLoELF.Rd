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
        mdir='~/',
        seq.file='sequences.txt',
        loc.file='methlocations.txt',
        call.file='methcalls.txt',
        crunch.too=T,
        process=T,
        read.length=c(10,3000),
        completeness=0.9,
        matching=0.9,
        target_fwd_auc=0.9,
        save.file='align.RData',
        load.file='align.RData',
        plot_title=NULL)
}
\arguments{
\item{parent}{A single character string indicating the substrate's FWD/parent strand sequence (5'->3'); DEFAULT -- none.}
\item{target}{A single character string indicating the substrate's REV/target strand sequence (5'->3'); DEFAULT -- none.}
\item{FWD.sites}{A numeric vector indicating indices for CpG methylation sites in the parent strand; DEFAULT -- none.}
\item{REV.sites}{A numeric vector indicating indices for CpG methylation sites in the target strand; DEFAULT -- none.}
\item{mdir}{A single character string indicating the working directory for files; DEFAULT -- base directory.}
\item{seq.file}{A single character string indicating the pre-processed filename for the txt file containing the Nanopore read sequences; DEFAULT -- 'sequences.txt'.}
\item{loc.file}{A single character string indicating the pre-processed filename for the txt file containing the Nanopore read CpG indices; DEFAULT -- 'methlocations.txt'.}
\item{call.file}{A single character string indicating the pre-processed filename for the txt file containing the Nanopore read methyl and hydroxy-methyl scores; DEFAULT -- 'methcalls.txt'.}
\item{crunch.too}{A logical indicating whether alignment data should be generated instead of loaded; DEFAULT -- TRUE.}
\item{process}{A logical indicating whether alignment data should be processed and analyzed; DEFAULT -- TRUE.}
\item{read.length}{A numeric vector of length 2 indicating minimum then maximum read length to process (>3000 base reads are likely to be data processing artifacts & <10 are likely to have little alignment accuracy); DEFAULT -- c(10,3000).}
\item{completeness}{A numeric value indicating the minimum proportion of reference sequence for a mapped fragment to be used in analysis (scales with substrate length, suggested ~0.95 for 40mer, ~0.8 for 250mer); DEFAULT -- 0.9.}
\item{matching}{A numeric value indicating the minimum base-call accuracy for a mapped fragment to be used in analysis (scales with substrate length, suggested ~0.95 for 40mer, ~0.8 for 250mer); DEFAULT -- 0.9.}
\item{target_fwd_auc}{A numeric value indicating the percentile confidence score cutoff below which methylation data is ignored; DEFAULT -- 0.9.}
\item{save.file}{A single character string indicating the name of the RData file to export mapping data to if crunch.too==T; DEFAULT -- 'align.RData'.}
\item{load.file}{A single character string indicating the name of the Rdata file to import mapping data from if crunch.too==F and process==T; DEFAULT -- 'align.RData'.}
\item{plot_title}{A single character string indicating the titles of outputted analysis graphs; DEFAULT -- 'NA'.}

}

\value{
Files:

\item{align.RData}{A compressed R data file containing all variables created as part of the pre-analysis fragment mapping and methylation score assignment.}
\item{methylation_matrices.RDS}{TBA}
\item{survival_summary.csv}{TBA}
\item{revbinary.csv}{TBA}

Graphs:

\item{methylation_survival.png}{TBA}
\item{methyl_distribution.png}{TBA}
}

\examples{

  # N/A

}
