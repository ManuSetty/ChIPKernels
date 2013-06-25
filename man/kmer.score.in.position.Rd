\name{kmer.score.in.position}
\alias{kmer.score.in.position}
\title{kmer scores in all positions for sequences}
\usage{
  kmer.score.in.position(dictionary.file, seqs, kmer)
}
\arguments{
  \item{dictionary.file}{File containing the dictionary
  built using \code{build.mismatch.dictionary}.}

  \item{seqs}{DNAStringSet representing the sequences to
  build the mismatch kernel. None of the sequences should
  contain \code{N} in the nucleotide composition}

  \item{kmer}{kmer to be scores.}
}
\value{
  Matrix of scores for each sequence and position
}
\description{
  kmer scores in all positions for sequences
}
\seealso{
  \code{\link{build.mismatch.dictionary}},
  \code{\link{kernelMatrix}}
}

