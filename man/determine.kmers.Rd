\name{determine.kmers}
\alias{determine.kmers}
\title{Determine all kmers with mismatches of a given length}
\usage{
  determine.kmers(alphabet, mismatches, kmer.len,
    range = NULL, consecutive.mismatches)
}
\arguments{
  \item{alphabet}{Nucleotide alphabet}

  \item{mismatches}{Number of mismatches allowed}

  \item{range}{Range of numbers to determine kmers}

  \item{kmer.len}{Length of kmers}
}
\value{
  A list of DNAStringSet elements: \code{kmers}
  containining kmers without mismatches and
  \code{kmers.mismatches} containing all the kmers that
  were generated'
}
\description{
  Builds all possible kmers in the given \code{range} for
  the given \code{alphabet}
}

