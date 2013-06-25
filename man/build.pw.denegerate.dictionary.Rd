\name{build.pw.denegerate.dictionary}
\alias{build.pw.denegerate.dictionary}
\title{Build dictionary for positional degenerate wildcard kernel}
\usage{
  build.pw.denegerate.dictionary(kmer.len, mismatches,
    out.dir,
    alphabet = c("A", "C", "G", "T", "N", "S", "M", "R", "Y", "K", "W"),
    consecutive.mismatches = TRUE, no.cores = 1)
}
\arguments{
  \item{kmer.len}{Length of kmers}

  \item{mismatches}{Number of mismatches}

  \item{out.dir}{Directory to write the diciontary in}

  \item{alphabet}{Nucleotide alphabet. By defualt using all
  single/double nucleotide combinations and N Note that the
  first four alphabets have to be \code{c('A', 'C', 'G',
  'T')}}

  \item{consecutive.mismatches}{Logical indicating whether
  only consecutive mismatches should be allowed}

  \item{no.cores}{Number of cores of parallel processing}
}
\value{
  Creates a file containing the following elements:
  \code{pairwise.kmers} Distance matrix and
  \code{kmer.mapping} Mapping of reverse complements and
  kmers to features The filename is constructed based on
  kmer length and kernel used and has the structure:
  positional_wildcard_dict_kmer<\code{kmer.len}_mismatches\code{mismatches}_alpha\code{length
  (alphabet)}.Rdata
}
\description{
  This build a matrix of distances between all kmers and
  unique kmers, considering reverse complements and matches
  to degenerate nucleotide alphabet.
}

