\name{build.mismatch.dictionary}
\alias{build.mismatch.dictionary}
\title{Build mismatch dictionary}
\usage{
  build.mismatch.dictionary(kmer.len, mismatches, out.dir,
    dimismatch = FALSE)
}
\arguments{
  \item{kmer.len}{Desired length of the kmer}

  \item{mismatches}{Number of mismatches to be allowed}

  \item{out.dir}{Output directory}

  \item{dimismatch}{Boolean indicating if dimismatches
  should be used. Default is \code{FALSE}}
}
\value{
  Creates a file containing the following elements:
  \code{pairwise.kmers} Distance matrix and
  \code{kmer.mapping} Mapping of reverse complements and
  kmers to features The filename is constructed based on
  mismatches, kmer length and kernel used and has the
  structure:
  \code{<di>}mismatch_dict_kmer<\code{kmer.len}_mismatches\code{mismatches}.Rdata
}
\description{
  This build a matrix of distances between all kmers and
  unique kmers (considering reverse complements).
}

