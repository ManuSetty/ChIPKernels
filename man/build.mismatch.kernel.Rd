\name{build.mismatch.kernel}
\alias{build.mismatch.kernel}
\title{Build features and kernels}
\usage{
  build.mismatch.kernel(dictionary.file, seqs,
    kmers = NULL, kernel = TRUE)
}
\arguments{
  \item{dictionary.file}{File containing the dictionary
  built using \code{build.mismatch.dictionary} or
  \code{build.pw.denegerate.dictionary}.}

  \item{seqs}{DNAStringSet representing the sequences to
  build the mismatch kernel. None of the sequences should
  contain \code{N} in the nucleotide composition}

  \item{kmers}{Subset of kmers to build the kernel.
  \code{NULL} setting uses all kmers.}

  \item{kernel}{Boolean indicating if the kernel should be
  computed. If \code{FALSE}, only the feature
  representation is returned.}
}
\value{
  List with following elements: \code{features} Sparse
  feature matrix (seqs X kmers) and \code{kernel} Linear
  mismatch kernel.
}
\description{
  This builds the feature matrix and kernel using the
  dictionary. Linear kernel is computed using
  \code{kernlab} package.
}
\seealso{
  \code{\link{build.mismatch.dictionary}},
  \code{\link{kernelMatrix}}
}

