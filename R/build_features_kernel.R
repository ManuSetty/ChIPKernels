

#' Build features and kernels
#'
#' @param dictionary.file File containing the dictionary built using \code{build.mismatch.dictionary} or \code{build.pw.denegerate.dictionary}.
#' @param seqs DNAStringSet representing the sequences to build the mismatch kernel.
#' None of the sequences should contain \code{N} in the nucleotide composition
#' @param kmers Subset of kmers to build the kernel. \code{NULL} setting uses all kmers.
#' @param kernel Boolean indicating if the kernel should be computed. If \code{FALSE}, only the feature
#' representation is returned.
#' @description This builds the feature matrix and kernel using the dictionary. Linear kernel is computed
#' using \code{kernlab} package.
#' @return List with following elements: \code{features} Sparse feature matrix (seqs X kmers)
#' and \code{kernel} Linear mismatch kernel.
#' @seealso \code{\link{build.mismatch.dictionary}}, \code{\link{kernelMatrix}}
#' @export

build.features.kernels <- function (dictionary.file, seqs,
                                   kmers=NULL, kernel=TRUE){

  if (length (grep ('N', seqs)) > 0 )
    stop ('Some of the sequences contain N - please remove these and re run.')

  block.size <- 20
  
  ## Load dictionary
  load (dictionary.file)

  ## Set up
  if (is.null (kmers))
    kmers <- colnames (pairwise.kmers)
  kmer.len <- nchar (kmers[1])
  seq.len <- nchar (seqs[1])
  
  ## Determine distances
  features <- sparseMatrix (length (seqs), length (kmers), x=0)
  features.iteration <- features

  time.start <- get.time ()
  for (i in 1:(seq.len - kmer.len + 1)) {
    if (i %% block.size == 1) {
      show (sprintf ("Position %d of %d", i, seq.len - kmer.len + 1))
      features <- features + features.iteration
      features.iteration <- sparseMatrix (length (seqs), length (kmers), x=0)
    }

    sub <- as.character (substring (seqs, i, i + kmer.len -1))
    features.iteration <- features.iteration + kmer.pairwise.features (pairwise.kmers, sub, kmer.column.mapping[kmers])
    gc ()
    
  }
  features <- features + features.iteration
  time.end <- get.time ()
  show (sprintf ("Time for determining features : %.2f", (time.end - time.start) / 60 ))

  results <- list ()
  results$features <- features
  colnames (results$features) <- kmers
  
  ## Normalize features
  if (kernel) {
    
    n.factors <- sqrt (rowSums (features ^ 2)) + 1e-10
    features <- features / n.factors
    
    ## Kernel
    kernel <- kernelMatrix(polydot (d=1), as.matrix (features)) - 1
    diag (kernel) <- 0
    results$kernel <- kernel
    dimnames (results$kernel) <- list (kmers, kmers)
  }

  return (results)
}
