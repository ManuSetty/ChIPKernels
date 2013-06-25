
# Revised           Comments
# 02/03/2013        Added kmer.score.in.position
# 02/21/2013        Updates to use create.kmers


#' Build mismatch dictionary
#'
#' @param kmer.len Desired length of the kmer
#' @param mismatches Number of mismatches to be allowed
#' @param out.dir Output directory
#' @param dimismatch Boolean indicating if dimismatches should be used. Default is \code{FALSE}
#' @description This build a matrix of distances between all kmers and unique kmers
#' (considering reverse complements).
#' @return Creates a file containing the following elements: \code{pairwise.kmers} Distance matrix
#' and \code{kmer.mapping} Mapping of reverse complements and kmers to features
#' The filename is constructed based on mismatches, kmer length and kernel used and has the structure:
#' \code{<di>}mismatch_dict_kmer<\code{kmer.len}_mismatches\code{mismatches}.Rdata
#' @export

build.mismatch.dictionary <- function (kmer.len, mismatches, out.dir, dimismatch = FALSE) {

  time.start <- get.time ()

  # Generate all kmers
  kmers.list <- create.kmers (kmer.len)
  kmer.rows <- kmers.list$all.kmers
  kmer.cols <- kmers.list$unique.kmers
  rkmer.cols <- kmers.list$unique.kmers.rc
  
  ## Dimismatch if necessary
  if (dimismatch) {
    orig.kmer.rows <- kmer.rows
    kmer.rows <- dinuc.alphabet (kmer.rows)
    kmer.cols <- dinuc.alphabet (kmer.cols)
    rkmer.cols <- dinuc.alphabet (rkmer.cols)
  }
  
  show ('Find distances ...')
  ## Find pairwise distances between rows and columns
  show ('Fwd kmers....')
  dist.fwd <- hammingWithMismatches (kmer.rows, kmer.cols, mismatches)
  show ('Rev kmers....')
  dist.rev <- hammingWithMismatches (kmer.rows, rkmer.cols, mismatches)

  show ('Set up matrix and save....')
  ## Find minimum mismatches for overlapping scores in forward and reverse
  rownames (dist.fwd) <- sprintf ("%d_%d", dist.fwd[,1], dist.fwd[,2])
  rownames (dist.rev) <- sprintf ("%d_%d", dist.rev[,1], dist.rev[,2])
  common <- intersect (rownames (dist.fwd), rownames (dist.rev))
  dist.fwd[common,3] <- pmin (dist.fwd[common,3], dist.rev[common,3])

  ## Merge
  dist.fwd <- rbind (dist.fwd, dist.rev[setdiff (rownames (dist.rev), rownames (dist.fwd)),])
  dist.fwd[,3] <- 1 - dist.fwd[,3]/(kmer.len-1)

  ## Sparse matrix of scores
  pairwise.kmers <- sparseMatrix(dist.fwd[,1], dist.fwd[,2], x=dist.fwd[,3],
                                 dims=c(length (kmer.rows), length (kmer.cols)),
                                 dimnames=list (as.character (kmer.rows), as.character (kmer.cols)))
  # kmer.column.mapping <- rep (as.character (kmer.cols), 2)
  # names (kmer.column.mapping) <- c(as.character (kmer.cols), as.character (rkmer.cols))

  out.file <- sprintf ("mismatch_dict_kmer%d_mismatches%d.Rdata", kmer.len, mismatches)
  if (dimismatch)
    out.file <- sprintf ("di%s", out.file)
  out.file <- sprintf ("%s/%s", out.dir, out.file)
  
  # Convert back if dimismatch
  if (dimismatch) {
    mapping <- as.character (orig.kmer.rows)
    names (mapping) <- as.character (kmer.rows)
    rownames (pairwise.kmers) <- mapping[rownames (pairwise.kmers)]
    colnames (pairwise.kmers) <- mapping[colnames (pairwise.kmers)]
  }

  # Mapping for columns
  kmer.column.mapping <- rep (colnames (pairwise.kmers), 2)
  names (kmer.column.mapping) <- c(colnames (pairwise.kmers), 
    as.character( reverseComplement (DNAStringSet (colnames (pairwise.kmers)))))

  save (pairwise.kmers, kmer.column.mapping,
        dimismatch, file=out.file)
  gc ()
  
  time.end <- get.time ()                        
  show (sprintf ("Total time for building dictionary: %.2f", (time.end - time.start)/60))
}



#' kmer scores in all positions for sequences
#'
#' @param dictionary.file File containing the dictionary built using \code{build.mismatch.dictionary}.
#' @param seqs DNAStringSet representing the sequences to build the mismatch kernel.
#' None of the sequences should contain \code{N} in the nucleotide composition
#' @param kmer kmer to be scores.
#' @return Matrix of scores for each sequence and position
#' @seealso \code{\link{build.mismatch.dictionary}}, \code{\link{kernelMatrix}}
#' @export

kmer.score.in.position <- function (dictionary.file, seqs, kmer) {

  if (length (grep ('N', seqs)) > 0 )
    stop ('Some of the sequences contain N - please remove these and re run.')

  if (length (unique (width (seqs))) != 1)
    stop ('All sequences have to be of the same width')

  # Load dictionary file
  load (dictionary.file)

  # Scores matrix
  kmer.len <- nchar (kmer)
  seq.len <- width (seqs[1])
  scores <- sparseMatrix (length (seqs), seq.len, x=0)

  # Score all positions
  for (i in 1:(seq.len - kmer.len + 1)) {
    sub <- as.character (substring (seqs, i, i + kmer.len -1))
    scores[,i] <- kmer.pairwise.features (pairwise.kmers, sub, kmer.column.mapping[kmer])
  }
  gc ()

  return (scores)
}

  



if (FALSE) {
    if (strand.info) {
    show ('Not collapsing forward and reverse complements...')
    ## Convert to distance
    dist.fwd[,3] <- 1 - dist.fwd[,3]/(kmer.len-1)
    dist.rev[,3] <- 1 - dist.rev[,3]/(kmer.len-1)
    
    ## Pool the data together
    pairwise.kmers <- sparseMatrix (c(dist.fwd[,1], dist.rev[,1]),
                                    c(dist.fwd[,2], dist.rev[,2] + length (kmer.cols)),
                                    x= c(dist.fwd[,3], dist.rev[,3]),
                                    dims=c(length (kmer.rows), length (kmer.cols) + length (rkmer.cols)),
                                    dimnames=list (as.character (kmer.rows),
                                      c(as.character (kmer.cols), as.character (rkmer.cols))))

    ## Remove duplicates
    pairwise.kmers <- pairwise.kmers[,!duplicated (colnames (pairwise.kmers))]
    
    out.file=sprintf ("mismatch_dict_strand_kmer%d_mismatches%d.Rdata", kmer.len, mismatches)
    if (dimismatch)
      out.file <- sprintf ("di%s", out.file)
    save (pairwise.kmers, dimismatch, file=sprintf ("%s/%s", out.dir, out.file))
    return (0)
  }


}