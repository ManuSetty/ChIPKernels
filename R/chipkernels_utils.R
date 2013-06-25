

## Revised        Comments
## 10/09/2012     Changing the way dimnames are set to NULL in kmer.pairwise.features
## 12/12/2012     Added hammingWithMismatches
## 02/21/2013     Added create.kmers


#' Convert to dinucleotide alphabet
#'
#' @param seqs DNAStringSet object of sequences to be converted to dinucleotide
#' @return BStringSet object of dinucleotides
#' @export

dinuc.alphabet <- function (seqs){

  ## Create mapping from nucleotides to dinucleotide alphabet
  Anew <- strsplit ('xbydefzhijklmnop---------', '')[[1]]
  nucs <- c('A', 'C', 'G', 'T')
  p <- permutations (4, 2, repeats.allowed=TRUE)
  dinucs <- paste (nucs[p[,1]], nucs[p[,2]], sep='')
  ## Add in alphabets with N nucleotide
  dinucs <- c(dinucs, "NN", sprintf ("N%s", nucs), sprintf ("%sN", nucs))
  Anew <- Anew[1:length (dinucs)]
  names (Anew) <- dinucs

  ## Convert DNA sequence to dinucleotide sequence
  seq.length <- width(seqs[1]) - 1
  dinuc.seqs <- rep ('', length (seqs))
  for (i in 1:seq.length) {
    show (sprintf ("Position %d",  i))
    dinuc.seqs <- paste (dinuc.seqs,
                             Anew[as.character (subseq (seqs, i, width=2))], sep='')
  }

  return (BStringSet (dinuc.seqs))
}


#' Subset sequences and kmer features
kmer.pairwise.features <- function (pairwise.kmers, seq.kmers, mapped.kmers) {
  scores <- pairwise.kmers[seq.kmers, mapped.kmers, drop=FALSE]
  dimnames (scores) <- list (NULL, NULL)
  return (scores)
}


#' hamming distance with mismatches
hammingWithMismatches <- function (x, y, max.nmis, fixed=TRUE, no.cores=1) {

  # Build function to apply
  toapply <- function (y.inds) {
    answer <- .Call ("XStringSet_dist_hamming_xy", x, y[y.inds], as.integer (max.nmis), fixed,
      PACKAGE="ChIPKernels")
    class (answer) <- 'numeric'
    answer <- matrix (answer, ncol=3, byrow=TRUE)
    return (answer)
  }

  if (no.cores == 1) {
    answer <- toapply (y.inds = 1:length (y))
  } else {
    y.list <- as.integer (seq (0, length (y), length.out=no.cores + 1))
    y.list <- lapply (2:length (y.list), function (x) {(y.list[x-1]+1):y.list[x]})
    res.list <- mclapply (y.list, FUN=toapply, mc.cores=no.cores)

    # Merge results from multiple threads
    answer <- c()
    cs <- c(0, cumsum (sapply (y.list, length)))
    for (i in 1:no.cores) {
      r <- res.list[[i]]
      r[,2] <- r[,2] + cs[i]
      answer <- rbind (answer, r)
    }
  }
  return (answer)
}


#' Function to create all kmers
create.kmers <- function (kmer.len) {

  show ('Generating kmers...')
  nucs <- c('A', 'C', 'G', 'T')
  
  ## Define all possible kmers
  kmers.mat <- permutations (length (nucs), kmer.len, nucs, repeats.allowed=TRUE)
  kmers <- rep ("", nrow (kmers.mat))
  for (i in 1:ncol (kmers.mat))
    kmers <- sprintf ("%s%s", kmers, kmers.mat[,i])
  kmers <- DNAStringSet (kmers)

  # Unique kmers by eliminating reverse complements
  tmp <- cbind (as.character (kmers), as.character (reverseComplement (kmers)))
  tmp[tmp[,1] > tmp[,2],1] <- tmp[tmp[,1] > tmp[,2],2]
  unique.kmers <- DNAStringSet (unique (tmp[,1]))
  unique.kmers.rc <- reverseComplement (unique.kmers)

  # REturn
  return (list (all.kmers = kmers,
    unique.kmers = unique.kmers,
    unique.kmers.rc = unique.kmers.rc))

}

