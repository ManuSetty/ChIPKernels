
# Script for creating positional wildcard kernel with degenerate nucleotide representation
# Manu Setty
# 06/21/2013

#' Build dictionary for wildcard kernel
#'
#' @param kmer.len Length of kmers
#' @param mismatches Number of mismatches
#' @param out.dir Directory to write the diciontary in
#' @param alphabet Nucleotide alphabet. By defualt using all single/double nucleotide combinations and N
#' Note that the first four alphabets have to be \code{c('A', 'C', 'G', 'T')}
#' @param consecutive.mismatches Logical indicating whether only consecutive mismatches should be allowed
#' @param no.cores Number of cores of parallel processing
#' @param verbose Logical to show logs
#' @description This build a matrix of distances between all kmers and unique kmers,
#' considering reverse complements and matches to degenerate nucleotide alphabet.
#' @return Creates a file containing the following elements: \code{pairwise.kmers} Distance matrix
#' and \code{kmer.mapping} Mapping of reverse complements and kmers to features
#' The filename is constructed based on kmer length and kernel used and has the structure:
#' positional_wildcard_dict_kmer<\code{kmer.len}_mismatches\code{mismatches}_alpha\code{length (alphabet)}.Rdata
#' @export

build.wildcard.dictionary <- function (kmer.len, mismatches, out.dir,
		alphabet = c('A', 'C', 'G', 'T', 'N', 'S', 'M', 'R', 'Y', 'K', 'W'),
		consecutive.mismatches=TRUE, no.cores=1, verbose=TRUE) {

	start.time <- get.time ()
	
	# Some book keeping constants
	max <- length (alphabet) ^ kmer.len
	batch.size <- 5000000
	seq <- c(seq (-1, max-1, batch.size), max-1)

	# Build kmers
	time.start <- get.time ()
	res.list <- mclapply (2:length (seq), mc.cores=no.cores, mc.preschedule=FALSE,
		FUN=function (x) { determine.kmers (alphabet, mismatches, kmer.len, 
			(seq[x-1]+1):seq[x], consecutive.mismatches, verbose)})

	# Combine results from different cores
	kmers <- kmers.mismatches <- DNAStringSet ()
	for (i in 1:length (res.list)) {
		if (is.null (res.list[[i]]))
		next

		kmers <- c(kmers, res.list[[i]]$kmers)
		kmers.mismatches <- c(kmers.mismatches, res.list[[i]]$kmers.mismatches)
	}
	time.end <- get.time ()
	if (verbose)
   	    show (sprintf ("Time for finding kmers: %.2f", (time.end - time.start)/60))

	# Remove redundant kmers 
	kmers.mismatches <- remove.reverse.complements (kmers.mismatches)

	# Hamming distances
	time.start <- get.time ()
	if (verbose)
	    show ('Determining hamming distances...')
	dist.fwd <- hammingWithMismatches (kmers, kmers.mismatches, 0, FALSE, no.cores)
	if (verbose)
    	show ('Determining hamming distances for reverse complements...')
	dist.rev <- hammingWithMismatches (kmers, reverseComplement (kmers.mismatches), 0, FALSE, no.cores)
	time.end <- get.time ()
	if (verbose)
    	show (sprintf ("Time for determining distances: %.2f", (time.end - time.start) / 60))

	# Build features matrix
	if (verbose)
    	show ('Creating dictionary and cleaning up...')
	pairwise.kmers <- sparseMatrix (c(dist.fwd[,1], dist.rev[,1]),
		c(dist.fwd[,2], dist.rev[,2]), x=1, 
		dims=c(length (kmers), length (kmers.mismatches)),
		dimnames=list (as.character (kmers), as.character (kmers.mismatches)))
	gc ()

	# Eliminate repeated values
	summ <- summary (pairwise.kmers)
	pairwise.kmers <- sparseMatrix (summ[,1], summ[,2], x=1, 
		dims=c(nrow (pairwise.kmers), ncol (pairwise.kmers)),
		dimnames=list (as.character (kmers), as.character (kmers.mismatches)))
	gc ()

	# Column mapping
	kmer.column.mapping <- c(colnames (pairwise.kmers), colnames (pairwise.kmers))
	names (kmer.column.mapping) <- c(colnames (pairwise.kmers), 
		as.character (reverseComplement (DNAStringSet (colnames (pairwise.kmers)))))

	# File name
	file.name <- sprintf ("%s/wildcard_dict_kmer%d_mismatches%d_alpha%d", 
			out.dir, kmer.len, mismatches, length (alphabet))
	if (consecutive.mismatches)
	    file.name <- sprintf ("%s_consecutive_mis", file.name)
    file.name <- sprintf ("%s.Rdata", file.name)


	# Save results
	save (pairwise.kmers, kmer.column.mapping,
		file=file.name)
	end.time <- get.time ()
	show (sprintf ("Total time for build dictionary: %.2f", (end.time - start.time)/60))
}




#' Determine all kmers with mismatches of a given length
#'
#' @param alphabet Nucleotide alphabet
#' @param mismatches Number of mismatches allowed
#' @param range Range of numbers to determine kmers
#' @param kmer.len Length of kmers
#' @param verbose Logical to show logs
#' @description Builds all possible kmers in the given \code{range} for the given \code{alphabet}
#' @return A list of DNAStringSet elements: \code{kmers} containining kmers without mismatches
#' and \code{kmers.mismatches} containing all the kmers that were generated'
#' @export

determine.kmers <- function (alphabet, mismatches, kmer.len, 
	range=NULL, consecutive.mismatches, verbose=TRUE) {
	time.start <- get.time ()

	# Set up range if necessary 
	if (is.null (range)) {
		max <- length (alphabet) ^ kmer.len
		range <- 0:(max-1)
	}

	if (verbose)
	     show (sprintf ("Determining kmers for range: %d-%d", range[1], range[length (range)]))

	# Determine matrix of numbers
	mat <- digitsBase (range, length (alphabet), kmer.len)

	# Count mismatches
	temp <- mat; temp[temp < 4] <- 0;  temp[temp > 0] <- 1;
	seq.mismatches <- colSums (temp)

	# Reduce set to contain only allowed mismatches
	include.inds <- which (seq.mismatches <= mismatches)
	if (length (include.inds) == 0 )
	    return (NULL)

	temp <- temp[,include.inds]
	mat <- mat[,include.inds]
	seq.mismatches <- seq.mismatches[include.inds]	

	# Filter for consecutive mismatches
	if (consecutive.mismatches) {
		test.inds <- which (seq.mismatches > 1)
		temp.with.mismatches <- rbind (temp, seq.mismatches)[, test.inds]
		include.subset <- apply (temp.with.mismatches, 2, 
			function (x) { mis.ind <- kmer.len + 1;
				cs <- cumsum (x[1:kmer.len]); length (which (cs > 0 & cs < x[mis.ind])) == (x[mis.ind] - 1)})

		# Eliminate non consecutive mismatchs
		include.inds <- c(which (seq.mismatches <= 1), test.inds[include.subset])
		mat <- mat[, include.inds]
		seq.mismatches <- seq.mismatches[include.inds]
	}

	# Extract seqs 
	seqs <- DNAStringSet (apply (mat, 2, function (x) {
		paste (alphabet[x + 1], collapse = "")
		}))
	seqs <- matrix (alphabet[mat + 1], nrow=kmer.len)
	seqs <- DNAStringSet (apply (seqs, 2, paste, collapse=""))

	# Extract kmers with correct mismatches
	kmers.mismatches <- seqs
	# Store kmers without mismatches for row references
	kmers <- seqs[seq.mismatches == 0]

	time.end <- get.time ()
	if (verbose)
	    show (sprintf ("Time for %d-%d: %.2f", range[1], range[length (range)], (time.end - time.start)/60))

	return (list (kmers=kmers, kmers.mismatches=kmers.mismatches))
}


# Function to remove reverse complements
remove.reverse.complements <- function (kmers) {
	kmers.rc <- reverseComplement (kmers)
	inds <- kmers > kmers.rc
	temp <- kmers
	temp[inds] <- kmers.rc[inds]
	return (temp[!duplicated (temp)])
}
