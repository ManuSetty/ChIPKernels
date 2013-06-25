#' @useDynLib ChIPKernels
#' @import Biostrings
#' @import methods
#' @import BiocGenerics

get.time <- function () 
as.numeric(format(Sys.time(), "%s"))
