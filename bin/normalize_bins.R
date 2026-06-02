library(data.table)

# read in MULTIBAMSUMMARY output, sample names (sample sheet), and output file name (interally generated)
args     <- commandArgs(trailingOnly=TRUE)
bed      <- args[1]
samples  <- strsplit(args[2], ",")[[1]]
outfile  <- args[3]

# Function for BPM normalization from the chromDicts R package:
NormalizeBins <- function(input_file, output_file=NULL, sample_names=NULL){
  res <- utils::read.delim(input_file, comment.char = "#", header=F)
  
  bs <- as.numeric(res[[2]][2]-res[[2]][1])
  
  res[4:ncol(res)] <- lapply(res[4:ncol(res)], as.numeric)
  # normalize:
  res[4:ncol(res)] <-
    apply(res[4:ncol(res)], 2, function(x){
      return((as.numeric(x)/sum(as.numeric(x)))*1e06)
    })
  
  # Add column names to res:
  if(is.null(sample_names)){
    colnames(res) <- c("chr", "start", "end", colnames(res)[4:ncol(res)])
  } else {
    colnames(res) <- c("chr", "start", "end", sample_names)
  }
  
  # add a single position for each bin at the centerpoint:
  res$pos <- res$start + as.integer(bs*0.5)
  
  # remove start and end from res:
  res <- res[c(1,ncol(res),4:(ncol(res)-1))]
  
  if(is.null(output_file)){
    return(res)
  } else {
    message("saving output file...")
    utils::write.table(x = res, file = output_file, sep="\t", quote=F, col.names = T, row.names = F)
    return(res)
  }
}

NormalizeBins(
  input_file   = bed,
  output_file  = outfile,
  sample_names = samples
)