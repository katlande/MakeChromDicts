library(data.table)

# read in sample names (sample sheet), sample conditions (sample sheet), and BPM table (NORMALIZEBINS output)
args       <- commandArgs(trailingOnly=TRUE)
samples    <- strsplit(args[1], ",")[[1]]
conditions <- strsplit(args[2], ",")[[1]]
bpm_file   <- args[3]

# make meta data data.frame:
df  <- data.frame(sample=samples, condition=conditions, stringsAsFactors=FALSE)
# read in BPM table (computed in workflow.2)
bpm <- read.delim(bpm_file) 

# function to make a chromdict, from the chromDicts R package
chromDictAny <- function(mat, IDs=NULL, multiple.samples=T, remove.extra=T){
  
  # set first two column names to "chr" and "pos"
  colnames(mat)[1:2] <- c("chr", "pos")
  
  if(multiple.samples){
    message("Calculating mean value differences...")
    if(length(IDs) == 1){
      mat$value <- rowMeans(mat[which(grepl(IDs, colnames(mat)))])
    } else {
      mat$value <- rowMeans(mat[which(colnames(mat) %in% IDs)])
    }
    
  } else {
    message(paste("Extracting values for sample:", IDs))
    if(length(IDs) > 1){
      warning("Cannot have multiple IDs when multiple.samples=FALSE!")
      stop()
    }
    # VAL column is just the column of interest
    mat$value <- mat[[which(colnames(mat) == IDs)]]
  }
  
  # create a master data.frame of bp chrom, pos, and value:
  mat$cpgID <- paste0(mat$chr, ":", mat$pos)
  row.names(mat) <- mat$cpgID
  diffdf <- mat[, c("chr", "pos", "value")]
  
  # if specified, remove any contig or chromosome other than autosomes and sex chromosomes. Note that this will NOT WORK in genome versions that use abnormal chromosome nomenclature, and is designed for human and mouse:
  if(remove.extra){
    allchrnames <- unique(diffdf$chr)
    allchrnames_strip <- gsub("chr", "", allchrnames, ignore.case = TRUE)
    keep.chrnames <- allchrnames[which(allchrnames_strip %in% c(1:1000, "x", "y", "X", "Y"))]
    # only subset if there are actually chromosomes present to remove:
    if(! identical(allchrnames, keep.chrnames)){
      message(paste("Removing", (length(allchrnames)-length(keep.chrnames)), "abnormal chromosomes/contigs..."))
      diffdf <- subset(diffdf, chr %in% keep.chrnames)
    } else {
      message("No abnormal chromosomes/contigs found!")
    }
  }
  
  # split by chromosome and order+index with data.table
  outlist <- list()
  message("Splitting by chromosome...")
  for(i in unique(as.character(diffdf$chr))){
    message(i)
    d <- diffdf[diffdf$chr==i,]
    d$pos <- as.numeric(d$pos)
    d <- d[order(d$pos),]
    data.table::setDT(d)
    data.table::setkey(d,pos)
    outlist <- append(outlist, list(d))
  }
  
  names(outlist) <- unique(as.character(diffdf$chr))
  return(outlist)
}

# apply chromDictAny() to all conditions:
lapply(unique(df$condition), function(cond) {
  tmpdict <- chromDictAny(bpm, IDs=df$sample[df$condition==cond],
                          multiple.samples=TRUE, remove.extra=TRUE)
  saveRDS(tmpdict, paste0(cond, ".rds"))
  rm(tmpdict); gc() # clear memory
})