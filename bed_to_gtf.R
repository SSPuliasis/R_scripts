library(rtracklayer)
bed <- import.bed('C:/Users/SP43416/Documents/scratch/filtered_UP000005640_9606_proteome.bed')

### add artifical trans id
bed <- sort(bed,by = ~ name)
idx <- rle(bed$name)
idx <- sequence(idx$lengths)
trans <- paste0(bed$name,'.',idx)
bed$name <- paste0(bed$name,';',trans)

### get gtf
gtf <- unlist(blocks(bed))
gtf$gene_id <- gsub(";.*", "", names(gtf))
gtf$transcript_id <- gsub(".*;", "", names(gtf))
gtf$type <- 'exon'
names(gtf) <- NULL

## get cds
gr <- split(gtf,gtf$transcript_id)

thick <-GRanges(seqnames(bed),bed$thick)
thick$gene_id <- gsub(";.*", "", bed$name)
thick$transcript_id <- gsub(".*;", "", bed$name)
names(thick) <- thick$transcript_id
thick <- thick[names(gr)]

## get
cds <- pintersect(x = gr,y = thick,drop.nohit.ranges=T)
cds <- unlist(cds)
cds$type <- 'CDS'
names(cds) <- NULL

gtf$phase <- '.'
cds$phase <- '0'

output <- c(gtf,cds)
output <- sort(output,by = ~ seqnames + start + end)

export_gtf(output,file2save  = 'C:/Users/SP43416/Documents/scratch/filtered_UP000005640_9606_proteome.gtf')


#' Export Granges object to gtf file
#' @param gr Grange object 

export_gtf <- function(gr,file2save,
                       source.col='source',
                       feature.col=ifelse('type' %in% colnames(mcols(gr)),'type','feature'),
                       score='.',
                       phase='.',
                       mandatory = c("gene_id", "transcript_id"),
                       return.obj=F,
                       write.obj=T,
                       do_sort=T){
  options(scipen=100)
  if(NROW(gr) == 0){
    write.table(NULL, file=file2save, sep="\t", quote=F, row.names=F, col.names=F)
    message('The gtf is empty')
  } else {
    if(do_sort)
      gr <- sort(gr, by = ~ seqnames + start + end)
    ###check factor
    check <- mcols(gr)
    idx <- sapply(check@listData, is.factor)
    idx <- names(idx)[idx]
    for(i in idx){
      mcols(gr)[,i] <- as.character(mcols(gr)[,i])
    }
    
    idx <- unique(c('gene_id','transcript_id',feature.col,source.col,mandatory))
    if(!(source.col %in% colnames(mcols(gr))))
      mcols(gr)[,source.col] <- 'source'
    mcols(gr) <- mcols(gr)[,idx]
    ## change strands
    strd <- as.character(strand(gr))
    strd[strd=='*'] <- '.'
    score[is.na(score)] <- '.'
    phase[is.na(phase)] <- '.'
    
    gtf <- data.frame(chr=as.character(seqnames(gr)),
                      source=as.character(mcols(gr)[,source.col]),
                      feature=as.character(mcols(gr)[,feature.col]),
                      start=as.numeric(start(gr)),
                      end=as.numeric(end(gr)),
                      score=score,
                      strand=strd,
                      phase=phase,
                      attribute=makeGtfAttributes(df = as.data.frame(mcols(gr))),
                      stringsAsFactors = F)
    if(write.obj)
      write.table(gtf, file=file2save, sep="\t", quote=F, row.names=F, col.names=F)
    if(return.obj)
      return(gtf)
  }
  
  
}

makeGtfAttributes <- function(df, cols=NULL) {
  if (is.null(cols))
    cols = colnames(df)
  # make sure that gene_id and transcript_id are the first two columns
  mandatory = c("gene_id", "transcript_id")
  o = match(c(mandatory, setdiff(cols, mandatory)), cols)
  if (any(is.na(o[1:length(mandatory)]))) {
    o = o[!is.na(o)]
  }
  cols = cols[o]
  return(paste(apply(sapply(cols, function(s) {
    content = df[,s]
    if (is.character(content) | is.factor(content)) {
      content = paste('"', content, '"', sep="")
    }
    paste(gsub(".", "_", s, fixed=T), content, sep=" ")
  }), 1, paste, collapse="; "), ";", sep=""))
}
