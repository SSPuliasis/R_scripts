count_seqs <- function(motif_list){
  for (motif in motif_list){
    motif_subset <- subset(top_results, motif_alt_id == motif)#subset the results for just the rows containing a certain motif id
    unique_genes <- unique(motif_subset$sequence_name)#count how many sequences that motif appears in
    print(motif)
    print(length(unique_genes))
  }
}

count_seqs(motifs)
