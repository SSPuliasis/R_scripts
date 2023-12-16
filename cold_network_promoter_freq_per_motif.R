library(data.table)

fimo_results <- top_cluster_42_results

promoter_df <- data.frame(promoter_list)

unique_motif_ids <- na.omit(unique(fimo_results$motif_alt_id)) 
for (motif_name in unique_motif_ids){
  filtered_by_motif <- subset(fimo_results, motif_alt_id == motif_name)
  freqlist <- list()
  for (promoter in promoter_df$promoter_list){
    freqlist[[promoter]] <- nrow(subset(filtered_by_motif, sequence_name == promoter))
  }
    freqlist <- data.frame(freqlist)
    freqlist <- transpose(freqlist)
    colnames(freqlist) <- motif_name
    promoter_df <- cbind(promoter_df, motif_name = freqlist)
}

#change output file name here
write.csv(promoter_df, "cl42 promoter motif frequencies.csv")
