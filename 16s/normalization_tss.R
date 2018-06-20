normalization_filter_tss<-function(biom_file, norm_method){
  if (norm_method == 'TSS'){
  return(transform_sample_counts(biom_file, function(x){x / sum(x)}))
  }
}