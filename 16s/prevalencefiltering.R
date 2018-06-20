prevalence_filter<-function(biom_file,prev_th){

  print("
#Prevalence  filtering
#phyloseq provides useful tools for  ltering, subsetting, and agglomerating taxa – a task that is often appropriate or even necessary for effective analysis of microbiome count data. In this subsection, we graphically explore the prevalence of taxa in the example dataset, and demonstrate how this can be used as a  filtering criteria. One of the reasons to  filter in this way is to avoid spending much time analyzing taxa that were only rarely seen. This also turns out to be a useful  filter of noise (taxa that are actually just artifacts of the data collection process), a step that should probably be considered essential for datasets constructed via heuristic OTU-clustering methods, which are notoriously prone to generating spurious taxa.
# Define prevalence of each taxa
#(in how many samples did each taxa appear at least once))"
    )

  prev0 = apply(X = otu_table(biom_file),MARGIN = ifelse(taxa_are_rows(biom_file), yes = 1, no = 2),FUN = function(x){sum(x > 0)})
  prevdf = data.frame(Prevalence = prev0,TotalAbundance = taxa_sums(biom_file),tax_table(biom_file))
  keepPhyla <<-table(prevdf$Phylum)[(table(prevdf$Phylum) > prev_th)]
  print(names(keepPhyla))
  prevdf1 = subset(prevdf, Phylum %in% names(keepPhyla))
  # Define prevalence threshold as 5% of total samples 
  prevalenceThreshold = prev_th * 0.01 * nsamples(biom_file) 
  prevalenceThreshold
  # Execute prevalence  lter, using `prune_taxa()` function
  merged_mapping_biom_flt= prune_taxa((prev0 > prevalenceThreshold), biom_file) 
  merged_mapping_biom_flt
  
  # Filter entries with unidenti ed Phylum.
  merged_mapping_biom_flt_unidentified_Phylum = subset_taxa(merged_mapping_biom_flt, Phylum %in% names(keepPhyla))

  
  
  outlist<-list("flt_b"=merged_mapping_biom_flt_unidentified_Phylum, "prevdf1"=prevdf1, "prev_th"=prevalenceThreshold)
 
  return(outlist)
}
  