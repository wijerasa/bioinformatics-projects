normalization_filter_css<-function(obj,present, depeth, norm_method){
  if (norm_method == 'CSS'){
    obj<-filterData(obj, present=present, depth=depth)
    p<-cumNormStatFast(obj)
    obj<-cumNorm(obj, p = p)
  }
  return(obj)
}
