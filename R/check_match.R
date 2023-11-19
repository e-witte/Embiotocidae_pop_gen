check_match <- function(locus, thresh, table) {
  vals <- ((locus-thresh):(locus+thresh))
  x<- any(vals %in% table)
  return(x)
}
