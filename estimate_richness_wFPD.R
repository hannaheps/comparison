##The following code is for the estimate_richness() function from phyloseq (McMurdie & Holmes 2013)
##Modified to add a step for calculating Faith's phylogenetic diversity as "FaithPD"

estimate_richness_wPD <- function(physeq, split=TRUE, measures=NULL){
  
  if( !any(otu_table(physeq)==1) ){
    # Check for singletons, and then warning if they are missing.
    # These metrics only really meaningful if singletons are included.
    warning(
      "The data you have provided does not have\n",
      "any singletons. This is highly suspicious. Results of richness\n",
      "estimates (for example) are probably unreliable, or wrong, if you have already\n",
      "trimmed low-abundance taxa from the data.\n",
      "\n",
      "We recommended that you find the un-trimmed data and retry."
    )
  }
  
  # If we are not splitting sample-wise, sum the species. Else, enforce orientation.
  if( !split ){
    OTU <- taxa_sums(physeq)		
  } else if( split ){
    OTU <- as(otu_table(physeq), "matrix")
    if( taxa_are_rows(physeq) ){ OTU <- t(OTU) }
  }
  
  # Define renaming vector:
  renamevec = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher", "FaithPD")
  names(renamevec) <- c("S.obs", "S.chao1", "S.ACE", "shannon", "simpson", "invsimpson", "fisher")
  # If measures was not explicitly provided (is NULL), set to all supported methods
  if( is.null(measures) ){
    measures = as.character(renamevec)
  }
  # Rename measures if they are in the old-style
  if( any(measures %in% names(renamevec)) ){
    measures[measures %in% names(renamevec)] <- renamevec[names(renamevec) %in% measures]
  }
  
  # Stop with error if no measures are supported
  if( !any(measures %in% renamevec) ){
    stop("None of the `measures` you provided are supported. Try default `NULL` instead.")
  }
  
  # Initialize to NULL
  outlist = vector("list")
  # Some standard diversity indices
  estimRmeas = c("Chao1", "Observed", "ACE")
  if( any(estimRmeas %in% measures) ){ 
    outlist <- c(outlist, list(t(data.frame(estimateR(OTU)))))
  }
  if( "Shannon" %in% measures ){
    outlist <- c(outlist, list(shannon = diversity(OTU, index="shannon")))
  }
  if( "Simpson" %in% measures ){
    outlist <- c(outlist, list(simpson = diversity(OTU, index="simpson")))
  }
  if( "InvSimpson" %in% measures ){
    outlist <- c(outlist, list(invsimpson = diversity(OTU, index="invsimpson")))
  }
  if( "FaithPD" %in% measures){
    outlist <- c(outlist, list(FaithPD = t(picante::pd(samp = OTU, tree = phy_tree(physeq), include.root = F))[1,] ))
  }
  if( "Fisher" %in% measures ){
    fisher = tryCatch(fisher.alpha(OTU, se=TRUE), 
                      warning=function(w){
                        warning("phyloseq::estimate_richness: Warning in fisher.alpha(). See `?fisher.fit` or ?`fisher.alpha`. Treat fisher results with caution")
                        suppressWarnings(fisher.alpha(OTU, se=TRUE)[, c("alpha", "se")])
                      }
    )
    if(!is.null(dim(fisher))){
      colnames(fisher)[1:2] <- c("Fisher", "se.fisher")
      outlist <- c(outlist, list(fisher))
    } else {
      outlist <- c(outlist, Fisher=list(fisher))
    }
  }
  out = do.call("cbind", outlist)
  # Rename columns per renamevec
  namechange = intersect(colnames(out), names(renamevec))
  colnames(out)[colnames(out) %in% namechange] <- renamevec[namechange]
  # Final prune to just those columns related to "measures". Use grep.
  colkeep = sapply(paste0("(se\\.){0,}", measures), grep, colnames(out), ignore.case=TRUE)
  out = out[, sort(unique(unlist(colkeep))), drop=FALSE]
  # Make sure that you return a data.frame for reliable performance.
  out <- as.data.frame(out)
  return(out)
}
