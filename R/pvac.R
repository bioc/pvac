

pvacFilter <- function(abatch,pct=0.99){
  if(class(abatch) != "AffyBatch")
    stop("Wrong input, the first parameter should be an AffyBatch object")

  if(length(sampleNames(abatch))<6)
    warning("Sample size is small (<6), PCA-filtering may NOT be effective!")

  cat("Making absent/present calls, preprocessing data ... \n")
  calls  = mas5calls(abatch,verbose=FALSE)
  abatch = bg.correct.rma(abatch)
  abatch = normalize.AffyBatch.quantiles(abatch,type="pmonly")
  pm.mat  = log2(pm(abatch))
  pn.list = probeNames(abatch)
  pn.list = split(1:(length(pn.list)), pn.list)
  pca.func <- function(gid){
    xx = tryCatch(
      {
        e = pm.mat[pn.list[[gid]],]
        # o <- abs(e) == Inf | is.na(e)
        # e[o] <- min(e[!o])
        d = (svd(scale(t(e)))$d)^2
        d[1]/sum(d)
        # pca = prcomp(t(e), retx=FALSE, center=TRUE, scale =TRUE)
        # summary(pca)$importance[3,1]
      },
        error=function(e){
          0
        })
    xx
  }
  gn = featureNames(abatch)

  cat("Computing the PVAC scores ... \n")
  if(suppressWarnings(require(pbapply,quietly=TRUE))){
    # show progress bar
    pvac = pbsapply(gn,pca.func)
  }else{
    pvac = sapply(gn,pca.func)
    names(pvac) = gn
  }

  cat("Deriving the filtering threshold ... \n")

  gs = apply(exprs(calls),1,function(x) all(x=="A"))
  nullset = names(gs)[gs]
  cutoff  = min(quantile(pvac[nullset], pct),0.5)
  aset = names(pvac[pvac>cutoff])
  cat(paste("PVAC cutoff score (<=0.5): ",round(cutoff,5),"\n\n"))

  return(list(aset=aset,nullset=nullset,pvac=pvac,cutoff=cutoff))
}
