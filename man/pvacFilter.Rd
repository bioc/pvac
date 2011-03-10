\name{pvacFilter}
\alias{pvacFilter}
\title{
  Filter genes by the proportion of variation accounted for by the first
  principal component (PVAC)
}
\description{
  Compute the PVAC scores, derive a filtering threshold value, and
  return the names of probesets that have passed the filter
}
\usage{
  pvacFilter(abatch, pct=0.99)
}
\arguments{
  \item{abatch}{an instance of
    \code{\link[affy:AffyBatch-class]{AffyBatch}} from the function call \code{\link[affy]{ReadAffy}}
  }
  \item{pct}{the percentile value of the emperical distribution of PVAC
    scores of a set of ``non-expressed'' genes. Used to select the
    filtering threshold. The default value is 0.99.
  }

}
\value{
A list with the following components,

  \item{aset }{Names of the probesets that have passed the filter}
  \item{nullset }{Names of the presumably ``non-expressed'' probesets
    (those with absent calls across all the study samples)}
  \item{pvac }{A named vector containing the PVAC scores of all probesets}
  \item{cutoff }{The PVAC cutoff value. The maximum is set to 0.5 (which
  corresponds to 50\% of the total variation in a probeset)}

}
\details{
  This function implements a new filtering method for Affymetrix GeneChips, based
  on principal component analysis (PCA) on the probe-level expression
  data. Given that all the probes in a probeset are designed to target one or a common
  cluster of transcripts, the measurements of probes in a probeset
  should be correlated. The degree of concordance of gene expression
  among probes can be approximated by the proportion of variation
  accounted by the first principal component (PVAC). Using
  a wholly defined spike-in dataset, we have shown that
  filtering by PVAC provides increased sensitivity in detecting
  truly differentially expressed genes while controlling the false
  discoveries. The filtering threshold value is chosen from the PVAC
  score distribution in a set of ``non-expressed'' gene (those with absent calls in all samples).
}

\author{
Jun Lu
}
\examples{
    if ( require(affydata) ) {
    data(Dilution)
    res = pvacFilter(Dilution)
    res$aset[1:5]  # 5 probesets that have passed the filter
 }
}
\keyword{ filter }

