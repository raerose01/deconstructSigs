\name{deconstructSigs-package}
\alias{deconstructSigs-package}
\alias{deconstructSigs}
\docType{package}
\title{
 Identifies pre-defined signatures present in a tumor sample when given mutation context frequencies
}
\description{
  Takes sample information in the form of the fraction of mutations
  in each of 96 trinucleotide contexts and identifies the weighted combination 
  of published signatures that, when summed, most closely reconstructs the 
  mutational profile.
}
\details{
\tabular{ll}{
Package: \tab deconstructSigs\cr
Type: \tab Package\cr
Version: \tab 1.5\cr
Date: \tab 2015-12-08\cr
License: \tab GPL (>= 2)\cr
}
}
\author{

Rachel Rosenthal rachel.rosenthal.14@ucl.ac.uk
}
\examples{
test = whichSignatures(tumor.ref = randomly.generated.tumors,
                       signatures.ref = signatures,
                       sample.id = "2",
                       contexts.needed = FALSE,
                       trimer.counts.loc = tri.counts.exome)                     
plotSignatures(example.output, sub = "example")
makePie(example.output, sub = "example")
}
