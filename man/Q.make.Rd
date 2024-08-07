% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Q.make.R
\name{Q.make}
\alias{Q.make}
\title{Generate Q Matrix}
\usage{
Q.make(
  QTL,
  marker,
  geno = NULL,
  interval = FALSE,
  type = "RI",
  ng = 2,
  cM = TRUE
)
}
\arguments{
\item{QTL}{matrix. A q*2 matrix contains the QTL information, where
the row dimension q is the number of QTLs in the chromosomes. The
first column labels the chromosomes where the QTLs are located, and
the second column labels the positions of QTLs (in morgan (M) or
centimorgan (cM)). Note that chromosome and position must be divided
in order.}

\item{marker}{matrix. A k*2 matrix contains the marker information,
where the row dimension k is the number of markers in the chromosomes.
The first column labels the chromosomes where the markers are located,
and the second column labels the positions of QTLs (in morgan (M) or
centimorgan (cM)). Note that chromosomes and positions must be divided
in order.}

\item{geno}{matrix. A n*k matrix contains the k markers of the n
individuals. The marker genotypes of P1 homozygote (MM),
heterozygote (Mm) and P2 homozygote (mm) are coded as 2, 1 and 0,
respectively, and NA for missing value.}

\item{interval}{logical. When the QTL with the same position of a
marker, whether the marker will be skipped and not be regarded as a
flanking marker or not. interval=TRUE presents the marker will be
skipped.}

\item{type}{character. The population type of the dataset. Include
backcross (type="BC"), advanced intercross population (type="AI"), and
recombinant inbred population (type="RI").}

\item{ng}{integer. The generation number of the population type. For
example, the BC1 population is type="BC" with ng=1; the AI F3
population is type="AI" with ng=3.}

\item{cM}{logical. Specify the unit of marker position. cM=TRUE for
centi-Morgan. Or cM=FALSE for Morgan.}
}
\value{
The output contains k conditional probability matrices for the k
flanking marker pairs (the k Q-matrices) and a conditional
probability matrix of each QTL for all individuals (the cp-matrix,
if the genotype data of testing population is input).
}
\description{
Generate the conditional probability matrix by the information of QTL
and marker and the genotype data.
}
\note{
If geno=NULL, the function can be run too and the output will
contain k Q-matrices but no cp-matrix.
}
\examples{
# load the example data
load(system.file("extdata", "exampledata.RDATA", package = "QTLEMM"))

# run and result
result <- Q.make(QTL, marker, geno)
head(result$cp.matrix)
}
\references{
KAO, C.-H. and Z.-B. ZENG 1997 General formulas for obtaining the maximum
likelihood estimates and the asymptotic variance-covariance matrix in QTL
mapping when using the EM algorithm. Biometrics 53, 653-665.

KAO, C.-H., Z.-B. ZENG and R. D. TEASDALE 1999 Multiple interval mapping
for Quantitative Trait Loci. Genetics 152: 1203-1216.
}
