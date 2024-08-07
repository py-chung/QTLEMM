% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MIM.search.R
\name{MIM.search}
\alias{MIM.search}
\title{QTL search by MIM}
\usage{
MIM.search(
  QTL,
  marker,
  geno,
  y,
  method = "EM",
  type = "RI",
  D.matrix = NULL,
  ng = 2,
  cM = TRUE,
  speed = 1,
  QTLdist = 15,
  conv = 10^-3,
  console = TRUE
)
}
\arguments{
\item{QTL}{matrix. A q*2 matrix contains the known QTL information,
where the row dimension q is the number of QTLs in the chromosomes. The
first column labels the chromosomes where the QTLs are located, and
the second column labels the positions of QTLs (in morgan (M) or
centimorgan (cM)). Note that chromosome and position must be divided
in order.}

\item{marker}{matrix. A k*2 matrix contains the marker information,
where the row dimension k is the number of markers in the chromosomes.
The first column labels the chromosomes where the markers are located,
and the second column labels the positions of QTLs (in morgan (M) or
centimorgan (cM)). Note that chromosome and position must be divided
in order.}

\item{geno}{matrix. A n*k matrix contains the k markers of the n
individuals. The marker genotypes of P1 homozygote (MM),
heterozygote (Mm) and P2 homozygote (mm) are coded as 2, 1 and 0,
respectively, and NA for missing value.}

\item{y}{vector. An vector with n elements that contains the phenotype
values of individuals.}

\item{method}{character. method="EM" means the interval mapping method
by Lander and Botstein (1989) is used in the analysis, while
method="REG" means  the approximate regression interval mapping method
by Haley and Knott (1992) is considered in the analysis.}

\item{type}{character. The population type of the dataset. Include
backcross (type="BC"), advanced intercross population (type="AI"), and
recombinant inbred population (type="RI").}

\item{D.matrix}{matrix. The design matrix of QTL effects which is a
g*p matrix, where g is the number of possible QTL genotypes, and p
is the number of effects considered in the MIM model. Note that the QTL
number of the design matrix must be the original QTL number plus one.
The design matrix can be easily generated by the function D.make(). If
being NULL, it will automatically generate a design matrix with all
additive and dominant effect and without any epistasis effect.}

\item{ng}{integer. The generation number of the population type. For
example, the BC1 population is type="BC" with ng=1; the AI F3
population is type="AI" with ng=3.}

\item{cM}{logical. Specify the unit of marker position. cM=TRUE for
centi-Morgan. Or cM=FALSE for Morgan.}

\item{speed}{numeric. The walking speed of the QTL search (in cM).}

\item{QTLdist}{numeric. The minimum distance (cM) among different
linked significant QTL. The position near the position of the known
QTLs under this distance will not be consider as the candidate position
in the search process.}

\item{conv}{numeric. The convergence criterion of EM algorithm.
The E and M steps will be iterated until a convergence criterion
is satisfied.}

\item{console}{logical. To decide whether the process of algorithm will
be shown in the R console or not.}
}
\value{
\item{effect}{The estimated effects, log likelihood value, and LRT
statistics of all searched positions.}
\item{QTL.best}{The positions of the best QTL combination.}
\item{effect.best}{The estimated effects and LRT statistics of the best
QTL combination.}
}
\description{
Expectation-maximization algorithm for QTL multiple interval mapping.
Find one more QTL in the presence of some known QTLs.
}
\examples{
# load the example data
load(system.file("extdata", "exampledata.RDATA", package = "QTLEMM"))

# run and result
QTL <- c(1, 23)
result <- MIM.search(QTL, marker, geno, y, type = "RI", ng = 2, speed = 15, QTLdist = 50)
result$QTL.best
result$effect.best
}
\references{
KAO, C.-H. and Z.-B. ZENG 1997 General formulas for obtaining the maximum
likelihood estimates and the asymptotic variance-covariance matrix in QTL
mapping when using the EM algorithm. Biometrics 53, 653-665.

KAO, C.-H., Z.-B. ZENG and R. D. TEASDALE 1999 Multiple interval mapping
for Quantitative Trait Loci. Genetics 152: 1203-1216.
}
\seealso{
\code{\link[QTLEMM]{EM.MIM}}
\code{\link[QTLEMM]{MIM.search2}}
}
