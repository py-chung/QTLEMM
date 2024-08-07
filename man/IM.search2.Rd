% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IM.search2.R
\name{IM.search2}
\alias{IM.search2}
\title{QTL search by IM with Selective Genotyping}
\usage{
IM.search2(
  marker,
  geno,
  y,
  yu = NULL,
  sele.g = "n",
  tL = NULL,
  tR = NULL,
  method = "EM",
  type = "RI",
  D.matrix = NULL,
  ng = 2,
  cM = TRUE,
  speed = 1,
  conv = 10^-5,
  d.eff = FALSE,
  LRT.thre = TRUE,
  simu = 1000,
  alpha = 0.05,
  detect = TRUE,
  QTLdist = 15,
  plot.all = TRUE,
  plot.chr = TRUE,
  console = TRUE
)
}
\arguments{
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

\item{y}{vector. A vector that contains the phenotype values of
individuals with genotyped.}

\item{yu}{vector. A vector that contains the phenotype value
of the individuals without genotyped.}

\item{sele.g}{character. If sele.g="n", it will consider that the
data is not a selective genotyping data but the complete genotyping
data. If sele.g="p", it will consider that the data is a selective
genotyping data, and use the proposed model (Lee 2014) to analyze.
If sele.g="t", it will consider that the data is a selective
genotyping data, and use the truncated model (Lee 2014) to analyze.
If sele.g="f, it will consider that the data is a selective
genotyping data, and use the frequency-based model (Lee 2014) to
analyze. Note that the yu must be input when sele.g="p" of "f".}

\item{tL}{numeric. The lower truncation point of phenotype value
when sele.g="t". Note that when sele.g="t" and tL=NULL, the yu
must be input and the function will consider the minimum of yu
as the lower truncation point.}

\item{tR}{numeric. The upper truncation point of phenotype value
when sele.g="t". Note that when sele.g="t" and tR=NULL, the yu
must be input and the function will consider the maximum of yu
as the upper truncation point.}

\item{method}{character. method="EM" means the interval mapping method
by Lander and Botstein (1989) is used in the analysis, while
method="REG" means  the approximate regression interval mapping method
by Haley and Knott (1992) is considered in the analysis.}

\item{type}{character. The population type of the dataset. Include
backcross (type="BC"), advanced intercross population (type="AI"), and
recombinant inbred population (type="RI").}

\item{D.matrix}{matrix. The design matrix of the IM model. If
D.matrix=NULL, the design matrix will be the constructed using the
Cockerham’s model. In BC population, it is a 2*1 matrix which contains
0.5, -0.5 for additive effect. In RI or AI population, it is a 3*2 matrix
which contains  1, 0, -1 for additive effect and 0.5, -0.5, 0.5 for
dominant effect.}

\item{ng}{integer. The generation number of the population type. For
example, the BC1 population is type="BC" with ng=1; the AI F3
population is type="AI" with ng=3.}

\item{cM}{logical. Specify the unit of marker position. cM=TRUE for
centi-Morgan. Or cM=FALSE for Morgan.}

\item{speed}{numeric. The walking speed of the QTL search (in cM).}

\item{conv}{numeric. The convergent criterion of EM algorithm.
The E and M steps will be iterated until a convergent criterion
is satisfied.}

\item{d.eff}{logical. Specify whether the dominant effect will be
considered in the parameter estimation or not for AI or RI population.}

\item{LRT.thre}{logical or numeric. If being TRUE, the LRT threshold
will be computed based on the Gaussian stochastic process
(Kao and Ho 2012). Or users can input a numerical value as the LRT
threshold to assessing the significance of QTL detection.}

\item{simu}{integer. To decide how many simulation samples will be used
to compute the LRT threshold using the Gaussian process.}

\item{alpha}{numeric. The type I error rate for the LRT threshold.}

\item{detect}{logical. Whether the significant QTL whose LRT statistic
is larger than the LRT threshold will be shown in the output dataset or
not.}

\item{QTLdist}{numeric. The minimum distance (cM) among different
linked significant QTL.}

\item{plot.all}{logical. If being TRUE, output the profile of LRT
statistics for the genome in one figure.}

\item{plot.chr}{logical. If being TRUE, output the profile of LRT
statistics for the chromosomes.}

\item{console}{logical. To decide whether the process of algorithm will
be shown in the R console or not.}
}
\value{
\item{effect}{The estimated effects and LRT statistics of all positions.}
\item{LRT.threshold}{The LRT threshold value computed for the data using the
Gaussian stochastic process (Kuo 2011; Kao and Ho 2012).}
\item{detect.QTL}{The positions, effects and LRT statistics of the detected
QTL significant using the obtained LRT threshold value.}
\item{model}{The model of selective genotyping data in this analyze.}

Graphical outputs including LOD value and effect of each position.
}
\description{
Expectation-maximization algorithm for QTL interval mapping to search
the possible position of QTL in all chromosome. This function can
handle the genotype data witch is selective genotyping.
}
\examples{
# load the example data
load(system.file("extdata", "exampledata.RDATA", package = "QTLEMM"))

# make the seletive genotyping data
ys <- y[y > quantile(y)[4] | y < quantile(y)[2]]
yu <- y[y >= quantile(y)[2] & y <= quantile(y)[4]]
geno.s <- geno[y > quantile(y)[4] | y < quantile(y)[2],]

# run and result
result <- IM.search2(marker, geno.s, ys, yu, sele.g = "p", type = "RI", ng = 2,
speed = 7.5, conv = 10^-3, LRT.thre = 10)
result$detect.QTL
}
\references{
KAO, C.-H. and Z.-B. ZENG 1997 General formulas for obtaining the maximum
likelihood estimates and the asymptotic variance-covariance matrix in QTL
mapping when using the EM algorithm. Biometrics 53, 653-665.

KAO, C.-H., Z.-B. ZENG and R. D. TEASDALE 1999 Multiple interval mapping
for Quantitative Trait Loci. Genetics 152: 1203-1216.

H.-I LEE, H.-A. HO and C.-H. KAO 2014 A new simple method for improving
QTL mapping under selective genotyping. Genetics 198: 1685-1698.

KAO, C.-H. and H.-A. Ho 2012 A score-statistic approach for determining
threshold values in QTL mapping. Frontiers in Bioscience. E4, 2670-2682.
}
\seealso{
\code{\link[QTLEMM]{EM.MIM2}}
\code{\link[QTLEMM]{IM.search}}
\code{\link[QTLEMM]{LRTthre}}
}
