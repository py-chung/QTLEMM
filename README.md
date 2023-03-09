# QTLEMM

In this package, we provide the tool for QTL mapping by EM (Expectation-maximization) algorithm. It can handle various populations from different breeding scheme, including backcross (BC), F2, recombinant inbred (RI) populations, advanced intercrossed (AI) populations. For each population, both the complete genotyping data and selective genotyping data are considered. The functions allow to use linear regression, interval mapping ([Lander and Botstein 1989](https://academic.oup.com/genetics/article/121/1/185/5997927)), and multiple interval mapping ([Kao and Zeng 1997](https://biostat.wisc.edu/~kbroman/teaching/statgen/2006/refs/kao_zeng.pdf); [Kao, Zeng, and Teasdale 1999](https://academic.oup.com/genetics/article/152/3/1203/6094249?login=true); [Zeng, Kao, and Basten 1999](https://www.cambridge.org/core/journals/genetics-research/article/estimating-the-genetic-architecture-of-quantitative-traits/D5C17B27152E0240558490E02355D417)) methods for the detection of QTL.
  
We also provide the tool for QTL hotspot detecting in this package. We ([Yang et al. 2019](https://academic.oup.com/g3journal/article/9/2/439/6026674?login=true); [Wu et al. 2021](https://academic.oup.com/g3journal/article/11/4/jkab056/6151767)) develop a statistical framework that can carry out the QTL hotspot detecting process by LOD scores. And the permutation algorithm is deployed to randomly shift the tightly linked and/or pleiotropic QTL together along the genome separately by trait group for detecting spurious hotspots in the detection analysis.
  
## Installation
  
QTLEMM can be installed form GitHub by the following command:  
```install_github
# library(devtools)  
install_github("py-chung/QTLEMM", dependencies = TRUE, force = TRUE)
```
  
And QTLEMM can be installed form CRAN by the following command:
```install.packages
install.packages("QTLEMM")
```
  
## Main functions
  
+ `progeny()` Generate the simulated phenotype and genotype data.
+ `D.make()` Generate the genetic design matrix. 
+ `Q.make()` Generate the conditional probability matrix. 
+ `LRTthre()` LRT threshold for QTL interval mapping.
+ `IM.search()` QTL interval mapping to search the possible position of QTL in all chromosome. 
+ `IM.search2()` QTL interval mapping to search the possible position of QTL in all chromosome with selective genotyping. 
+ `EM.MIM()` EM algorithm for QTL multiple interval mapping.
+ `EM.MIM2()` EM algorithm for QTL multiple interval mapping with selective genotyping.
+ `MIM.search()` EM algorithm for QTL multiple interval mapping to find one more QTL by known QTLs.
+ `MIM.search2()` EM algorithm for QTL multiple interval mapping to find one more QTL by known QTLs with selective genotyping.
+ `MIM.points()` EM algorithm for QTL multiple interval mapping to find the best QTL position near the designated QTL position.
+ `MIM.points2()` EM algorithm for QTL multiple interval mapping to find the best QTL position near the designated QTL position with selective genotyping.
+ `LOD.QTLdetect()` Detect QTL by likelihood of odds(LOD) matrix.
+ `EQF.permu()` The EQF matrix cluster permutation process for QTL hotspot detection.
+ `EQF.plot()` Depict the EQF plot by the result of permutation process to detect the QTL hotspot.
  
More information can be seen in the following file:  
[Package ‘QTLEMM’](https://cran.r-project.org/web/packages/QTLEMM/QTLEMM.pdf)
  
## Example dataset
  
The example dataset is provided for performing hotspot detection by cluster permutation process. The LOD matrix of yeast data ([Brem and
Kruglyak 2005](https://www.pnas.org/doi/full/10.1073/pnas.0408709102)) is divided into four parts. The data can be downloaded from GitHub and combined into a complete matrix for use by the following command:
  
```yeast.LOD
load(url("https://github.com/py-chung/QTLEMM/raw/main/inst/extdata/yeast.LOD.1.RDATA"))
load(url("https://github.com/py-chung/QTLEMM/raw/main/inst/extdata/yeast.LOD.2.RDATA"))
load(url("https://github.com/py-chung/QTLEMM/raw/main/inst/extdata/yeast.LOD.3.RDATA"))
load(url("https://github.com/py-chung/QTLEMM/raw/main/inst/extdata/yeast.LOD.4.RDATA"))
yeast.LOD <- rbind(yeast.LOD.1,yeast.LOD.2,yeast.LOD.3,yeast.LOD.4)
```
  
And the imformation of bins of the LOD matrix can be downloaded from GitHub by the following command:
  
```yeast.LOD.bin
load(url("https://github.com/py-chung/QTLEMM/raw/main/inst/extdata/yeast.LOD.bin.RDATA"))
```







