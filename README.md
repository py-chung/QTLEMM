# QTLEMM

In this package, we provide the tool for QTL mapping by EM algorithm. It can handle various populations from different breeding scheme, including backcross (BC), F2, recombinant inbred (RI) populations, advanced intercrossed (AI) populations populations. For each population, both the complete genotyping data and selective genotyping data are considered. The functions allow to use linear regression, interval mapping ([Lander and Botstein 1989](https://academic.oup.com/genetics/article/121/1/185/5997927)), and multiple interval mapping ([Kao and Zeng 1997](https://biostat.wisc.edu/~kbroman/teaching/statgen/2006/refs/kao_zeng.pdf); [Kao, Zeng, and Teasdale 1999](https://academic.oup.com/genetics/article/152/3/1203/6094249?login=true); [Zeng, Kao, and Basten 1999](https://www.cambridge.org/core/journals/genetics-research/article/estimating-the-genetic-architecture-of-quantitative-traits/D5C17B27152E0240558490E02355D417)) methods for the detection of QTL.
  
We also provide the tool for QTL hotspot detecting in this package. We ([Yang et al. 2019](https://academic.oup.com/g3journal/article/9/2/439/6026674?login=true); [Wu et al. 2021](https://academic.oup.com/g3journal/article/11/4/jkab056/6151767)) develop a statistical framework that can carry out the QTL hotspot detecting process by LOD scores. And the permutation algorithm is deployed to randomly shift the tightly linked and/or pleiotropic QTL together along the genome separately by trait group for detecting the spurious hotspots in the detection analysis.
