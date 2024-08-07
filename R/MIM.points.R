#' QTL Short Distance Correction by MIM
#'
#' Expectation-maximization algorithm for QTL multiple interval mapping.
#' Find the best QTL position near the designated QTL position.
#'
#' @param QTL matrix. A q*2 matrix contains the QTL information, where
#' the row dimension q is the number of QTLs in the chromosomes. The
#' first column labels the chromosomes where the QTLs are located, and
#' the second column labels the positions of QTLs (in morgan (M) or
#' centimorgan (cM)). Note that chromosome and position must be divided
#' in order.
#' @param marker matrix. A k*2 matrix contains the marker information,
#' where the row dimension k is the number of markers in the chromosomes.
#' The first column labels the chromosomes where the markers are located,
#' and the second column labels the positions of QTLs (in morgan (M) or
#' centimorgan (cM)). Note that chromosome and position must be divided
#' in order.
#' @param geno matrix. A n*k matrix contains the k markers of the n
#' individuals. The marker genotypes of P1 homozygote (MM),
#' heterozygote (Mm) and P2 homozygote (mm) are coded as 2, 1 and 0,
#' respectively, and NA for missing value.
#' @param y vector. An vector with n elements that contains the phenotype
#' values of individuals.
#' @param method character. method="EM" means the interval mapping method
#' by Lander and Botstein (1989) is used in the analysis, while
#' method="REG" means  the approximate regression interval mapping method
#' by Haley and Knott (1992) is considered in the analysis.
#' @param type character. The population type of the dataset. Include
#' backcross (type="BC"), advanced intercross population (type="AI"), and
#' recombinant inbred population (type="RI").
#' @param D.matrix matrix. The design matrix of QTL effects which is a
#' g*p matrix, where g is the number of possible QTL genotypes, and p
#' is the number of effects considered in the MIM model. The design
#' matrix can be easily generated by the function D.make(). If being NULL,
#' it Will automatically generate a design matrix with all additive and
#' dominant effect and without any epistasis effect.
#' @param ng integer. The generation number of the population type. For
#' example, the BC1 population is type="BC" with ng=1; the AI F3
#' population is type="AI" with ng=3.
#' @param cM logical. Specify the unit of marker position. cM=TRUE for
#' centi-Morgan. Or cM=FALSE for Morgan.
#' @param scope numeric vector. The search scope of each QTL. In the
#' MIM process, it will search forward and backward for the corresponding
#' cM. User can assign a numeric number for every QTL or a numeric vector
#' for each QTL. Note that 0 denote that the corresponding QTL position
#' is fixed, and the position of its surrounding positions will not be
#' searched.
#' @param speed numeric. The walking speed of the QTL search (in cM).
#' @param conv numeric. The convergence criterion of EM algorithm.
#' The E and M steps will be iterated until a convergence criterion
#' is satisfied.
#' @param console logical. To decide whether the process of algorithm will
#' be shown in the R console or not.
#'
#' @return
#' \item{effect}{The estimated effects, log likelihood value, and LRT
#' statistics of all searched positions.}
#' \item{QTL.best}{The positions of the best QTL combination.}
#' \item{effect.best}{The estimated effects and LRT statistics of the best
#' QTL combination.}
#'
#' @export
#'
#' @references
#'
#' KAO, C.-H. and Z.-B. ZENG 1997 General formulas for obtaining the maximum
#' likelihood estimates and the asymptotic variance-covariance matrix in QTL
#' mapping when using the EM algorithm. Biometrics 53, 653-665.
#'
#' KAO, C.-H., Z.-B. ZENG and R. D. TEASDALE 1999 Multiple interval mapping
#' for Quantitative Trait Loci. Genetics 152: 1203-1216.
#'
#' @seealso
#' \code{\link[QTLEMM]{EM.MIM}}
#' \code{\link[QTLEMM]{MIM.points2}}
#'
#' @examples
#' # load the example data
#' load(system.file("extdata", "exampledata.RDATA", package = "QTLEMM"))
#'
#' # run and result
#' result <- MIM.points(QTL, marker, geno, y, type = "RI", ng = 2, scope = c(0,3,0), speed = 2)
#' result$QTL.best
#' result$effect.best
MIM.points <- function(QTL, marker, geno, y, method = "EM", type = "RI", D.matrix = NULL, ng = 2, cM = TRUE,
                      scope = 5, speed = 1, conv = 10^-3, console = TRUE){

  if(is.null(QTL) | is.null(marker) | is.null(geno) |  is.null(y)){
    stop("Input data is missing, please cheak and fix.", call. = FALSE)
  }

  genotest <- table(geno)
  datatry <- try(geno*geno, silent=TRUE)
  if(class(datatry)[1] == "try-error" | FALSE %in% (names(genotest) %in% c(0, 1, 2))  | !is.matrix(geno)){
    stop("Genotype data error, please cheak your genotype data.", call. = FALSE)
  }

  marker <- as.matrix(marker)
  markertest <- c(ncol(marker) != 2, NA %in% marker, marker[,1] != sort(marker[,1]), nrow(marker) != ncol(geno))
  datatry <- try(marker*marker, silent=TRUE)
  if(class(datatry)[1] == "try-error" | TRUE %in% markertest){
    stop("Marker data error, or the number of marker does not match the genetype data.", call. = FALSE)
  }

  QTL <- as.matrix(QTL)
  datatry <- try(QTL*QTL, silent=TRUE)
  if(class(datatry)[1] == "try-error" | ncol(QTL) != 2 | NA %in% QTL | max(QTL[, 1]) > max(marker[, 1])){
    stop("QTL data error, please cheak your QTL data.", call. = FALSE)
  }

  for(i in 1:nrow(QTL)){
    ch0 <- QTL[i, 1]
    q0 <- QTL[i, 2]
    if(!ch0 %in% marker[, 1]){
      stop("The specified QTL is not in the position that can estimate by the marker data.", call. = FALSE)
    }
    if(q0 > max(marker[marker[, 1] == ch0, 2]) | q0 < min(marker[marker[, 1] == ch0, 2])){
      stop("The specified QTL is not in the position that can estimate by the marker data.", call. = FALSE)
    }
  }

  y[is.na(y)] <- mean(y,na.rm = TRUE)

  datatry <- try(y%*%geno, silent=TRUE)
  if(class(datatry)[1] == "try-error"){
    stop("Phenotype data error, or the number of individual does not match the genetype data.", call. = FALSE)
  }

  nq <- nrow(QTL)
  ch0 <- QTL[, 1]

  if(!type[1] %in% c("AI","RI","BC") | length(type) > 1){
    stop("Parameter type error, please input AI, RI, or BC.", call. = FALSE)
  }

  if(!is.null(D.matrix)){
    D.matrix <- as.matrix(D.matrix)
    datatry <- try(D.matrix*D.matrix, silent=TRUE)
    if(type == "BC"){dn0 <- 2
    } else {dn0 <- 3}
    if(class(datatry)[1] == "try-error" | NA %in% D.matrix | nrow(D.matrix) != dn0^nq){
      stop("Parameter D.matrix error, or the combination of genotypes in design matrix is error.",
           call. = FALSE)
    }
  }

  if(!is.numeric(ng) | length(ng) > 1 | min(ng) < 1){
    stop("Parameter ng error, please input a positive integer.", call. = FALSE)
  }
  ng <- round(ng)

  if(!cM[1] %in% c(0,1) | length(cM > 1)){cM <- TRUE}

  if(!is.numeric(speed) | length(speed) > 1 | min(speed) < 0){
    stop("Parameter speed error, please input a positive number.", call. = FALSE)
  }

  scope <- c(scope)
  datatry <- try(scope*scope, silent=TRUE)
  if(class(datatry)[1] == "try-error" | NA %in% scope | length(scope) > nq | min(scope) < 0){
    stop("Parameter scope error, please input the positive integers or 0.", call. = FALSE)
  }

  if(length(scope) < nq){
    scope <- rep(scope[1],nq)
  }

  if(!is.numeric(conv) | length(conv) > 1 | min(conv) < 0){
    stop("Parameter conv error, please input a positive number.", call. = FALSE)
  }

  if(!console[1] %in% c(0,1) | length(console) > 1){console <- TRUE}

  if(!cM){
    QTL <- QTL*100
    marker[, 2] <- marker[, 2]*100
  }

  sr <- list()
  q0 <- QTL[1, 2]
  r0 <- union(seq(q0-scope[1], q0, speed), seq(q0, q0+scope[1], speed))
  r0 <- r0[r0 >= min(marker[marker[, 1] == ch0[1], 2]) & r0 <= max(marker[marker[, 1] == ch0[1], 2])]
  sr[[1]] <- r0
  if(nq > 1){
    for(i in 2:nq){
      q0 <- QTL[i, 2]
      r0 <- union(seq(q0-scope[i], q0, speed), seq(q0, q0+scope[i], speed))
      r0 <- r0[r0 >= min(marker[marker[, 1] == ch0[1], 2]) & r0 <= max(marker[marker[, 1] == ch0[1], 2])]
      if(ch0[i] == ch0[i-1]){
        r0 <- r0[r0 > max(sr[[i-1]])]
      }
      sr[[i]] <- r0
    }
  }

  sl <- unlist(lapply(sr,length))
  sm <- matrix(NA, prod(sl), length(sl))
  sl2 <- c(1,sl,1)
  for(i in 1:length(sl)){
    sm[, i] <- rep(sr[[i]], each = prod(sl2[(i+2):length(sl2)]), time = prod(sl2[1:i]))
  }

  if(is.null(D.matrix)){
    D.matrix <- D.make(as.numeric(nq), type = type)
  }

  name0 <- cbind(paste("QTL", 1:nq, ".ch", sep = ""), paste("QTL", 1:nq, ".cM", sep = ""))
  if(console){cat("#", t(name0), "LRT", "log.likelihood", "\n", sep = "\t")}

  if(method == "EM"){
    meth <- function(D.matrix, cp.matrix, y, conv){
      EM <- EM.MIM(D.matrix, cp.matrix, y, conv = conv, console = FALSE)
      eff <- as.numeric(EM$E.vector)
      mu0 <- as.numeric(EM$beta)
      sigma <- sqrt(as.numeric(EM$variance))
      LRT <- EM$LRT
      like <- EM$log.likelihood
      R2 <- EM$R2
      result <- list(eff, mu0, sigma, LRT, like, R2)
      return(result)
    }
  } else if (method == "REG"){
    meth <- function(D.matrix, cp.matrix, y, conv){
      X <- cp.matrix%*%D.matrix
      fit <- stats::lm(y~X)
      eff <- as.numeric(fit$coefficients[-1])
      mu0 <- as.numeric(fit$coefficients[1])
      ms <- stats::anova(fit)$`Mean Sq`
      sigma <- ms[2]^0.5
      R2 <- summary(fit)$r.squared

      L0 <- c()
      L1 <- c()
      for(k in 1:nrow(cp0)){
        L00 <- c()
        L01 <- c()
        for(m in 1:nrow(D.matrix)){
          L00[m] <- cp0[k,m]*stats::dnorm((y[k]-mu0)/sigma)
          L01[m] <- cp0[k,m]*stats::dnorm((y[k]-(mu0+D.matrix[m,]%*%eff))/sigma)
        }
        L0[k] <- sum(L00)
        L1[k] <- sum(L01)
      }
      like0 <- sum(log(L0))
      like1 <- sum(log(L1))
      LRT <- 2*(like1-like0)

      result <- list(eff, mu0, sigma, LRT, like1, R2)
      return(result)
    }
  }

  effect <- matrix(NA, nrow(sm), 2*nq+ncol(D.matrix)+3)
  for(i in 1:nrow(sm)){
    QTL0 <- cbind(ch0, sm[i,])
    cp0 <- Q.make(QTL0, marker, geno, type = type, ng = ng)$cp.matrix
    fit0 <- meth(D.matrix, cp0, y, conv)
    effect0 <- c(t(QTL0), fit0[[1]], fit0[[4]], fit0[[5]], fit0[[6]])

    LRT <- round(effect0[length(effect0)-2], 3)
    like <- round(effect0[length(effect0)-1], 3)
    if(console){cat(paste(i, "/", nrow(sm), sep = ""), t(QTL0), LRT, like, "\n", sep = "\t")}
    effect[i,] <- effect0
  }
  colnames(effect) <-  c(t(name0), colnames(D.matrix), "LRT", "log.likelihood", "R2")
  row.names(effect) <- 1:nrow(effect)

  best <- effect[effect[,ncol(effect)-1] == max(effect[,ncol(effect)-1]),]
  QTL.best <- matrix(best[1:(2*nq)], nq, 2, byrow = TRUE)
  colnames(QTL.best) <- c("chromosome", "position(cM)")

  effect.best <- best[-(1:(2*nq))]

  result <- list(effect = effect, QTL.best = QTL.best, effect.best = effect.best)
  return(result)
}
