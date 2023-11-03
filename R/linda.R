winsor.fun <- function(Y, quan) {
  N <- colSums(Y)
  P <- t(t(Y) / N)
  cut <- apply(P, 1, quantile, quan)
  Cut <- matrix(rep(cut, ncol(Y)), nrow(Y))
  ind <- P > Cut
  P[ind] <- Cut[ind]
  Y <- round(t(t(P) * N))
  return(Y)
}

#' Linear (Lin) model for differential abundance (DA) analysis
#'
#' The function implements a simple, robust and highly scalable approach to tackle
#' the compositional effects in differential abundance analysis. It fits linear regression models
#' on the centered log2-ratio transformed data, identifies a bias term due to the transformation
#' and compositional effect, and corrects the bias using the mode of the regression coefficients.
#' It could fit mixed-effect models.
#'
#' @param otu.tab data frame or matrix representing observed OTU table. Row: taxa; column: samples.
#' NAs are not expected in OTU tables so are not allowed in function \code{linda}.
#' @param meta data frame of covariates. The rows of \code{meta} correspond to the columns of \code{otu.tab}.
#' NAs are allowed. If there are NAs, the corresponding samples will be removed in the analysis.
#' @param formula character. For example: \code{formula = '~x1*x2+x3+(1|id)'}. At least one fixed effect is required.
#' @param type character; the type of \code{otu.tab}. It could be "count" or "proportion". The default is "count".
#' @param adaptive TRUE or FALSE. The default is TRUE. If TRUE, the parameter \code{imputation} will be treated as FALSE no matter
#' what it is actually set to be. Then the significant correlations between the sequencing depth and explanatory variables will be tested
#' via the linear regression between the log of the sequencing depths and \code{formula}. If any p-value is smaller than or equal to
#' \code{corr.cut}, the imputation approach will be used; otherwise, the pseudo-count approach will be used.
#' The information of whether the imputation or pseudo-count approach is used will be printed.
#' @param imputation TRUE or FALSE. The default is FALSE. If TRUE, then we use the imputation approach, i.e., zeros in \code{otu.tab} will be
#' imputed using the formula in the referenced paper.
#' @param pseudo.cnt a positive real value. The default is 0.5. If \code{adaptive} and \code{imputation} are both FALSE,
#' then we use the pseudo-count approach, i.e., we add \code{pseudo.cnt} to each value in \code{otu.tab}.
#' @param corr.cut a real value between 0 and 1; significance level of correlations between the sequencing depth and
#' explanatory variables. The default is 0.1.
#' @param p.adj.method character; p-value adjusting approach. See R function \code{p.adjust}. The default is 'BH'.
#' @param alpha a real value between 0 and 1; significance level of differential abundance. The default is 0.05.
#' @param prev.cut a real value between 0 and 1; taxa with prevalence (percentage of nonzeros)
#' less than prev.cut are excluded. The default is 0 (no taxa will be excluded).
#' @param lib.cut a non-negative real value; samples with less than \code{lib.cut} read counts are excluded.
#' The default is 1 (no samples will be excluded).
#' @param winsor.quan a real value between 0 and 1; winsorization cutoff (quantile) for \code{otu.tab}, e.g., 0.97. The default is NULL.
#' If NULL, winsorization process will not be conducted.
#' @param n.cores a positive integer. If \code{n.cores > 1} and formula is in a form of mixed-effect model,
#' \code{n.cores} parallels will be conducted. The default is 1.
#'
#' @return A list with the elements
#' \item{variables}{A vector of variable names of all fixed effects in \code{formula}. For example: \code{formula = '~x1*x2+x3+(1|id)'}.
#' Suppose \code{x1} and \code{x2} are numerical, and \code{x3} is a categorical variable of three levels: a, b and c.
#' Then the elements of \code{variables} would be \code{('x1', 'x2', 'x3b', 'x3c', 'x1:x2')}.}
#' \item{bias}{numeric vector; each element corresponds to one variable in \code{variables};
#' the estimated bias of the regression coefficients due to the compositional effect.}
#' \item{output}{a list of data frames with columns 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj', 'reject',
#'  'df'; \code{names(output)} is equal to \code{variables}; the rows of the data frame corresponds to taxa.
#'  Note: if there are taxa being excluded due to \code{prev.cut}, the number of the rows of the output data frame
#'  will be not equal to the number of the rows of \code{otu.tab}. Taxa are identified by the rownames.
#'  If the rownames of \code{otu.tab} are NULL, then \code{1 : nrow(otu.tab)} is set as the rownames of \code{otu.tab}.
#'  \itemize{
#'  \item{baseMean:}{ 2 to the power of the intercept coefficients (normalized by one million)}
#'  \item{log2FoldChange:}{ bias-corrected coefficients}
#'  \item{lfcSE:}{ standard errors of the coefficients}
#'  \item{stat:}{ \code{log2FoldChange / lfcSE}}
#'  \item{pvalue:}{ \code{2 * pt(-abs(stat), df)}}
#'  \item{padj:}{ \code{p.adjust(pvalue, method = p.adj.method)}}
#'  \item{reject:}{ \code{padj <= alpha}}
#'  \item{df:}{ degrees of freedom. The number of samples minus the number of explanatory variables (intercept included) for
#'  fixed-effect models; estimates from R package \code{lmerTest} with Satterthwaite method of approximation for mixed-effect models.}
#'  }}
#' \item{covariance}{a list of data frames; the data frame records the covariances between a certain regression coefficient and other coefficients;
#'  \code{names(covariance)} is equal to \code{variables}; the rows of the data frame corresponds to taxa. If the length of \code{variables}
#'  is equal to 1, then the \code{covariance} is NULL.}
#' \item{otu.tab.use}{the OTU table used in the abundance analysis (the \code{otu.tab} after the preprocessing:
#' samples that have NAs in the variables in \code{formula} or have less than \code{lib.cut} read counts are removed;
#' taxa with prevalence less than \code{prev.cut} are removed and data is winsorized if \code{!is.null(winsor.quan)};
#' and zeros are treated, i.e., imputed or pseudo-count added).}
#' \item{meta.use}{the meta data used in the abundance analysis (only variables in \code{formula} are stored; samples that have NAs
#' or have less than \code{lib.cut} read counts are removed; numerical variables are scaled).}
#' \item{wald}{a list for use in Wald test. If the fitting model is a linear model, then it includes
#' \itemize{
#' \item{beta:}{ a matrix of the biased regression coefficients including intercept and all fixed effects; the culumns correspond to taxa}
#' \item{sig:}{ the standard errors; the elements corresponding to taxa}
#' \item{X:}{ the design matrix of the fitting model}
#' \item{bias:}{ the estimated biases of the regression coefficients including intercept and all fixed effects}
#' }
#' If the fitting model is a linear mixed-effect model, then it includes
#' \itemize{
#' \item{beta:}{ a matrix of the biased regression coefficients including intercept and all fixed effects; the culumns correspond to taxa}
#' \item{beta.cov:}{ a list of covariance matrices of \code{beta}; the elements corresponding to taxa}
#' \item{rand.cov:}{ a list with covariance matrices of variance parameters of random effects; the elements corresponding to taxa; see more details in the paper of 'lmerTest'}
#' \item{Joc.beta.cov.rand:} { a list of a list of Jacobian matrices of \code{beta.cov} with respect to the variance parameters; the elements corresponding to taxa}
#' \item{bias:}{ the estimated biases of the regression coefficients including intercept and all fixed effects}
#' }}
#'
#' @author Huijuan Zhou \email{huijuanzhou2019@gmail.com}
#' Jun Chen \email{Chen.Jun2@mayo.edu}
#' Xianyang Zhang \email{zhangxiany@stat.tamu.edu}
#' @references Huijuan Zhou, Kejun He, Jun Chen, and Xianyang Zhang. LinDA: Linear Models for Differential Abundance
#' Analysis of Microbiome Compositional Data.
#' @importFrom modeest mlv
#' @importFrom lmerTest lmer
#' @importFrom lmerTest as_lmerModLmerTest
#' @import foreach
#' @import parallel
#' @examples
#'
#' #install package "phyloseq" for importing "smokers" dataset
#' ind <- smokers$meta$AIRWAYSITE == 'Throat'
#' otu.tab <- as.data.frame(smokers$otu[, ind])
#' meta <- cbind.data.frame(Smoke = factor(smokers$meta$SMOKER[ind]),
#'                          Sex = factor(smokers$meta$SEX[ind]),
#'                          Site = factor(smokers$meta$SIDEOFBODY[ind]),
#'                          SubjectID = factor(smokers$meta$HOST_SUBJECT_ID[ind]))
#' ind1 <- which(meta$Site == 'Left')
#' res.left <- linda(otu.tab[, ind1], meta[ind1, ], formula = '~Smoke+Sex', alpha = 0.1,
#'                   prev.cut = 0.1, lib.cut = 1000, winsor.quan = 0.97)
#' ind2 <- which(meta$Site == 'Right')
#' res.right <- linda(otu.tab[, ind2], meta[ind2, ], formula = '~Smoke+Sex', alpha = 0.1,
#'                   prev.cut = 0.1, lib.cut = 1000, winsor.quan = 0.97)
#' rownames(res.left$output[[1]])[which(res.left$output[[1]]$reject)]
#' rownames(res.right$output[[1]])[which(res.right$output[[1]]$reject)]
#'
#' linda.obj <- linda(otu.tab, meta, formula = '~Smoke+Sex+(1|SubjectID)', alpha = 0.1,
#'                    prev.cut = 0.1, lib.cut = 1000, winsor.quan = 0.97)
#' linda.plot(linda.obj, c('Smokey', 'Sexmale'),
#'            titles = c('Smoke: n v.s. y', 'Sex: female v.s. male'), alpha = 0.1, lfc.cut = 1,
#'            legend = TRUE, directory = NULL, width = 11, height = 8)
#'
#' @export

linda <- function(otu.tab, meta, formula, type = 'count',
                  adaptive = TRUE, imputation = FALSE, pseudo.cnt = 0.5, corr.cut = 0.1,
                  p.adj.method = 'BH', alpha = 0.05,
                  prev.cut = 0, lib.cut = 1, winsor.quan = NULL, n.cores = 1) {
  if(any(is.na(otu.tab))) {
    stop('The OTU table contains NAs! Please remove!\n')
  }
  allvars <- all.vars(as.formula(formula))
  Z <- as.data.frame(meta[, allvars])

  ## preprocessing
  keep.sam <- which(colSums(otu.tab) >= lib.cut & rowSums(is.na(Z)) == 0)
  Y <- otu.tab[, keep.sam]
  Z <- as.data.frame(Z[keep.sam, ])
  names(Z) <- allvars

  n <- ncol(Y)
  keep.tax <- which(rowSums(Y > 0) / n >= prev.cut)
  Y <- Y[keep.tax, ]
  m <- nrow(Y)

  ## some samples may have zero total counts after screening taxa
  if(any(colSums(Y) == 0)) {
    ind <- which(colSums(Y) > 0)
    Y <- Y[, ind]
    Z <- as.data.frame(Z[ind, ])
    names(Z) <- allvars
    keep.sam <- keep.sam[ind]
    n <- ncol(Y)
  }

  ## scaling numerical variables
  ind <- sapply(1 : ncol(Z), function(i) is.numeric(Z[, i]))
  Z[, ind] <- scale(Z[, ind])

  ## winsorization
  if(!is.null(winsor.quan)) {
    Y <- winsor.fun(Y, winsor.quan)
  }

  ##
  if(grepl('\\(', formula)) {
    random.effect <- TRUE
  } else {
    random.effect <- FALSE
  }

  if(is.null(rownames(otu.tab))) {
    taxa.name <- (1 : nrow(otu.tab))[keep.tax]
  } else {
    taxa.name <- rownames(otu.tab)[keep.tax]
  }
  if(is.null(rownames(meta))) {
    samp.name <- (1 : nrow(meta))[keep.sam]
  } else {
    samp.name <- rownames(meta)[keep.sam]
  }

  ## handling zeros
  if(type == 'count') {
    if(any(Y == 0)) {
      N <- colSums(Y)
      if(adaptive) {
        logN <- log(N)
        if(random.effect) {
          tmp <- lmer(as.formula(paste0('logN', formula)), Z)
        } else {
          tmp <- lm(as.formula(paste0('logN', formula)), Z)
        }
        corr.pval <- coef(summary(tmp))[-1, "Pr(>|t|)"]
        if(any(corr.pval <= corr.cut)) {
          cat('Imputation approach is used.\n')
          imputation <- TRUE
        } else {
          cat('Pseudo-count approach is used.\n')
          imputation <- FALSE
        }
      }
      if(imputation) {
        N.mat <- matrix(rep(N, m), nrow = m, byrow = TRUE)
        N.mat[Y > 0] <- 0
        tmp <- N[max.col(N.mat)]
        Y <- Y + N.mat / tmp
      } else {
        Y <- Y + pseudo.cnt
      }
    }
  }

  if(type == 'proportion') {
    if(any(Y == 0)) {
      ## Half minimum approach
      Y <- t(apply(Y, 1, function (x) {
        x[x == 0] <- 0.5 * min(x[x != 0])
        return(x)
      }))
    }
  }

  ## CLR transformation
  logY <- log2(Y)
  W <- t(logY) - colMeans(logY)

  ## linear regression
  oldw <- getOption('warn')
  options(warn = -1)
  if(!random.effect) {
    suppressMessages(fit <- lm(as.formula(paste0('W', formula)), Z))
    res <- do.call(rbind, coef(summary(fit)))
    d <- ncol(model.matrix(fit))
    df <- rep(n - d, m)
    tmp <- vcov(fit)
    res.cov <- foreach(i = 1 : m) %do% {tmp[((i-1)*d+1) : (i*d), ((i-1)*d+1) : (i*d)]}
    wald <- list(beta = coef(fit), sig = sigma(fit), X = model.matrix(fit))
    res.cov <- do.call(rbind, res.cov)
    rownames(res.cov) <- rownames(res)
    colnames(res.cov) <- rownames(res)[1 : d]
  } else {
    fun <- function(i) {
      w <- W[, i]
      fit <- lmer(as.formula(paste0('w', formula)), Z)
      a <- as_lmerModLmerTest(fit)
      rand.cov <- a@vcov_varpar
      Jac.beta.cov.rand <- a@Jac_list
      list(coef(summary(fit)), vcov(fit), rand.cov, Jac.beta.cov.rand)
    }
    if(n.cores > 1) {
      tmp <- mclapply(c(1 : m), function(i) fun(i), mc.cores = n.cores)
    } else {
      suppressMessages(tmp <- foreach(i = 1 : m) %do% fun(i))
    }
    res <- do.call(rbind, lapply(tmp, `[[`, 1))
    res.cov <- do.call(rbind, lapply(tmp, `[[`, 2))
    wald <- list(beta = do.call(cbind, lapply(lapply(tmp, `[[`, 1), function(x)x[,1])),
                 beta.cov = lapply(tmp, `[[`, 2),
                 rand.cov = lapply(tmp, `[[`, 3),
                 Jac.beta.cov.rand = lapply(tmp, `[[`, 4))
  }
  options(warn = oldw)

  res.intc <- res[which(rownames(res) == '(Intercept)'), ]
  rownames(res.intc) <- NULL
  options(warn = -1)
  suppressMessages(bias.intc <- mlv(sqrt(n) * res.intc[, 1],
                               method = 'meanshift', kernel = 'gaussian') / sqrt(n))
  options(warn = oldw)
  baseMean <- 2 ^ (res.intc[, 1] - bias.intc)
  baseMean <- baseMean / sum(baseMean) * 1e6

  output.fun <- function(x) {
    res.voi <- res[which(rownames(res) == x), ]
    rownames(res.voi) <- NULL

    if(random.effect) {
      df <- res.voi[, 3]
    }

    log2FoldChange <- res.voi[, 1]
    lfcSE <- res.voi[, 2]
    oldw <- getOption('warn')
    options(warn = -1)
    suppressMessages(bias <- mlv(sqrt(n) * log2FoldChange,
                                 method = 'meanshift', kernel = 'gaussian') / sqrt(n))
    options(warn = oldw)
    log2FoldChange <- log2FoldChange - bias
    stat <- log2FoldChange / lfcSE

    pvalue <- 2 * pt(-abs(stat), df)
    padj <- p.adjust(pvalue, method = p.adj.method)
    reject <- padj <= alpha
    output <- cbind.data.frame(baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, reject, df)
    rownames(output) <- taxa.name
    return(list(bias = bias, output = output))
  }

  cov.fun <- function(x) {
    tmp <- (1 : ncol(res.cov))[-c(1, which(colnames(res.cov) == x))]
    covariance <- as.data.frame(as.matrix(res.cov[which(rownames(res.cov) == x), tmp]))
    rownames(covariance) <- taxa.name
    colnames(covariance) <- colnames(res.cov)[tmp]
    return(covariance)
  }

  variables <- unique(rownames(res))[-1]
  variables.n <- length(variables)
  bias <- rep(NA, variables.n)
  output <- list()
  if(variables.n == 1) {
    covariance <- NULL
  } else {
    covariance <- list()
  }
  for(i in 1 : variables.n) {
    tmp <- output.fun(variables[i])
    output[[i]] <- tmp[[2]]
    bias[i] <- tmp[[1]]
    if(variables.n > 1) {
      covariance[[i]] <- cov.fun(variables[i])
    }
  }
  names(output) <- variables
  if(variables.n > 1) {
    names(covariance) <- variables
  }

  rownames(Y) <- taxa.name
  colnames(Y) <- samp.name
  rownames(Z) <- samp.name
  wald[['bias']] <- c(bias.intc, bias)
  return(list(variables = variables, bias = bias, output = output, covariance = covariance, otu.tab.use = Y, meta.use = Z, wald = wald))
}

#' Wald test for bias-corrected regression coefficient
#'
#' The function implements Wald test for bias-corrected regression coefficient learned from the \code{linda} function.
#' One can utilize the function to perform ANOVA-type analyses.
#'
#' @param linda.obj return from the \code{linda} function.
#' @param L A matrix for testing \code{Lb = 0}, where \code{b} includes the intercept and all fixed effects. Thus the number of columns of
#' L must be equal to \code{length(variables)+1}, where \code{variables} is from \code{linda.obj}, which does not include the intercept.
#' @param model \code{'LM'} or \code{'LMM'} indicating the model fitted in \{linda\} is linear model or linear mixed-effect model.
#' @param alpha significance level for testing \code{Lb = 0}.
#' @param p.adj.method P-value adjusting approach. See R function \code{p.adjust}. The default is 'BH'.
#'
#' @return A data frame with columns
#' \item{Fstat}{Wald statistics for each taxon}
#' \item{df1}{The numerator degrees of freedom}
#' \item{df2}{The denominator degrees of freedom}
#' \item{pvalue}{ \code{1 - pf(Fstat, df1, df2)}}
#' \item{padj}{ \code{p.adjust(pvalue, method = p.adj.method)}}
#' \item{reject}{ \code{padj <= alpha}}
#'
#' @author Huijuan Zhou \email{huijuanzhou2019@gmail.com}
#' Jun Chen \email{Chen.Jun2@mayo.edu}
#' Xianyang Zhang \email{zhangxiany@stat.tamu.edu}
#' @references Huijuan Zhou, Kejun He, Jun Chen, and Xianyang Zhang. LinDA: Linear Models for Differential Abundance
#' Analysis of Microbiome Compositional Data.
#' @examples
#'
#' #install package "phyloseq" for importing "smokers" dataset
#' ind <- smokers$meta$AIRWAYSITE == 'Throat'
#' otu.tab <- as.data.frame(smokers$otu[, ind])
#' meta <- cbind.data.frame(Smoke = factor(smokers$meta$SMOKER[ind]),
#'                          Sex = factor(smokers$meta$SEX[ind]),
#'                          Site = factor(smokers$meta$SIDEOFBODY[ind]),
#'                          SubjectID = factor(smokers$meta$HOST_SUBJECT_ID[ind]))
#' linda.obj <- linda(otu.tab, meta, formula = '~Smoke+Sex+(1|SubjectID)+(Smoke|Site)', alpha = 0.1,
#'                    prev.cut = 0.1, lib.cut = 1000, winsor.quan = 0.97)
#' L <- matrix(c(0, 1, 0, 0, 0, 1), nrow = 2, byrow = TRUE)
#' result <- linda.wald.test(linda.obj, L, 'LMM', alpha = 0.1)
#'
#' @export

linda.wald.test <- function(linda.obj, L, model = c('LM', 'LMM'), alpha = 0.05, p.adj.method = 'BH') {
  r <- qr(L)$rank
  wald <- linda.obj$wald
  bias <- wald$bias
  beta <- wald$beta - bias
  n <- ncol(linda.obj$otu.tab.use)
  m <- nrow(linda.obj$otu.tab.use)

  if(model == 'LM') {
    p <- nrow(beta)
    X <- wald$X
    Lbeta.cov.inv <- solve(L %*% solve(t(X)%*% X) %*% t(L))
    Lbeta <- L %*% beta
    Fstat <- sapply(1 : m, function(i) {
      lbeta <- Lbeta[, i]
      t(lbeta) %*% Lbeta.cov.inv %*% lbeta / (wald$sig[i] ^ 2 * r)
    })
    df2 <- rep(n - p, m)
  } else if(model == 'LMM') {
    beta.cov <- wald$beta.cov
    rand.cov <- wald$rand.cov
    Jac.beta.cov.rand <- wald$Jac.beta.cov.rand
    Lbeta <- L %*% beta
    fun <- function(i) {
      beta.cov.i <- beta.cov[[i]]
      rand.cov.i <- rand.cov[[i]]
      Jac.beta.cov.rand.i <- Jac.beta.cov.rand[[i]]
      Lbeta.cov.i <- L %*% beta.cov.i %*% t(L)
      if(nrow(L) == 1) {
        d <- as.vector(Lbeta.cov.i)
        g <- sapply(Jac.beta.cov.rand.i, function(x) L %*% x %*% t(L))
        nu <- 2 * d ^ 2 / as.vector(t(g) %*% rand.cov.i %*% g)
        Lbeta.cov.inv.i <- 1 / d
      } else {
        Lbeta.cov.i.eig <- eigen(Lbeta.cov.i)
        P <- Lbeta.cov.i.eig$vectors
        D <- Lbeta.cov.i.eig$values
        PtL <- t(P) %*% L
        nu <- rep(NA, r)
        for(j in 1 : r) {
          d <- D[j]
          l <- PtL[j, ]
          g <- sapply(Jac.beta.cov.rand.i, function(x) t(l) %*% x %*% l)
          nu[j] <- 2 * d ^ 2 / as.vector(t(g) %*% rand.cov.i %*% g)
        }
        Lbeta.cov.inv.i <- P %*% diag(D ^ (-1)) %*% t(P)
      }
      tmp <- sum(nu / (nu - 2))
      nu <- 2 * tmp / (tmp - r)
      lbeta <- Lbeta[, i]
      Fstat <- as.vector(t(lbeta) %*% Lbeta.cov.inv.i %*% lbeta) / r
      df2 <- nu
      return(c(Fstat, df2))
    }
    Fstat <- df2 <- rep(NA, m)
    for (i in 1 : m) {
      res <- fun(i)
      Fstat[i] <- res[1]
      df2[i] <- res[2]
    }
  }
  df1 <- rep(r, m)
  pvalue <- 1 - pf(Fstat, df1, df2)
  padj <- p.adjust(pvalue, method = p.adj.method)
  reject <- padj <= alpha
  res <- cbind.data.frame(Fstat, df1, df2, pvalue, padj, reject)
  rownames(res) <- rownames(linda.obj$otu.tab.use)
  return(res)
}

#' Plot linda results
#'
#' The function plots the effect size plot and volcano plot based on the output from \code{linda}.
#'
#' @param linda.obj return from function \code{linda}.
#' @param variables.plot vector; variables whose results are to be plotted. For example, suppose the return
#' value \code{variables} is equal to \code{('x1', 'x2', 'x3b', 'x3c', 'x1:x2')}, then one could set \code{variables.plot = c('x3b', 'x1:x2')}.
#' @param titles vector; titles of the effect size plot and volcano plot for each variable in \code{variables.plot}.
#' The default is NULL. If NULL, the titles will be set as \code{variables.plot}.
#' @param alpha a real value between 0 and 1; cutoff for \code{padj}.
#' @param lfc.cut a positive value; cutoff for \code{log2FoldChange}.
#' @param legend TRUE or FALSE; whether to show the legends of the effect size plot and volcano plot.
#' @param directory character; the directory to save the figures, e.g., \code{getwd()}. The default is NULL. If NULL, figures will not be saved.
#' @param width the width of the graphics region in inches. See R function \code{pdf}.
#' @param height the height of the graphics region in inches. See R function \code{pdf}.
#'
#' @return A list of \code{ggplot2} objects.
#' \item{plot.lfc}{a list of effect size plots. Each plot corresponds to one variable in \code{variables.plot}.}
#' \item{plot.volcano}{a list of volcano plots. Each plot corresponds to one variable in \code{variables.plot}.}
#'
#' @author Huijuan Zhou \email{huijuanzhou2019@gmail.com}
#' Jun Chen \email{Chen.Jun2@mayo.edu}
#' Xianyang Zhang \email{zhangxiany@stat.tamu.edu}
#' @references Huijuan Zhou, Kejun He, Jun Chen, and Xianyang Zhang. LinDA: Linear Models for Differential Abundance
#' Analysis of Microbiome Compositional Data.
#' @import ggplot2
#' @import ggrepel
#' @examples
#'
#' #install package "phyloseq" for importing "smokers" dataset
#' ind <- smokers$meta$AIRWAYSITE == 'Throat'
#' otu.tab <- as.data.frame(smokers$otu[, ind])
#' meta <- cbind.data.frame(Smoke = factor(smokers$meta$SMOKER[ind]),
#'                          Sex = factor(smokers$meta$SEX[ind]),
#'                          Site = factor(smokers$meta$SIDEOFBODY[ind]),
#'                          SubjectID = factor(smokers$meta$HOST_SUBJECT_ID[ind]))
#' ind1 <- which(meta$Site == 'Left')
#' res.left <- linda(otu.tab[, ind1], meta[ind1, ], formula = '~Smoke+Sex', alpha = 0.1,
#'                   prev.cut = 0.1, lib.cut = 1000, winsor.quan = 0.97)
#' ind2 <- which(meta$Site == 'Right')
#' res.right <- linda(otu.tab[, ind2], meta[ind2, ], formula = '~Smoke+Sex', alpha = 0.1,
#'                   prev.cut = 0.1, lib.cut = 1000, winsor.quan = 0.97)
#' rownames(res.left$output[[1]])[which(res.left$output[[1]]$reject)]
#' rownames(res.right$output[[1]])[which(res.right$output[[1]]$reject)]
#'
#' linda.obj <- linda(otu.tab, meta, formula = '~Smoke+Sex+(1|SubjectID)', alpha = 0.1,
#'                    prev.cut = 0.1, lib.cut = 1000, winsor.quan = 0.97)
#' linda.plot(linda.obj, c('Smokey', 'Sexmale'),
#'            titles = c('Smoke: n v.s. y', 'Sex: female v.s. male'), alpha = 0.1, lfc.cut = 1,
#'            legend = TRUE, directory = NULL, width = 11, height = 8)
#'
#' @export

linda.plot <- function(linda.obj, variables.plot, titles = NULL, alpha = 0.05, lfc.cut = 1,
                       legend = FALSE, directory = NULL, width = 11, height = 8) {
  bias <- linda.obj$bias
  output <- linda.obj$output
  otu.tab <- linda.obj$otu.tab.use
  meta <- linda.obj$meta.use
  variables <- linda.obj$variables
  if(is.null(titles)) titles <- variables.plot

  taxa <- rownames(otu.tab)
  m <- length(taxa)

  tmp <- match(variables, variables.plot)
  voi.ind <- order(tmp)[1 : sum(!is.na(tmp))]
  padj.mat <- foreach(i = voi.ind, .combine = 'cbind') %do% {
    output[[i]]$padj
  }

  ## effect size plot
  if(is.matrix(padj.mat)) {
    ind <- which(colSums(padj.mat <= alpha) > 0)
  } else if(is.vector(padj.mat)) {
    tmp <- which(padj.mat <= alpha)
    if(length(tmp) > 0){
      ind <- 1
    } else {
      ind <- integer(0)
    }
  }

  if(length(ind) == 0) {
    plot.lfc <- NULL
  } else {
    if(!is.null(directory)) pdf(paste0(directory, '/plot_lfc.pdf'),
                                width = width, height = height)
    plot.lfc <- list()
    j <- 1
    for(i in ind) {
      output.i <- output[[voi.ind[i]]]
      bias.i <- bias[voi.ind[i]]
      lfc <- output.i$log2FoldChange
      lfcSE <- output.i$lfcSE
      padj <- output.i$padj

      ind.rej <- which(padj <= alpha)
      n.rej <- length(ind.rej)
      taxa.rej <- taxa[ind.rej]
      taxa.rej <- factor(taxa.rej, levels = taxa.rej)
      data.plot.lfc <- cbind.data.frame(Taxa = rep(taxa.rej, 2),
                                        Log2FoldChange = c(lfc[ind.rej], lfc[ind.rej] + bias.i),
                                        lfcSE = c(lfcSE[ind.rej], rep(NA, n.rej)),
                                        bias = rep(c('Debiased', 'Non-debiased'), each = n.rej))
      plot.lfc.i <- ggplot(data.plot.lfc, aes(x = Log2FoldChange, y = Taxa)) +
        geom_point(aes(color = bias, shape = bias), size = 3) +
        geom_errorbar(aes(xmin = Log2FoldChange - 1.96 * lfcSE,
                          xmax = Log2FoldChange + 1.96 * lfcSE), width = .2) +
        geom_vline(xintercept = 0, color = 'gray', linetype = 'dashed') +
        ggtitle(titles[i]) +
        theme_bw(base_size = 18)
      if(legend) {
        plot.lfc.i <- plot.lfc.i +
          theme(legend.title = element_blank(),
                legend.key.width = unit(1, 'cm'), plot.margin = unit(c(1, 1, 1, 1.5), 'cm'))
      } else {
        plot.lfc.i <- plot.lfc.i +
          theme(legend.position = 'none', plot.margin = unit(c(1, 1, 1, 1.5), 'cm'))
      }
      plot.lfc[[j]] <- plot.lfc.i
      j <- j + 1
      if(!is.null(directory)) print(plot.lfc.i)
    }
    if(!is.null(directory)) dev.off()
  }

  ## volcano plot
  plot.volcano <- list()
  if(!is.null(directory)) pdf(paste0(directory, '/plot_volcano.pdf'),
                              width = width, height = height)
  leg1 <- paste0('padj>', alpha, ' & ', 'lfc<=', lfc.cut)
  leg2 <- paste0('padj>', alpha, ' & ', 'lfc>', lfc.cut)
  leg3 <- paste0('padj<=', alpha, ' & ', 'lfc<=', lfc.cut)
  leg4 <- paste0('padj<=', alpha, ' & ', 'lfc>', lfc.cut)

  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1 : n]
  }
  color <- gg_color_hue(3)

  for(i in 1 : length(voi.ind)) {
    output.i <- output[[voi.ind[i]]]
    bias.i <- bias[voi.ind[i]]
    lfc <- output.i$log2FoldChange
    padj <- output.i$padj

    ind1 <- padj > alpha & abs(lfc) <= lfc.cut
    ind2 <- padj > alpha & abs(lfc) > lfc.cut
    ind3 <- padj <= alpha & abs(lfc) <= lfc.cut
    ind4 <- padj <= alpha & abs(lfc) > lfc.cut

    leg <- rep(NA, m)
    leg[ind1] = leg1
    leg[ind2] = leg2
    leg[ind3] = leg3
    leg[ind4] = leg4
    leg <- factor(leg, levels = c(leg1, leg2, leg3, leg4))
    taxa.sig <- rep('', m)
    taxa.sig[ind3 | ind4] <- taxa[ind3 | ind4]

    data.volcano <- cbind.data.frame(taxa = taxa.sig, Log2FoldChange = lfc,
                                     Log10Padj = -log10(padj), leg = leg)
    plot.volcano.i <- ggplot(data.volcano, aes(x = Log2FoldChange, y = Log10Padj)) +
      geom_point(aes(color = leg), size = 2) +
      geom_text_repel(aes(label = taxa), max.overlaps = Inf) +
      scale_colour_manual(values = c('darkgray', color[c(2, 3, 1)])) +
      geom_hline(aes(yintercept = -log10(alpha)), color = 'gray', linetype = 'dashed') +
      geom_vline(aes(xintercept = -lfc.cut), color = 'gray', linetype = 'dashed') +
      geom_vline(aes(xintercept = lfc.cut), color = 'gray', linetype = 'dashed') +
      ylab('-Log10Padj') +
      ggtitle(titles[i]) +
      theme_bw(base_size = 18)
    if(legend) {
      plot.volcano.i <- plot.volcano.i +
        theme(legend.title = element_blank(),
              legend.key.width = unit(1, 'cm'), plot.margin = unit(c(1, 1, 1, 1.5), 'cm'))
    } else {
      plot.volcano.i <- plot.volcano.i +
        theme(legend.position = 'none', plot.margin = unit(c(1, 1, 1, 1.5), 'cm'))
    }
    plot.volcano[[i]] <- plot.volcano.i
    if(!is.null(directory)) print(plot.volcano.i)
  }
  if(!is.null(directory)) dev.off()

  return(list(plot.lfc = plot.lfc, plot.volcano = plot.volcano))
}
