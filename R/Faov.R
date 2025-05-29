#' Perform a Functional ANOVA Over a Signal

#' @param dta A numeric matrix of dimension n × T, where n is the number of samples and T the number of time points.
#' @param design A design matrix for the full model (under H1).
#' @param design0 An optional design matrix for the reduced model (under H0). Default is a column of ones.
#' @param edesign An optional contrast matrix to test specific linear hypotheses.
#' @param nbf An integer or vector of integers indicating the number of factors for the FA approximation.
#' @param pvalue Method to compute p-values. One of `"none"`, `"Satterthwaite"`, or `"MC"` (Monte Carlo). Default is `"none"`.
#' @param nsamples Number of Monte Carlo samples if `pvalue = "MC"`. Default is 200.
#' @param min.err Convergence threshold for the FA estimation algorithm. Default is 0.01.
#' @param verbose Logical. If TRUE, convergence messages are printed. Default is FALSE.
#' @param parallel Logical. If TRUE, computations are done in parallel. Default is TRUE.
#' @param nbcores Number of cores to use if `parallel = TRUE`. If NULL, all available cores are used.
#' @param sd Logical. If TRUE, residual standard deviations are returned. Default is FALSE.
#'
#' @return A list with components:
#' \describe{
#'   \item{Fols}{F-statistic using ordinary least squares (OLS).}
#'   \item{Fgls0}{F-statistic using generalized least squares (GLS) without FA modeling.}
#'   \item{Fgls}{F-statistics using GLS with FA modeling (if `nbf > 0`).}
#'   \item{F}{Pointwise OLS F-statistics (vector of length T).}
#'   \item{wF}{Pointwise GLS F-statistics (matrix of size `length(nbf)` × T).}
#'   \item{b.ols}{OLS beta coefficients.}
#'   \item{coef}{All coefficients from OLS (including non-tested effects).}
#'   \item{beta}{Coefficients for the tested effects only.}
#'   \item{sdres}{Residual standard deviations (length T).}
#' }
#'
#' @examples
#' # Simulated data
#' set.seed(123)
#' n <- 50  # number of subjects
#' T <- 100  # number of time points
#' x <- rnorm(n)
#' y <- rbinom(n, 1, 0.5)
#' X <- model.matrix(~ x + y)
#' data <- matrix(rnorm(n * T), nrow = n)
#'
#' # Run FAov with no factor modeling
#' result <- Faov(dta = data, design = X, nbf = 0, nbcores = 2)
#' result$Fols  # global OLS F-statistic
#'
#' @export
#' @importFrom parallel mclapply detectCores
#' @importFrom Matrix chol


Faov = function (dta, design, design0 = NULL, edesign = NULL, nbf = 0,
                    pvalue = c("none", "Satterthwaite", "MC"), nsamples = 200,
                    min.err = 0.01, verbose = FALSE, parallel = TRUE, nbcores = NULL,
                    sd = FALSE)
{

  # Variance of each row
  rowVars <- function(x) {
    n <- nrow(x)
    p <- ncol(x)
    m <- rowMeans(x)
    cx <- x - tcrossprod(m, rep(1, p))
    (p/(p - 1)) * rowMeans(cx^2)
  }


  F0 = function(s, u, uw = NULL, y, vid, cZ, sqrtcz, idsignal,
                rdf, rdfw = NULL, lab, nbf = 0, min.err = 1e-05, verbose = FALSE) {



    emfa <- function(dta, nbf, min.err = 1e-05, verbose = FALSE) {
      ifa <- function(Psi, B) {
        m = nrow(B)
        nbf = ncol(B)

        # Eq section 4.1
        phi = 1/sqrt(Psi)
        Phi = outer(phi, rep(1, nbf))

        # ??
        beta = B * Phi

        #L in section 4.1
        svdBeta = corpcor::fast.svd(beta)


        #theta in section 4.1
        theta = svdBeta$u * outer(rep(1, m), svdBeta$d/sqrt(1 +
                                                              svdBeta$d^2))

        phitheta = theta * Phi

        # t(beta)%*%theta
        thetabeta = crossprod(theta, beta)
        theta2beta = tcrossprod(theta, t(thetabeta))
        aux = beta - theta2beta
        iSB = Phi * aux
        return(list(iSB = iSB, Phi = phi, Theta = theta))
      }
      m = ncol(dta)
      n = nrow(dta)
      mdta = colMeans(dta)
      vdta = colMeans(dta^2) - mdta^2

      # sd at each time point
      sddta = sqrt(n/(n - 1)) * sqrt(vdta)

      #centered values
      cdta = dta - outer(rep(1, n), mdta)

      lPsi = as.list(rep(0, length(nbf)))
      lB = as.list(rep(0, length(nbf)))

      # ??
      lPsi[nbf == 0] = lapply(1:sum(nbf == 0), function(i,
                                                        v) v, v = sddta^2)
      # ??
      lB[nbf == 0] = lapply(1:sum(nbf == 0), function(i,
                                                      v) v, v = NULL)
      if (any(nbf > 0)) {
        svddta = corpcor::fast.svd(cdta/sqrt(n - 1))

        # squared singular values
        levalues = lapply(nbf[nbf > 0], function(n,
                                                 d) (d[1:n])^2, d = svddta$d)

        # left singular vectors
        levectors = lapply(nbf[nbf > 0], function(n,
                                                  v) v[, 1:n, drop = FALSE], v = svddta$v)


        # ??

        lB[nbf > 0] =  Map(function(val, vec, m) vec *
                             outer(rep(1, m), sqrt(val)), levalues, levectors,
                           m = m)
        lPsi[nbf > 0] = lapply(lB[nbf > 0], function(B,
                                                     v) v - rowSums(B^2), v = sddta^2)
        lPsi[nbf > 0] = lapply(lPsi[nbf > 0], function(Psi) {
          tmp <- Psi
          tmp[tmp <= 1e-16] = 1e-16
          tmp
        })
        crit = rep(1, sum(nbf > 0))
        lBnew = lB
        lPsinew = lPsi
        while (any(crit > min.err)) {
          not_convgd = crit > min.err
          liS = Map(ifa, lPsi[nbf > 0][not_convgd],
                    lB[nbf > 0][not_convgd])
          liSB = lapply(liS, function(x) x$iSB)
          lxiSB = lapply(liSB, function(iSB, x) x %*%
                           iSB, x = cdta)
          lCyz = lapply(lxiSB, function(xiSB, x, n) crossprod(x,
                                                              xiSB)/(n - 1), x = cdta, n = n)
          liSBCyz = Map(crossprod, liSB, lCyz)


          # Eq (7)
          ltmp = Map(function(B, iSB, nbf) diag(nbf) -
                       crossprod(B, iSB), lB[nbf > 0][not_convgd],
                     liSB, as.list(nbf[nbf > 0][not_convgd]))
          lCzz = Map("+", liSBCyz, ltmp)
          lBnew[nbf > 0][not_convgd] = Map(function(Cyz,
                                                    Czz) Cyz %*% solve(Czz), lCyz, lCzz)
          lPsinew[nbf > 0][not_convgd] = lapply(lBnew[nbf >
                                                        0][not_convgd], function(B, v) v - rowSums(B^2),
                                                v = sddta^2)
          lPsinew[nbf > 0][not_convgd] = lapply(lPsinew[nbf >
                                                          0][not_convgd], function(Psi) {
                                                            tmp <- Psi
                                                            tmp[tmp <= 1e-16] = 1e-16
                                                            tmp
                                                          })
          crit = unlist(Map(function(Psi, Psinew) mean((Psi -
                                                          Psinew)^2), lPsi[nbf > 0][not_convgd], lPsinew[nbf >
                                                                                                           0][not_convgd]))
          lB[nbf > 0] = lBnew[nbf > 0]
          lPsi[nbf > 0] = lPsinew[nbf > 0]
          if (verbose)
            print(paste("Convergence criterion: ", signif(max(crit),
                                                          digits = ceiling(-log10(min.err))), sep = ""))
        }
      }
      res = list(B = lB, Psi = lPsi, Objective = max(crit))
      return(res)
    }
    isqrtfa <- function(lPsi, lB, v) {
      m = length(lPsi[[1]])
      vnbf = unlist(lapply(lB, ncol))
      sing = Map(function(B, Psi, nbf) B/tcrossprod(sqrt(Psi),
                                                    rep(1, nbf)), lB, lPsi, vnbf)
      sing = lapply(sing, corpcor::fast.svd)
      ld = lapply(sing, function(x) x$d)
      lu = lapply(sing, function(x) x$u)
      lvphi = lapply(lPsi, function(Psi, v, k) v/tcrossprod(rep(1,
                                                                k), sqrt(Psi)), v = v, k = nrow(v))
      ltmp = Map(function(nbf, d) rep(1, nbf) - 1/sqrt(rep(1,
                                                           nbf) + d^2), vnbf, ld)
      ltmp = Map(function(u, tmp, m) u * tcrossprod(rep(1,
                                                        m), tmp), lu, ltmp, m = m)
      ltmp = Map("%*%", lvphi, ltmp)
      ltmp = Map(function(vphi, tmp, u) vphi - tcrossprod(tmp,
                                                          u), lvphi, ltmp, lu)
      return(ltmp)
    }
    n = nrow(y)
    m = ncol(y)
    p = length(idsignal)
    us = u[s, , drop = FALSE]
    tuy = crossprod(us, y)
    beta.ols = vid %*% tuy
    rownames(beta.ols) <- lab
    fit = us %*% tuy
    coeff.ols = beta.ols
    beta.ols = beta.ols[idsignal, , drop = FALSE]
    res = y - fit
    rss = colSums(res^2)
    sdres = sqrt(rss/rdf)
    scres = (res/tcrossprod(rep(1, n), sdres)) * sqrt((n -
                                                         1)/rdf)
    if (!is.null(uw)) {
      uws = uw[s, , drop = FALSE]
      tuyw = crossprod(uws, y)
      fitw = uws %*% tuyw
      resw = y - fitw
      rssw = colSums(resw^2)
      sdres = sqrt((rss - rssw)/(rdf - rdfw))
    }
    sigma = sqrt(mean(sdres^2))
    Phi = rep(1/sigma, length(sdres))
    Phibeta = beta.ols * tcrossprod(rep(1, p), Phi)
    Phibetatcz = crossprod(Phibeta, t(cZ))
    Fols = sum(Phibetatcz^2)/(T * p)
    Fgls = rep(NA, length(nbf))
    pointwise_w_F = matrix(NA, nrow = length(nbf), ncol = m)
    b.ols = (sqrtcz %*% beta.ols)/tcrossprod(rep(1, p),
                                             sdres)
    pointwise_F = colSums(b.ols^2)
    pointwise_w_F[nbf == 0, ] = tcrossprod(rep(1, sum(nbf ==
                                                        0)), pointwise_F)
    Fgls0 = sum(pointwise_F)
    Fgls[nbf == 0] = Fgls0
    if (length(idsignal) == 1) {
      pointwise_F = b.ols[1, ]
      pointwise_w_F[nbf == 0, ] = tcrossprod(rep(1, sum(nbf ==
                                                          0)), b.ols[1, ])
    }
    if (any(nbf > 0)) {
      lfa = emfa(scres, nbf = nbf[nbf > 0], min.err = min.err,
                 verbose = verbose)
      lOmega = isqrtfa(lfa$Psi, lfa$B, b.ols)
      lpointwise_w_F = lapply(lOmega, function(Omega) colSums(Omega^2))
      pointwise_w_F[nbf > 0, ] = matrix(unlist(lpointwise_w_F),
                                        nrow = sum(nbf > 0), byrow = TRUE)
      lFgls = lapply(lpointwise_w_F, sum)
      Fgls[nbf > 0] = unlist(lFgls)
      if (length(idsignal) == 1) {
        lpointwise_w_F = lapply(lOmega, function(x) x[1,
        ])
        pointwise_w_F[nbf > 0, ] = matrix(unlist(lpointwise_w_F),
                                          nrow = sum(nbf > 0), byrow = TRUE)
      }
    }
    return(list(Fols = Fols, Fgls0 = Fgls0, Fgls = Fgls,
                b.ols = b.ols, F = pointwise_F/p, wF = pointwise_w_F/p,
                beta = beta.ols, b.ols = b.ols, coef = coeff.ols,
                sdres = sdres))
  }
  if (is.null(design0))
    design0 = matrix(1, nrow = nrow(dta), ncol = 1)
  if (!is.logical(verbose))
    stop("verbose should be logical")
  erpdta_faov = as.matrix(dta)
  design = as.matrix(design)
  design0 = as.matrix(design0)
  if (!is.null(edesign))
    edesign = as.matrix(edesign)
  pvalue = match.arg(pvalue, choices = c("none", "Satterthwaite",
                                         "MC"))
  if (typeof(nsamples) != "double")
    stop("nsamples sould be an integer, usually larger than 200.")
  if (typeof(erpdta_faov) != "double")
    stop("ERPs should be of type double")
  if (nrow(erpdta_faov) != nrow(design))
    stop("dta and design should have the same number of rows")
  if (nrow(erpdta_faov) != nrow(design0))
    stop("dta and design0 should have the same number of rows")
  if (!is.null(edesign)) {
    if (nrow(erpdta_faov) != nrow(edesign))
      stop("dta and edesign should have the same number of rows")
  }
  if (ncol(design) <= ncol(design0))
    stop("design0 should have fewer columns than design")
  if (!is.null(edesign)) {
    if (ncol(edesign) <= max(ncol(design0), ncol(design)))
      stop("edesign should have more columns than design0 and design")
  }
  idsignal = NULL
  for (j in 1:ncol(design)) {
    cj = apply(design0, 2, function(x, y) all(x == y), y = design[,
                                                                  j])
    if (all(!cj))
      idsignal = c(idsignal, j)
  }
  if (length(idsignal) < (ncol(design) - ncol(design0)))
    stop("the null model design0 should be nested into the non-null model design")
  idsignalw0 = NULL
  if (!is.null(edesign)) {
    for (j in 1:ncol(edesign)) {
      cj = apply(design0, 2, function(x, y) all(x == y),
                 y = edesign[, j])
      if (all(!cj))
        idsignalw0 = c(idsignalw0, j)
    }
    if (length(idsignalw0) < (ncol(edesign) - ncol(design0)))
      stop("the null model design0 should be nested into model edesign")
  }
  idsignalw = NULL
  if (!is.null(edesign)) {
    for (j in 1:ncol(edesign)) {
      cj = apply(design, 2, function(x, y) all(x == y),
                 y = edesign[, j])
      if (all(!cj))
        idsignalw = c(idsignalw, j)
    }
    if (length(idsignalw) < (ncol(edesign) - ncol(design)))
      stop("the non-null model design should be nested into model edesign")
  }
  if ((pvalue == "Satterthwaite") & (nsamples < 200))
    stop("Since pvalue=Satterthwaite, the number of MC samples should be at least 200.")
  if (parallel & is.null(nbcores))
    nbcores = parallel::detectCores() - 1
  nbcores = min(nbcores, parallel::detectCores() - 1)
  if (parallel)
    cl = parallel::makeCluster(getOption("cl.cores", nbcores))
  n = nrow(erpdta_faov)
  T = ncol(erpdta_faov)
  svd.design = corpcor::fast.svd(design)
  svd.design0 = corpcor::fast.svd(design0)
  rdf1 = nrow(design) - length(svd.design$d)
  rdf0 = nrow(design0) - length(svd.design0$d)
  P0 = diag(n) - tcrossprod(svd.design0$u)
  pdesign = svd.design$v/tcrossprod(rep(1, ncol(design)),
                                    svd.design$d)
  pdesign = tcrossprod(pdesign, svd.design$u)
  rdfw <- NULL
  if (!is.null(edesign)) {
    svd.edesign = corpcor::fast.svd(edesign)
    rdfw = nrow(edesign) - length(svd.edesign$d)
    pedesign = svd.edesign$v/tcrossprod(rep(1, ncol(edesign)),
                                        svd.edesign$d)
    pedesign = tcrossprod(pedesign, svd.edesign$u)
  }
  Z = design[, idsignal, drop = FALSE]
  cZ = P0 %*% Z
  Szz = crossprod(cZ)/n
  svdcz = corpcor::fast.svd(cZ)
  sqrtcz = svdcz$v %*% diag(svdcz$d, nrow = length(svdcz$d),
                            ncol = length(svdcz$d)) %*% t(svdcz$v)
  vid = svd.design$v/tcrossprod(rep(1, ncol(design)), svd.design$d)
  mb.ols = NULL
  if (pvalue != "none")
    lsamples = lapply(1:nsamples, function(i, n) sample(1:n),
                      n = n)
  if (is.null(edesign)) {
    Fgls = F0(1:n, u = svd.design$u, uw = NULL, y = erpdta_faov,
              vid = vid, cZ = cZ, sqrtcz = sqrtcz, idsignal = idsignal,
              rdf = rdf1, rdfw = NULL, lab = colnames(design),
              nbf = nbf, min.err = min.err, verbose = verbose)

  }
  if (!is.null(edesign)) {
    Fgls = F0(1:n, u = svd.design$u, uw = svd.edesign$u,
              y = erpdta_faov, vid = vid, cZ = cZ, sqrtcz = sqrtcz,
              idsignal = idsignal, rdf = rdf1, rdfw = rdfw, lab = colnames(design),
              nbf = nbf, min.err = min.err, verbose = verbose)
  }
  sd_beta <- NULL
  if (sd) {
    sd_beta <- tcrossprod(sqrt(diag(solve(Szz))), Fgls$sdres)/sqrt(n)
    rownames(sd_beta) <- colnames(design)[idsignal]
  }
  pval.Fols <- NULL
  pval.Fgls <- NULL
  f0.gls <- NULL
  pval_aggregate <- NULL
  pval_aggregate_approx <- NULL
  if (pvalue != "none") {
    if (verbose)
      print("Starting Monte-Carlo estimation of the p-value")
    if (!parallel) {
      if (is.null(edesign)) {
        f0.gls = lapply(lsamples, F0, u = svd.design$u,
                        y = erpdta_faov, vid = vid, cZ = cZ, sqrtcz = sqrtcz,
                        idsignal = idsignal, rdf = rdf1, uw = NULL,
                        lab = colnames(design), rdfw = NULL, nbf = nbf,
                        min.err = min.err, verbose = FALSE)
      }
      if (!is.null(edesign)) {
        f0.gls = lapply(lsamples, F0, u = svd.design$u,
                        y = erpdta_faov, vid = vid, cZ = cZ, sqrtcz = sqrtcz,
                        idsignal = idsignal, rdf = rdf1, lab = colnames(design),
                        uw = svd.edesign$u, rdfw = rdfw, nbf = nbf,
                        min.err = min.err, verbose = FALSE)
      }
    }
    if (parallel) {
      if (is.null(edesign)) {
        f0.gls = parallel::parLapply(cl = cl, lsamples,
                                     F0, u = svd.design$u, y = erpdta_faov, vid = vid,
                                     cZ = cZ, sqrtcz = sqrtcz, idsignal = idsignal,
                                     rdf = rdf1, uw = NULL, lab = colnames(design),
                                     rdfw = NULL, nbf = nbf, min.err = min.err)
      }
      if (!is.null(edesign)) {
        f0.gls = parallel::parLapply(cl = cl, lsamples,
                                     F0, u = svd.design$u, y = erpdta_faov, vid = vid,
                                     cZ = cZ, sqrtcz = sqrtcz, idsignal = idsignal,
                                     rdf = rdf1, uw = svd.edesign$u, lab = colnames(design),
                                     rdfw = rdfw, nbf = nbf, min.err = min.err)
      }
    }
    lwF <- lapply(f0.gls, function(x) x$wF)
    f0.ols = lapply(f0.gls, function(x) x$Fols)
    f0.ols = unlist(f0.ols)
    f0.gls = lapply(f0.gls, function(x) x$Fgls)
    f0.gls = matrix(unlist(f0.gls), nrow = length(nbf))
    mf0 <- rowMeans(f0.gls)
    vf0 <- rowVars(f0.gls)
    const <- vf0/(2 * mf0)
    nu <- 2 * mf0^2/vf0
    pval_f0 <- Map(function(f, const, nu) {
      pchisq(f/const, df = nu, lower.tail = FALSE)
    }, data.frame(t(f0.gls)), const, nu)
    pval_f0 <- matrix(unlist(pval_f0), nrow = length(nbf),
                      byrow = TRUE)
    if (length(nbf) > 1) {
      sim_min_p_0 <- apply(pval_f0, 2, min)
    } else {
      ## Ajout remi
      sim_min_p_0 <- min(pval_f0)
    }
    if (pvalue == "MC") {
      pval.Fgls = rowMeans(f0.gls >= tcrossprod(Fgls$Fgls,
                                                rep(1, nsamples)))
      pval.Fols = mean(f0.ols >= Fgls$Fols)
    }
    if (pvalue == "Satterthwaite") {
      pval.Fgls = pchisq(Fgls$Fgls/const, df = nu, lower.tail = FALSE)
      const = var(f0.ols)/(2 * mean(f0.ols))
      nu = 2 * mean(f0.ols)^2/var(f0.ols)
      pval.Fols = pchisq(Fgls$Fols/const, df = nu, lower.tail = FALSE)
    }
    if (length(nbf) > 1) {
      pval_aggregate <- mean(sim_min_p_0 <= min(pval.Fgls))
      phi_min_p_0 <- qnorm(sim_min_p_0)
      pval_aggregate_approx <- pnorm(qnorm(min(pval.Fgls)),
                                     mean = mean(phi_min_p_0), sd = sd(phi_min_p_0))
    }    else { ## Ajout remi
      pval_aggregate <- mean(sim_min_p_0 <= min(pval.Fgls))
      phi_min_p_0 <- qnorm(sim_min_p_0)
      pval_aggregate_approx <- pnorm(qnorm(min(pval.Fgls)),
                                     mean = mean(phi_min_p_0), sd = sd(phi_min_p_0))
    }
  }
  if (parallel)
    parallel::stopCluster(cl)
  return(list(Fgls = Fgls$Fgls/T, Fols = Fgls$Fols, pval.Fgls = pval.Fgls,
              pval.Fols = pval.Fols, nbf = nbf, F = Fgls$F, wF = Fgls$wF,
              f0 = f0.gls/T, sd = Fgls$sdres, beta = Fgls$coef, sdbeta = sd_beta,
              pval = pval_aggregate, pval_approx = pval_aggregate_approx,
              rdf0 = rdf0, rdf1 = rdf1, erdf = rdfw))
}
