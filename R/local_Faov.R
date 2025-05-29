#' Perform Local Functional ANOVA Over a Signal
#'
#' This function identifies localized time intervals within a multivariate signal where a functional ANOVA (fANOVA) reveals statistically significant differences between experimental conditions.
#'
#' @param dta A numeric matrix of dimension \eqn{n \times T}, where each row corresponds to a subject and each column to a time point.
#' @param design A model matrix for the full (non-null) model (dimensions: \eqn{n \times p}).
#' @param design0 A model matrix for the null model (dimensions: \eqn{n \times q}, with \eqn{q < p}). If `NULL`, defaults to an intercept-only model.
#' @param edesign Optional: a model matrix for an extended model (typically for testing nested model assumptions). Must have more columns than both `design` and `design0`.
#' @param nbf A vector of indices indicating which effects in the model correspond to the factor of interest. These indices are used for F-statistic computation.
#' @param pvalue Method for approximating p-values: `"none"`, `"Satterthwaite"` (default), or `"MC"` for Monte Carlo approximation.
#' @param nsamples Integer. Number of Monte Carlo samples used when computing p-values (for `"Satterthwaite"` or `"MC"` methods).
#' @param min.err Minimum error for iterative computations (relevant when approximating p-values).
#' @param verbose Logical. If `TRUE`, prints progress messages during the computation.
#' @param parallel Logical. If `TRUE`, computations are parallelized over time intervals.
#' @param nbcores Number of CPU cores to use for parallelization. If `NULL`, defaults to the number of available cores minus one.
#' @param window_size Integer. The number of time points included in each local interval (i.e., the bandwidth).
#' @param number_intervals Integer. Number of intervals to test along the time axis.
#' @param sig.level Significance level used to determine whether a time interval is considered significant.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{`interval_centers`}{Vector of center indices for each time interval tested.}
#'   \item{`pvalues`}{Vector of p-values (one per tested interval).}
#'   \item{`pvalues_agg`}{Vector of aggregated p-values from global tests after removing stepwise significant intervals.}
#'   \item{`significant`}{Vector of time indices considered significant across intervals.}
#'   \item{`pval_out`}{p-value of the global test excluding significant intervals.}
#'   \item{`pval_in`}{p-value of the global test restricted to significant intervals.}
#' }
#'
#' @details
#' The function applies a sliding-window fANOVA along the time dimension of the signal to identify regions with statistically significant differences. It controls for multiple testing by performing a stepwise removal of significant regions and evaluating the global null hypothesis on the remaining signal.
#'
#' The function assumes that `design0` is nested within `design`, and that if provided, `edesign` nests both `design0` and `design`.
#'
#' Parallel computation can substantially speed up the local testing across time intervals.
#'
#' @seealso [Faov()] for the underlying fANOVA test.
#'
#' @examples
#' \dontrun{
#' # Simulated example with 100 subjects and 200 time points
#' n <- 100
#' T <- 200
#' signal <- matrix(rnorm(n * T), nrow = n)
#' group <- rep(1:2, each = n/2)
#' design <- model.matrix(~ factor(group))
#' res <- local_Faov(signal, design = design, pvalue = "Satterthwaite", nbcores = 2)
#' plot(res$interval_centers, res$pvalues, type = "b")
#' }
#'
#' @export
#' @importFrom parallel mclapply detectCores

local_Faov <- function (dta, design, design0 = NULL, edesign = NULL, nbf = 0,
                           pvalue = c("none", "Satterthwaite", "MC"), nsamples = 200,
                           min.err = 0.01, verbose = FALSE, parallel = TRUE, nbcores = NULL,
                           window_size = round(ncol(dta)/10), number_intervals = round(ncol(dta)/5),
                           sig.level = 0.05)
{
  if (is.null(design0))
    design0 = matrix(1, nrow = nrow(dta), ncol = 1)
  if (!is.logical(verbose))
    stop("verbose should be logical")
  erpdta = as.matrix(dta)
  design = as.matrix(design)
  design0 = as.matrix(design0)
  if (!is.null(edesign))
    edesign = as.matrix(edesign)
  pvalue = match.arg(pvalue, choices = c("none", "Satterthwaite",
                                         "MC"))
  if (typeof(nsamples) != "double")
    stop("nsamples sould be an integer, usually larger than 200.")
  if (typeof(erpdta) != "double")
    stop("ERPs should be of type double")
  if (nrow(erpdta) != nrow(design))
    stop("dta and design should have the same number of rows")
  if (nrow(erpdta) != nrow(design0))
    stop("dta and design0 should have the same number of rows")
  if (!is.null(edesign)) {
    if (nrow(erpdta) != nrow(edesign))
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
  n = nrow(erpdta)
  length_grid = ncol(erpdta)

  # Why do we start at 1 instead of 1+window_size ?
  # ex. seq(1 + 0.5*window_size, length_grid - 0.5*window_size, length = number_intervals)
  # interval_centers <- seq(from = 1, to = length_grid, length = number_intervals)

  interval_centers <- seq(from = 1 + 0.5*window_size, to = length_grid - 0.5*window_size,
                          length = number_intervals)

  if(number_intervals * window_size < length_grid) {
    warning("The whole time frame will not be covered (number_intervals * window_size < length_grid),
consider increasing number_intervals or window_size")
  }

  pvalues <- rep(NA, length(interval_centers))

  message()
  if (verbose)
    print("Starting local fANOVA signal identification")
  for (k in 1:number_intervals) {

    # Distances of each point to the centers of the intervals
    distances = abs((1:length_grid) - interval_centers[k])

    # l closest points to the center of the current interval
    interval = sort(order(distances)[1:window_size])

    # perform Functional ANOVA on this interval
    F = Faov(erpdta[, interval], design = design, design0 = design0,
                edesign = edesign, nbf = nbf, nsamples = nsamples,
                pvalue = "Satterthwaite", min.err = min.err, verbose = FALSE,
                parallel = parallel, nbcores = nbcores)
    pvalues[k] = F$pval_approx
    if (verbose & k %% 20 == 0)
      print(paste(k, "th local fANOVA p-value calculated over ",
                  number_intervals, sep = ""))
  }
  if (min(pvalues) >= sig.level) {
    significant <- integer(0)
    pval_agg <- NULL # Modification code David
    pval_out <- NULL
    pval_in <- NULL
  }
  if (min(pvalues) <= sig.level) {
    nb_sig <- sum(pvalues <= sig.level)


    ord <- order(pvalues, decreasing = FALSE)
    whole_interval <- integer(0)
    stepwise_intervals <- as.list(rep(0, nb_sig))
    j <- 0
    for (k in ord[1:nb_sig]) {
      j <- j + 1
      distances = abs((1:length_grid) - interval_centers[k])

      # take the points around the center of the interval
      # bandwith = window_size
      interval = sort(order(distances)[1:window_size])
      whole_interval <- union(whole_interval, interval)
      stepwise_intervals[[j]] <- whole_interval
    }
    stepwise_intervals <- stepwise_intervals[nb_sig:1]

    # Ajout / modification Rémi: avoid the fact that if all intervals
    # are significant, then pb with Faov (because of the removal of all the interval)
    if(nb_sig < number_intervals){



      # pb here if legnth(stepwise_intervals[[1]]) == ncol(erpdta)
      # + avoid pb if first interval is the wxhole time frame
      first_interval = stepwise_intervals[[which(unlist(lapply(stepwise_intervals,
                                                               length) )!= ncol(erpdta))[1]]]

      # Add Remi: drop = FALSE (avoid pb if only one column to select)
      # + avoid pb if nb factors dependance greater than actual nb of columns to consider

      nbf_reduced = nbf[which(nbf < ncol(erpdta[, -first_interval, drop = FALSE]))]

      # F = Faov(erpdta[, -stepwise_intervals[[1]], drop = FALSE], design = design,
      #        design0 = design0, edesign = edesign, nbf = nbf_reduced,
      #        nsamples = nsamples, pvalue = pvalue, min.err = min.err,
      #        verbose = FALSE, parallel = parallel, nbcores = nbcores)


      # Perform Fanova on the whole timeframe minus the union of
      # all significant intervals --> should be non significant !!
      F = Faov(erpdta[, -first_interval, drop = FALSE], design = design,
                  design0 = design0, edesign = edesign, nbf = nbf_reduced,
                  nsamples = nsamples, pvalue = pvalue, min.err = min.err,
                  verbose = FALSE, parallel = parallel, nbcores = nbcores)
    }
    else {
      F = Faov(erpdta, design = design,
                  design0 = design0, edesign = edesign, nbf = nbf,
                  nsamples = nsamples, pvalue = pvalue, min.err = min.err,
                  verbose = FALSE, parallel = parallel, nbcores = nbcores)
      first_interval = stepwise_intervals[[which(unlist(lapply(stepwise_intervals, length) )!= ncol(erpdta))[1]]]
      nbf_reduced = nbf[which(nbf < ncol(erpdta[, -first_interval, drop = FALSE]))]

    }
    pval_agg <- rep(NA, nb_sig)
    pval_agg[1] <- F$pval
    if (F$pval < sig.level) {

      #add remi
      significant <- first_interval
      pval_out <- pval_agg[1]
    }


    # Perform Fanova on the whole timeframe minus some significant intervals,
    #ranked by "significance level"
    # I.e intervals at the border of sig.level (ex. p = 0.04) are added first
    if (F$pval >= sig.level) {
      crit <- TRUE
      j <- 1
      while (crit & (j <= (nb_sig - 1))) {
        # while (crit & (j <= (nb_sig - 1))) {
        j <- j + 1
        interval <- stepwise_intervals[[j]]

        # Idea: keep enough points in the time frame considered ?
        # if ((length_grid - length(interval)) >= 20) {


          # + avoid pb if nb factors dependance greater than actual nb of columns to consider
          nbf_reduced = nbf[which(nbf < ncol(erpdta[, -interval, drop = FALSE]))]
          F = Faov(erpdta[, -interval], design = design,
                      design0 = design0, edesign = edesign, nbf = nbf_reduced,
                      nsamples = nsamples, pvalue = pvalue, min.err = min.err,
                      verbose = FALSE, parallel = parallel, nbcores = nbcores)
          pval_agg[j] <- F$pval
        # }
        crit <- (F$pval >= sig.level)
        if (crit & verbose)
          print(paste("Number of significant time points: ",
                      length(interval), sep = ""))
      }
      if(nb_sig == 1) {j = 2} # AJOT REMI
      significant <- stepwise_intervals[[j - 1]]
      pval_out <- pval_agg[j - 1]
    }
    # + avoid pb if nb factors dependance greater than actual nb of columns to consider
    F = Faov(erpdta[, significant], design = design, design0 = design0,
                edesign = edesign, nbf = nbf_reduced, nsamples = nsamples,
                pvalue = pvalue, min.err = min.err, verbose = FALSE,
                parallel = parallel, nbcores = nbcores)
    pval_in <- F$pval
  }
  return(list(interval_centers = interval_centers, pvalues = pvalues,
              pvalues_agg = pval_agg, significant = significant, pval_out = pval_out,
              pval_in = pval_in))
}
