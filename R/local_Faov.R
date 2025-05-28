
# Function taken from David with custom messages to identify problems

# data_fanova_sample = data_fanova %>% slice(1:800)
# dta <- data_fanova %>% select(-c(1:3))
# dta <- log_curves
# dta <- curves
# pvalue="Satterthwaite";nsamples=1000;nbf=0:4;
# window_size=window_size;
# number_intervals=nb_intervals
# edesign =  NULL
# parallel = TRUE
# verbose = TRUE
# min.err = 0.01
# sig.level = 0.05
# design0 = NULL
# nbcores = NULL
# pvalue = "Satterthwaite"
# nsamples=1000;
# nbf=0:1
#
# dta = as.matrix(data_test[,-c(1:2)]);
# design = design_test; window_size = 5; number_intervals = 10;verbose = TRUE

my_local_Faov <- function (dta, design, design0 = NULL, edesign = NULL, nbf = 0,
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
  T = ncol(erpdta)

  # Why do we start at 1 instead of 1+window_size ?
  # ex. seq(1 + 0.5*window_size, T - 0.5*window_size, length = number_intervals)
  interval_centers <- seq(from = 1, to = T, length = number_intervals)
  pvalues <- rep(NA, length(interval_centers))

  message()
  if (verbose)
    print("Starting local fANOVA signal identification")
  for (k in 1:number_intervals) {

    # Distances of each point to the centers of the intervals
    distances = abs((1:T) - interval_centers[k])

    # l closest points to the center of the current interval
    interval = sort(order(distances)[1:window_size])

    # perform Functional ANOVA on this interval
    F = my_Faov(erpdta[, interval], design = design, design0 = design0,
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
      distances = abs((1:T) - interval_centers[k])

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

      # Add Remi: drop = FALSE (avoid pb if only one column to select)

      # + avoid pb if nb factors dependance greater than actual nb of columns to consider

      # pb here if legnth(stepwise_intervals[[1]]) == ncol(erpdta)
      # + avoid pb if first interval is the wxhole time frame
      first_interval = stepwise_intervals[[which(unlist(lapply(stepwise_intervals, length) )!= ncol(erpdta))[1]]]
      nbf_reduced = nbf[which(nbf < ncol(erpdta[, -first_interval, drop = FALSE]))]

      # F = Faov(erpdta[, -stepwise_intervals[[1]], drop = FALSE], design = design,
      #        design0 = design0, edesign = edesign, nbf = nbf_reduced,
      #        nsamples = nsamples, pvalue = pvalue, min.err = min.err,
      #        verbose = FALSE, parallel = parallel, nbcores = nbcores)


      # Perform Fanova on the whole timeframe minus the union of all significant intervals --> should be non significant !!
      F = my_Faov(erpdta[, -first_interval, drop = FALSE], design = design,
                  design0 = design0, edesign = edesign, nbf = nbf_reduced,
                  nsamples = nsamples, pvalue = pvalue, min.err = min.err,
                  verbose = FALSE, parallel = parallel, nbcores = nbcores)
      # F = Faov(erpdta[, c(21:23), drop = FALSE], design = design,
      #          design0 = design0, edesign = edesign, nbf = nbf,
      #          nsamples = nsamples, pvalue = pvalue, min.err = min.err,
      #          verbose = FALSE, parallel = parallel, nbcores = nbcores)
    }
    else {
      F = my_Faov(erpdta, design = design,
                  design0 = design0, edesign = edesign, nbf = nbf,
                  nsamples = nsamples, pvalue = pvalue, min.err = min.err,
                  verbose = FALSE, parallel = parallel, nbcores = nbcores)
      first_interval = stepwise_intervals[[which(unlist(lapply(stepwise_intervals, length) )!= ncol(erpdta))[1]]]
      nbf_reduced = nbf[which(nbf < ncol(erpdta[, -first_interval, drop = FALSE]))]

    }
    pval_agg <- rep(NA, nb_sig)
    pval_agg[1] <- F$pval
    if (F$pval < sig.level) {
      # significant <- stepwise_intervals[[1]]

      #aadd remi
      significant <- first_interval
      pval_out <- pval_agg[1]
    }


    # Perform Fanova on the whole timeframe minus some significant intervals, ranked by "significance level"
    # I.e intervals at the border of sig.level (ex. p = 0.04) are added first
    if (F$pval >= sig.level) {
      crit <- TRUE
      j <- 1
      while (crit & (j <= (nb_sig - 1))) {
        # while (crit & (j <= (nb_sig - 1))) {
        j <- j + 1
        interval <- stepwise_intervals[[j]]

        # Idea: keep enough points in the time frame considered ?
        if ((T - length(interval)) >= 20) {


          # + avoid pb if nb factors dependance greater than actual nb of columns to consider
          nbf_reduced = nbf[which(nbf < ncol(erpdta[, -interval, drop = FALSE]))]
          F = my_Faov(erpdta[, -interval], design = design,
                      design0 = design0, edesign = edesign, nbf = nbf_reduced,
                      nsamples = nsamples, pvalue = pvalue, min.err = min.err,
                      verbose = FALSE, parallel = parallel, nbcores = nbcores)
          pval_agg[j] <- F$pval
        }
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
    F = my_Faov(erpdta[, significant], design = design, design0 = design0,
                edesign = edesign, nbf = nbf_reduced, nsamples = nsamples,
                pvalue = pvalue, min.err = min.err, verbose = FALSE,
                parallel = parallel, nbcores = nbcores)
    pval_in <- F$pval
  }
  return(list(interval_centers = interval_centers, pvalues = pvalues,
              pvalues_agg = pval_agg, significant = significant, pval_out = pval_out,
              pval_in = pval_in))
}
