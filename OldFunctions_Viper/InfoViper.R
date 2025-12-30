> viper
function (eset, regulon, dnull = NULL, pleiotropy = FALSE, nes = TRUE, 
          method = c("none", "scale", "rank", "mad", "ttest"), bootstraps = 0, 
          minsize = 25, adaptive.size = FALSE, eset.filter = TRUE, 
          mvws = 1, pleiotropyArgs = list(regulators = 0.05, shadow = 0.05, 
                                          targets = 10, penalty = 20, method = "adaptive"), cores = 1, 
          verbose = TRUE) 
{
  method <- match.arg(method)
  pdata <- NULL
  if (is(eset, "viperSignature")) {
    dnull <- eset$nullmodel
    eset <- eset$signature
    method = "none"
    if (bootstraps > 0) {
      bootstraps <- 0
      warning("Using a null model, bootstraps iterations are ignored.", 
              call. = FALSE)
    }
  }
  if (pleiotropy & bootstraps > 0) {
    bootstraps <- 0
    warning("Using pleiotropic correction, bootstraps iterations are ignored.", 
            call. = FALSE)
  }
  if (is(eset, "ExpressionSet")) {
    pdata <- phenoData(eset)
    eset <- exprs(eset)
  }
  else if (is.data.frame(eset)) {
    eset <- as.matrix(eset)
  }
  if (is.null(nrow(eset))) 
    eset <- matrix(eset, length(eset), 1, dimnames = list(names(eset), 
                                                          NULL))
  if (verbose) 
    message("\nComputing the association scores")
  if (names(regulon[[1]])[1] == "tfmode") 
    regulon <- list(regulon = regulon)
  if (bootstraps > 0) {
    return(bootstrapViper(eset = eset, regulon = regulon, 
                          nes = nes, bootstraps = bootstraps, eset.filter = eset.filter, 
                          adaptive.size = adaptive.size, minsize = minsize, 
                          cores = cores, verbose = verbose))
  }
  cores1 <- 1
  if (length(regulon) > 1) {
    cores1 <- min(cores, length(regulon))
    cores <- 1
    nes <- TRUE
  }
  if (cores > 1 | cores1 > 1) 
    verbose <- FALSE
  switch(method, scale = {
    tt <- t(scale(t(eset)))
  }, rank = {
    tt <- t(apply(eset, 1, rank)) * punif(length(eset), -0.1, 
                                          0.1)
  }, mad = {
    tt <- t(apply(eset, 1, function(x) (x - median(x))/mad(x)))
  }, ttest = {
    tt <- sapply(1:ncol(eset), function(i, eset) rowTtest(eset[, 
                                                               i] - eset[, -i])$statistic, eset = eset)
    colnames(tt) <- colnames(eset)
    rownames(tt) <- rownames(eset)
  }, none = {
    tt <- eset
  })
  if (verbose) 
    message("Computing regulons enrichment with aREA")
  res <- mclapply(regulon, function(regulon, dnull, pleiotropy, 
                                    nes, tt, eset.filter, adaptive.size, cores, pleiotropyArgs, 
                                    verbose) {
    if (eset.filter) {
      tmp <- c(names(regulon), unlist(lapply(regulon, function(x) names(x$tfmode)), 
                                      use.names = FALSE))
      tt <- tt[rownames(tt) %in% unique(tmp), ]
    }
    regulon <- lapply(regulon, function(x, genes) {
      filtro <- names(x$tfmode) %in% genes
      x$tfmode <- x$tfmode[filtro]
      if (length(x$likelihood) == length(filtro)) 
        x$likelihood <- x$likelihood[filtro]
      return(x)
    }, genes = rownames(tt))
    if (adaptive.size) 
      regulon <- regulon[sapply(regulon, function(x) {
        sum((x$likelihood/max(x$likelihood))^2)
      }) >= minsize]
    else regulon <- regulon[sapply(regulon, function(x) length(x$tfmode)) >= 
                              minsize]
    es <- aREA(tt, regulon, cores = cores, minsize = 0, verbose = verbose)
    if (!nes) {
      if (pleiotropy) 
        warning("No pleiotropy correction implemented when raw es is returned.", 
                call. = FALSE)
      return(es$es)
    }
    if (is.null(dnull)) 
      nes <- es$nes
    else {
      if (verbose) 
        message("\nEstimating NES with null model")
      tmp <- aREA(dnull, regulon, cores = cores, minsize = 0, 
                  verbose = verbose)$es
      if (ncol(tmp) > 499) {
        nes <- t(sapply(1:nrow(tmp), function(i, tmp, 
                                              es) {
          aecdf(tmp[i, ], symmetric = TRUE)(es[i, ])$nes
        }, tmp = tmp, es = es$es))
        rownames(nes) <- rownames(es$nes)
      }
      else {
        nes <- es$es/sqrt(frvarna(tmp)[, 1])
      }
    }
    if (pleiotropy) {
      pb <- NULL
      if (verbose) {
        message("\nComputing pleiotropy for ", ncol(nes), 
                " samples.")
        message("\nProcess started at ", date())
      }
      if (cores > 1) {
        nes <- mclapply(1:ncol(nes), function(i, ss, 
                                              nes, regulon, args, dnull) {
          nes <- nes[, i]
          sreg <- shadowRegulon(ss[, i], nes, regulon, 
                                regulators = args[[1]], shadow = args[[2]], 
                                targets = args[[3]], penalty = args[[4]], 
                                method = args[[5]])
          if (!is.null(sreg)) {
            if (is.null(dnull)) 
              tmp <- aREA(ss[, i], sreg, minsize = 5, 
                          cores = 1)$nes[, 1]
            else {
              tmp <- aREA(cbind(ss[, i], dnull), sreg, 
                          minsize = 5, cores = 1)$es
              tmp <- apply(tmp, 1, function(x) aecdf(x[-1], 
                                                     symmetric = TRUE)(x[1])$nes)
            }
            nes[match(names(tmp), names(nes))] <- tmp
          }
          return(nes)
        }, ss = tt, nes = nes, regulon = regulon, args = pleiotropyArgs, 
        dnull = dnull, mc.cores = cores)
        nes <- sapply(nes, function(x) x)
      }
      else {
        if (verbose) 
          pb <- txtProgressBar(max = ncol(nes), style = 3)
        nes <- sapply(1:ncol(nes), function(i, ss, nes, 
                                            regulon, args, dnull, pb) {
          nes <- nes[, i]
          sreg <- shadowRegulon(ss[, i], nes, regulon, 
                                regulators = args[[1]], shadow = args[[2]], 
                                targets = args[[3]], penalty = args[[4]], 
                                method = args[[5]])
          if (!is.null(sreg)) {
            if (is.null(dnull)) 
              tmp <- aREA(ss[, i], sreg, minsize = 5)$nes[, 
                                                          1]
            else {
              tmp <- aREA(cbind(ss[, i], dnull), sreg, 
                          minsize = 5)$es
              tmp <- apply(tmp, 1, function(x) aecdf(x[-1], 
                                                     symmetric = TRUE)(x[1])$nes)
            }
            nes[match(names(tmp), names(nes))] <- tmp
          }
          if (is(pb, "txtProgressBar")) 
            setTxtProgressBar(pb, i)
          return(nes)
        }, ss = tt, nes = nes, regulon = regulon, args = pleiotropyArgs, 
        dnull = dnull, pb = pb)
      }
      if (verbose) 
        message("\nProcess ended at ", date(), "\n")
      if (is.null(nrow(nes))) 
        nes <- matrix(nes, length(nes), 1, dimnames = list(names(nes), 
                                                           NULL))
      colnames(nes) <- colnames(eset)
      return(nes)
    }
    return(nes)
  }, dnull = dnull, pleiotropy = pleiotropy, nes = nes, tt = tt, 
  eset.filter = eset.filter, adaptive.size = adaptive.size, 
  cores = cores, pleiotropyArgs = pleiotropyArgs, verbose = verbose, 
  mc.cores = cores1)
  if (length(res) == 1) 
    nes <- res[[1]]
  else {
    genes <- unique(unlist(lapply(res, rownames), use.names = FALSE))
    nes <- sapply(res, function(x, genes) as.vector(x[match(genes, 
                                                            rownames(x)), ]), genes = genes)
    nes[is.na(nes)] <- 0
    if (length(mvws) == 1) {
      ws <- abs(nes)^mvws
    }
    else {
      ws <- sigT(abs(nes), mvws[2], mvws[1])
    }
    nes <- matrix(rowSums(nes * ws)/rowSums(ws), length(genes), 
                  ncol(res[[1]]), dimnames = list(genes, colnames(res[[1]])))
  }
  if (is.null(pdata)) 
    return(nes)
  return(ExpressionSet(assayData = nes, phenoData = pdata))
}
<bytecode: 0x30c0a1350>
  <environment: namespace:viper>