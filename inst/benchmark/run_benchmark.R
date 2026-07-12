#!/usr/bin/env Rscript
## =====================================================================
## run_benchmark.R
## Reproducible benchmark: fastsimplexreg vs simplexreg vs betareg.
##
## Produced by "Dr. Bench". Fair, cold-fit, guarded benchmark covering
## five scenarios (accuracy, scaling, links, realdata, mixed).
##
## Fairness:
##  * identical simulated / real data fed to every package in a cell;
##  * cold fits (fresh call each replication); median/min/IQR from reps;
##  * fastsimplexreg timed at n_threads = 1 AND n_threads = 4 separately;
##  * matched mean + dispersion/precision model structure across packages;
##  * every competitor call guarded by requireNamespace() + tryCatch();
##  * a per-fit wall-clock budget skips a slow package at larger n.
##
## simplexreg and fastsimplexreg fit the SAME (simplex) distribution, so
## they are compared on speed AND accuracy (coefficients / logLik).
## betareg fits the BETA distribution: a fair SPEED comparator only; its
## logLik is NOT comparable across families (still recorded per package).
##
## NOTE ON NAMESPACES: simplexreg exports its own rsimplex()/simplex
## helpers that MASK fastsimplexreg's when attached. This script NEVER
## attaches simplexreg/betareg; it always calls them as pkg::fun() and
## always calls the data generators as fastsimplexreg::rsimplex() etc.
## =====================================================================

suppressPackageStartupMessages(library(fastsimplexreg))

## ---- output location -------------------------------------------------
SCRATCH <- Sys.getenv(
  "BENCH_OUTDIR",
  "/tmp/claude-1000/-home-jlopes-Dropbox-Pacotes-Dev-fastsimplexreg/f885e592-d7ea-4441-bbf9-1953954a8f52/scratchpad"
)
dir.create(SCRATCH, showWarnings = FALSE, recursive = TRUE)

## ---- reproducibility -------------------------------------------------
set.seed(20260712)
BUDGET_SEC <- 45          # per-single-fit wall-clock budget (scaling)
N_THREADS_MULTI <- 4L     # realistic multi-thread setting
LAB_FAST1 <- "fastsimplexreg"
LAB_FAST4 <- "fastsimplexreg (4 threads)"

HAS_SR <- requireNamespace("simplexreg", quietly = TRUE)
HAS_BR <- requireNamespace("betareg",   quietly = TRUE)
message("simplexreg available: ", HAS_SR, " ; betareg available: ", HAS_BR)

t_start_all <- proc.time()[["elapsed"]]

## =====================================================================
## Helpers
## =====================================================================

## Extractors: return list(loglik, npar, converged, coefs=df(term,estimate))
extract_fast <- function(fit) {
  m <- fit$coefficients$mean
  d <- fit$coefficients$dispersion
  list(
    loglik    = as.numeric(fit$logLik),
    npar      = length(fit$par),
    converged = isTRUE(fit$convergence == 0),
    coefs     = rbind(
      data.frame(term = paste0("mean:", names(m)), estimate = unname(m),
                 stringsAsFactors = FALSE),
      data.frame(term = paste0("disp:", names(d)), estimate = unname(d),
                 stringsAsFactors = FALSE)
    )
  )
}

extract_sr <- function(fit) {
  b  <- fit$fixef[, "beta"]
  bn <- rownames(fit$fixef)
  if (!is.null(fit$dispar)) {
    ## variable dispersion: log-link gamma, matches fastsimplexreg
    a  <- fit$dispar[, "alpha"]
    an <- rownames(fit$dispar)
  } else {
    ## constant dispersion: simplexreg stores the (natural-scale) dispersion in
    ## $Dispersion; convert to log scale to line up with fastsimplexreg's log
    ## link intercept. NB simplexreg estimates the constant dispersion by a
    ## separate (moment/deviance) route, so it need NOT equal the fast MLE.
    a  <- log(as.numeric(fit$Dispersion))
    an <- "(Intercept)"
  }
  ll <- tryCatch(as.numeric(stats::logLik(fit)), error = function(e) NA_real_)
  list(
    loglik    = ll,
    npar      = length(b) + length(a),
    converged = all(is.finite(c(b, a))) && is.finite(ll),
    coefs     = rbind(
      data.frame(term = paste0("mean:", bn), estimate = unname(b),
                 stringsAsFactors = FALSE),
      data.frame(term = paste0("disp:", an), estimate = unname(a),
                 stringsAsFactors = FALSE)
    )
  )
}

extract_br <- function(fit) {
  m  <- fit$coefficients$mean
  p  <- fit$coefficients$precision
  ll <- tryCatch(as.numeric(stats::logLik(fit)), error = function(e) NA_real_)
  list(
    loglik    = ll,
    npar      = length(m) + length(p),
    converged = isTRUE(fit$converged),
    coefs     = rbind(
      data.frame(term = paste0("mean:", names(m)), estimate = unname(m),
                 stringsAsFactors = FALSE),
      data.frame(term = paste0("prec:", names(p)), estimate = unname(p),
                 stringsAsFactors = FALSE)
    )
  )
}

extract_mixed <- function(fit) {
  list(
    loglik    = as.numeric(fit$logLik),
    npar      = length(fit$par),
    converged = isTRUE(fit$convergence == 0),
    coefs     = NULL
  )
}

## Cold-fit timing over `reps`; stops early if a single fit exceeds budget.
## Returns list(times = numeric with NA on failure, fit = last good fit,
## note = "", over = logical).
time_and_fit <- function(expr_fun, reps, budget = Inf) {
  times  <- rep(NA_real_, reps)
  fitobj <- NULL
  note   <- ""
  over   <- FALSE
  for (r in seq_len(reps)) {
    fr <- NULL
    tt <- tryCatch(
      system.time(fr <- expr_fun())[["elapsed"]],
      error = function(e) { note <<- paste0("error: ", conditionMessage(e)); NA_real_ }
    )
    if (is.na(tt)) break
    times[r] <- tt
    fitobj   <- fr
    if (tt > budget) {
      over <- TRUE
      if (!nzchar(note))
        note <- sprintf("over budget: %.1fs > %.0fs", tt, budget)
      break
    }
  }
  list(times = times, fit = fitobj, note = note, over = over)
}

## Build tidy fit-level rows for one cell.
make_rows <- function(scenario, package, dataset, n, link, nAGQ, J,
                      times, ext, note) {
  reps <- length(times)
  data.frame(
    scenario  = scenario,
    package   = package,
    dataset   = dataset,
    n         = as.integer(n),
    link      = link,
    nAGQ      = as.integer(nAGQ),
    J         = as.integer(J),
    rep       = seq_len(reps),
    time_sec  = as.numeric(times),
    loglik    = if (is.null(ext)) NA_real_    else as.numeric(ext$loglik),
    npar      = if (is.null(ext)) NA_integer_ else as.integer(ext$npar),
    converged = if (is.null(ext)) FALSE       else isTRUE(ext$converged),
    note      = note,
    stringsAsFactors = FALSE
  )
}

## Run a cell and return list(rows = df, ext = extraction or NULL, over).
run_cell <- function(scenario, package, dataset, n, link, nAGQ, J,
                     fit_fun, extract_fun, reps, budget = Inf) {
  res   <- time_and_fit(fit_fun, reps, budget)
  valid <- res$times[!is.na(res$times)]
  if (length(valid) == 0L || is.null(res$fit)) {
    note <- if (nzchar(res$note)) res$note else "no successful fit"
    rows <- make_rows(scenario, package, dataset, n, link, nAGQ, J,
                      NA_real_, NULL, note)
    return(list(rows = rows, ext = NULL, over = res$over))
  }
  ext  <- tryCatch(extract_fun(res$fit), error = function(e) {
    message("  extraction failed [", package, " ", dataset, " n=", n,
            "]: ", conditionMessage(e))
    NULL
  })
  note <- res$note
  if (is.null(ext) && !nzchar(note)) note <- "extraction failed"
  rows <- make_rows(scenario, package, dataset, n, link, nAGQ, J,
                    valid, ext, note)
  list(rows = rows, ext = ext, over = res$over)
}

## Emit a "skipped: over budget" placeholder row.
skip_row <- function(scenario, package, dataset, n, link, nAGQ = NA, J = NA) {
  make_rows(scenario, package, dataset, n, link, nAGQ, J,
            NA_real_, NULL, "skipped: over budget")
}

reps_for_n <- function(n) {
  if (n <= 5000L)      10L
  else if (n <= 20000L) 5L
  else if (n <= 1e5)    3L
  else                  2L
}

## Simulate the standard simplex data set: logit mean (2 covariates),
## log dispersion (1 covariate). Uses fastsimplexreg's generators.
sim_data <- function(n, link = "logit", seed = NULL,
                     b = c(-0.4, 0.8, -0.5), g = c(-1.0, 0.6)) {
  if (!is.null(seed)) set.seed(seed)
  x1 <- stats::rnorm(n)
  x2 <- stats::rbinom(n, 1, 0.4)
  z1 <- stats::rnorm(n)
  eta <- b[1] + b[2] * x1 + b[3] * x2
  mu  <- fastsimplexreg::simplex_linkinv(eta, link = link)
  phi <- exp(g[1] + g[2] * z1)
  y   <- fastsimplexreg::rsimplex(n, mu, phi)
  data.frame(y = y, x1 = x1, x2 = x2, z1 = z1)
}

RES_LIST   <- list()   # fit-level rows
COEF_LIST  <- list()   # coefficient rows
add_res  <- function(df) RES_LIST[[length(RES_LIST) + 1L]]  <<- df
add_coef <- function(df) COEF_LIST[[length(COEF_LIST) + 1L]] <<- df

coef_rows <- function(scenario, dataset, n, link, package, ext) {
  if (is.null(ext) || is.null(ext$coefs)) return(NULL)
  data.frame(
    scenario = scenario, dataset = dataset, n = as.integer(n),
    link = link, package = package,
    term = ext$coefs$term, estimate = ext$coefs$estimate,
    stringsAsFactors = FALSE
  )
}

## =====================================================================
## SCENARIO 1 — accuracy (fastsimplexreg vs simplexreg, same simplex model)
## =====================================================================
message("== Scenario: accuracy ==")
for (n in c(500L, 2000L, 8000L)) {
  reps <- reps_for_n(n)
  d <- sim_data(n, "logit", seed = 100L + n)

  ## fastsimplexreg, 1 thread
  c1 <- run_cell("accuracy", LAB_FAST1, "simulated", n, "logit", NA, NA,
                 function() fastsimplexreg(y ~ x1 + x2 | z1, data = d,
                                           link = "logit", n_threads = 1L),
                 extract_fast, reps)
  add_res(c1$rows)
  add_coef(coef_rows("accuracy", "simulated", n, "logit", LAB_FAST1, c1$ext))

  ## fastsimplexreg, 4 threads
  c4 <- run_cell("accuracy", LAB_FAST4, "simulated", n, "logit", NA, NA,
                 function() fastsimplexreg(y ~ x1 + x2 | z1, data = d,
                                           link = "logit",
                                           n_threads = N_THREADS_MULTI),
                 extract_fast, reps)
  add_res(c4$rows)

  ## simplexreg
  if (HAS_SR) {
    cs <- run_cell("accuracy", "simplexreg", "simulated", n, "logit", NA, NA,
                   function() simplexreg::simplexreg(y ~ x1 + x2 | z1,
                                                     data = d, link = "logit"),
                   extract_sr, reps)
    add_res(cs$rows)
    add_coef(coef_rows("accuracy", "simulated", n, "logit", "simplexreg", cs$ext))
  }
}

## =====================================================================
## SCENARIO 2 — scaling (workhorse): fast[1], fast[4], simplexreg, betareg
## =====================================================================
message("== Scenario: scaling ==")
scaling_n <- c(200L, 1000L, 5000L, 20000L, 100000L, 500000L)

## per-package "dead" state (skip at larger n once over budget / failed)
dead <- c(fast1 = FALSE, fast4 = FALSE, sr = FALSE, br = FALSE)

for (n in scaling_n) {
  reps <- reps_for_n(n)
  d <- sim_data(n, "logit", seed = 200L + (n %% 100000L))

  ## ---- fastsimplexreg 1 thread (never skipped; must run through) ----
  c1 <- run_cell("scaling", LAB_FAST1, "simulated", n, "logit", NA, NA,
                 function() fastsimplexreg(y ~ x1 + x2 | z1, data = d,
                                           link = "logit", n_threads = 1L),
                 extract_fast, reps, budget = Inf)
  add_res(c1$rows)

  ## ---- fastsimplexreg 4 threads (never skipped) ----
  c4 <- run_cell("scaling", LAB_FAST4, "simulated", n, "logit", NA, NA,
                 function() fastsimplexreg(y ~ x1 + x2 | z1, data = d,
                                           link = "logit",
                                           n_threads = N_THREADS_MULTI),
                 extract_fast, reps, budget = Inf)
  add_res(c4$rows)

  ## ---- simplexreg (budgeted) ----
  if (HAS_SR) {
    if (dead["sr"]) {
      add_res(skip_row("scaling", "simplexreg", "simulated", n, "logit"))
    } else {
      cs <- run_cell("scaling", "simplexreg", "simulated", n, "logit", NA, NA,
                     function() simplexreg::simplexreg(y ~ x1 + x2 | z1,
                                                       data = d, link = "logit"),
                     extract_sr, reps, budget = BUDGET_SEC)
      add_res(cs$rows)
      if (cs$over || is.null(cs$ext)) dead["sr"] <- TRUE
    }
  }

  ## ---- betareg (BETA model; speed comparator only; budgeted) ----
  if (HAS_BR) {
    if (dead["br"]) {
      add_res(skip_row("scaling", "betareg", "simulated", n, "logit"))
    } else {
      cb <- run_cell("scaling", "betareg", "simulated", n, "logit", NA, NA,
                     function() betareg::betareg(y ~ x1 + x2 | z1, data = d,
                                                 link = "logit", link.phi = "log"),
                     extract_br, reps, budget = BUDGET_SEC)
      add_res(cb$rows)
      if (cb$over || is.null(cb$ext)) dead["br"] <- TRUE
    }
  }
}

## =====================================================================
## SCENARIO 3 — links (fastsimplexreg vs simplexreg): all four mean links
## =====================================================================
message("== Scenario: links ==")
n_link <- 3000L
for (lk in c("logit", "probit", "cloglog", "neglog")) {
  ## data simulated under the matching link (fair, same for both packages)
  d <- sim_data(n_link, lk, seed = 300L + which(c("logit","probit","cloglog","neglog") == lk),
                b = c(0.2, 0.6, -0.4))
  reps <- 10L

  cl1 <- run_cell("links", LAB_FAST1, "simulated", n_link, lk, NA, NA,
                  function() fastsimplexreg(y ~ x1 + x2 | z1, data = d,
                                            link = lk, n_threads = 1L),
                  extract_fast, reps)
  add_res(cl1$rows)

  if (HAS_SR) {
    cls <- run_cell("links", "simplexreg", "simulated", n_link, lk, NA, NA,
                    function() simplexreg::simplexreg(y ~ x1 + x2 | z1,
                                                      data = d, link = lk),
                    extract_sr, reps)
    add_res(cls$rows)
  }
}

## =====================================================================
## SCENARIO 4 — realdata: GasolineYield, ReadingSkills, FoodExpenditure
##   constant-dispersion 3-way match (fast, simplexreg, betareg).
## =====================================================================
message("== Scenario: realdata ==")
realdata_cor <- list()  # cor(fitted_mu) fast vs simplexreg, stored in meta

## squeeze boundary (0/1) values into the open interval if present.
squeeze01 <- function(y) {
  N <- length(y)
  if (any(y <= 0 | y >= 1)) (y * (N - 1) + 0.5) / N else y
}

realdata_specs <- list(
  list(name = "GasolineYield",
       build = function() {
         data("GasolineYield", package = "betareg", envir = environment())
         df <- get("GasolineYield")
         data.frame(y = squeeze01(df$yield), temp = df$temp)
       },
       fmean = y ~ temp),
  list(name = "ReadingSkills",
       build = function() {
         data("ReadingSkills", package = "betareg", envir = environment())
         df <- get("ReadingSkills")
         data.frame(y = squeeze01(df$accuracy),
                    dyslexia = df$dyslexia, iq = df$iq)
       },
       fmean = y ~ dyslexia + iq),
  list(name = "FoodExpenditure",
       build = function() {
         data("FoodExpenditure", package = "betareg", envir = environment())
         df <- get("FoodExpenditure")
         data.frame(y = squeeze01(df$food / df$income),
                    income = df$income, persons = df$persons)
       },
       fmean = y ~ income + persons)
)

for (spec in realdata_specs) {
  dnm <- spec$name
  d   <- spec$build()
  n   <- nrow(d)
  fm  <- spec$fmean                          # mean-only; constant dispersion
  reps <- 10L

  ## fastsimplexreg (1 thread)
  cf <- run_cell("realdata", LAB_FAST1, dnm, n, "logit", NA, NA,
                 function() fastsimplexreg(fm, data = d, link = "logit",
                                           n_threads = 1L),
                 extract_fast, reps)
  add_res(cf$rows)
  add_coef(coef_rows("realdata", dnm, n, "logit", LAB_FAST1, cf$ext))
  fast_fit <- tryCatch(fastsimplexreg(fm, data = d, link = "logit",
                                      n_threads = 1L), error = function(e) NULL)

  ## simplexreg
  sr_fit <- NULL
  if (HAS_SR) {
    cs <- run_cell("realdata", "simplexreg", dnm, n, "logit", NA, NA,
                   function() simplexreg::simplexreg(fm, data = d, link = "logit"),
                   extract_sr, reps)
    add_res(cs$rows)
    add_coef(coef_rows("realdata", dnm, n, "logit", "simplexreg", cs$ext))
    sr_fit <- tryCatch(simplexreg::simplexreg(fm, data = d, link = "logit"),
                       error = function(e) NULL)
  }

  ## betareg (BETA; speed / capability comparator; logLik not comparable)
  if (HAS_BR) {
    cb <- run_cell("realdata", "betareg", dnm, n, "logit", NA, NA,
                   function() betareg::betareg(fm, data = d, link = "logit"),
                   extract_br, reps)
    add_res(cb$rows)
    add_coef(coef_rows("realdata", dnm, n, "logit", "betareg", cb$ext))
  }

  ## cor(fitted mean) fast vs simplexreg (should be ~1)
  if (!is.null(fast_fit) && !is.null(sr_fit)) {
    mu_fast <- as.numeric(fast_fit$fitted.values)
    mu_sr   <- tryCatch(as.numeric(sr_fit$meanmu), error = function(e) NA_real_)
    realdata_cor[[dnm]] <- suppressWarnings(
      tryCatch(stats::cor(mu_fast, mu_sr), error = function(e) NA_real_))
  }
}

## =====================================================================
## SCENARIO 5 — mixed (fastsimplexregmixed only; no CRAN competitor)
##   scaling vs J clusters and nAGQ, at n_threads 1 and 4.
## =====================================================================
message("== Scenario: mixed ==")
sim_mixed <- function(J, nj = 8L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n  <- J * nj
  g  <- factor(rep(seq_len(J), each = nj))
  x1 <- stats::rnorm(n)
  z1 <- stats::rnorm(n)
  b  <- stats::rnorm(J, 0, 0.7)[g]
  mu <- fastsimplexreg::simplex_linkinv(0.3 - 0.6 * x1 + b, "logit")
  y  <- fastsimplexreg::rsimplex(n, mu, exp(-0.4 + 0.3 * z1))
  data.frame(y = y, x1 = x1, z1 = z1, g = g)
}

for (J in c(50L, 200L, 1000L)) {
  d <- sim_mixed(J, seed = 500L + J)
  n <- nrow(d)
  reps <- if (J <= 200L) 5L else 3L
  for (q in c(1L, 7L, 15L)) {
    ## 1 thread
    cm1 <- run_cell("mixed", "fastsimplexregmixed", "simulated", n, "logit", q, J,
                    function() fastsimplexregmixed(y ~ x1 | z1, random = ~ 1 | g,
                                                   data = d, nAGQ = q, n_threads = 1L),
                    extract_mixed, reps)
    add_res(cm1$rows)
    ## 4 threads
    cm4 <- run_cell("mixed", "fastsimplexregmixed (4 threads)", "simulated", n,
                    "logit", q, J,
                    function() fastsimplexregmixed(y ~ x1 | z1, random = ~ 1 | g,
                                                   data = d, nAGQ = q,
                                                   n_threads = N_THREADS_MULTI),
                    extract_mixed, reps)
    add_res(cm4$rows)
  }
}

## =====================================================================
## Assemble, coerce types, and save
## =====================================================================
benchmark_results <- do.call(rbind, RES_LIST)
rownames(benchmark_results) <- NULL
benchmark_results <- within(benchmark_results, {
  scenario  <- as.character(scenario)
  package   <- as.character(package)
  dataset   <- as.character(dataset)
  n         <- as.integer(n)
  link      <- as.character(link)
  nAGQ      <- as.integer(nAGQ)
  J         <- as.integer(J)
  rep       <- as.integer(rep)
  time_sec  <- as.numeric(time_sec)
  loglik    <- as.numeric(loglik)
  npar      <- as.integer(npar)
  converged <- as.logical(converged)
  note      <- as.character(note)
})
benchmark_results <- benchmark_results[, c(
  "scenario","package","dataset","n","link","nAGQ","J","rep",
  "time_sec","loglik","npar","converged","note")]

benchmark_coefs <- do.call(rbind, COEF_LIST[!vapply(COEF_LIST, is.null, logical(1))])
rownames(benchmark_coefs) <- NULL
benchmark_coefs <- within(benchmark_coefs, {
  scenario <- as.character(scenario)
  dataset  <- as.character(dataset)
  n        <- as.integer(n)
  link     <- as.character(link)
  package  <- as.character(package)
  term     <- as.character(term)
  estimate <- as.numeric(estimate)
})
benchmark_coefs <- benchmark_coefs[, c(
  "scenario","dataset","n","link","package","term","estimate")]

total_runtime_sec <- proc.time()[["elapsed"]] - t_start_all

benchmark_meta <- list(
  r_version        = R.version.string,
  platform         = R.version$platform,
  blas             = tryCatch(La_library(), error = function(e) NA_character_),
  blas_lapack      = tryCatch(extSoftVersion()[["BLAS"]], error = function(e) NA_character_),
  lapack           = tryCatch(La_version(), error = function(e) NA_character_),
  n_cores          = parallel::detectCores(),
  n_threads_multi  = N_THREADS_MULTI,
  budget_sec       = BUDGET_SEC,
  date             = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  seed             = 20260712,
  total_runtime_sec = total_runtime_sec,
  versions = list(
    fastsimplexreg = as.character(utils::packageVersion("fastsimplexreg")),
    simplexreg     = if (HAS_SR) as.character(utils::packageVersion("simplexreg")) else NA_character_,
    betareg        = if (HAS_BR) as.character(utils::packageVersion("betareg")) else NA_character_,
    microbenchmark = tryCatch(as.character(utils::packageVersion("microbenchmark")),
                              error = function(e) NA_character_)
  ),
  realdata_fitted_cor = realdata_cor,   # cor(fitted mu) fast vs simplexreg (~1)
  scenario_descriptions = list(
    accuracy = paste(
      "fastsimplexreg vs simplexreg on the SAME simplex model. Simulated logit",
      "mean (2 covariates) + log-dispersion (1 covariate); n in {500,2000,8000}.",
      "Compares coefficients and logLik (same family) and fit time. Note:",
      "simplexreg reports logLik under a different additive convention than",
      "fastsimplexreg/dsimplex, so compare COEFFICIENTS for accuracy."),
    scaling = paste(
      "Workhorse timing vs n in {200,1000,5000,20000,1e5,5e5} for",
      "fastsimplexreg(1 thread), fastsimplexreg(4 threads), simplexreg and",
      "betareg. Fixed mean(2)+dispersion(1) model. Per-fit budget of",
      paste0(BUDGET_SEC, "s"), "skips a package at larger n once exceeded.",
      "betareg fits the BETA model: SPEED comparator only."),
    links = paste(
      "fastsimplexreg vs simplexreg across all four mean links",
      "(logit/probit/cloglog/neglog) at n=3000; data simulated under the",
      "matching link. Reports time, logLik and convergence."),
    realdata = paste(
      "betareg real data sets GasolineYield (yield), ReadingSkills (accuracy),",
      "FoodExpenditure (food/income). Constant-dispersion 3-way fit",
      "(fastsimplexreg, simplexreg, betareg). Reports logLik, fit time,",
      "coefficients, and cor(fitted mu) fast-vs-simplexreg (in meta). betareg",
      "logLik not comparable across families. Boundary 0/1 values, if any, are",
      "squeezed via (y*(N-1)+0.5)/N identically for all packages."),
    mixed = paste(
      "fastsimplexregmixed only (no CRAN simplex/beta GLMM exists). Scaling of",
      "the mixed fit vs J clusters in {50,200,1000} and nAGQ in {1,7,15}, at",
      "n_threads 1 and 4, random intercept ~1|g with 8 obs/cluster.")
  )
)

saveRDS(benchmark_results, file.path(SCRATCH, "benchmark_results.rds"))
saveRDS(benchmark_coefs,   file.path(SCRATCH, "benchmark_coefs.rds"))
saveRDS(benchmark_meta,    file.path(SCRATCH, "benchmark_meta.rds"))

message(sprintf("Done. results=%d rows, coefs=%d rows, runtime=%.1f s",
                nrow(benchmark_results), nrow(benchmark_coefs), total_runtime_sec))
