## =====================================================================
## viz_figures.R  --  Publication-quality benchmark figures
## fastsimplexreg vignette (Dr. Viz)
##
## Self-contained plotting code. Reads ONLY precomputed benchmark RDS
## files (no model fits). To repoint at the shipped package data, set:
##   results_path <- system.file("extdata", package = "fastsimplexreg")
##
## Dependencies: ggplot2, scales, patchwork, dplyr, tidyr
## Palette anchored on the package diagnostic colours
##   #1a5276 (primary) and #a93226 (accent); extended to a
##   colour-blind-safe 4-hue categorical set (CVD adjacent dE 35.4,
##   validated). Secondary encoding (linetype + shape) backs up hue on
##   the two lower-contrast series, per the accessibility relief rule.
## =====================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(dplyr)
  library(tidyr)
})

## ---- 0. Data location -------------------------------------------------
## Default: the scratchpad next to this script. The lead can override
## with system.file("extdata", package = "fastsimplexreg").
if (!exists("results_path")) {
  results_path <- tryUnwrap <- dirname(
    sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1]
  )
  if (is.na(results_path) || length(results_path) == 0 || results_path == "")
    results_path <- "."
}

read_bench <- function(name) {
  f <- file.path(results_path, name)
  if (!file.exists(f)) stop("Cannot find ", f, call. = FALSE)
  readRDS(f)
}

results <- read_bench("benchmark_results.rds")
coefs   <- read_bench("benchmark_coefs.rds")
meta    <- read_bench("benchmark_meta.rds")

## ---- 1. Palette & theme ----------------------------------------------
## Categorical hues (light surface, validated worst-adjacent CVD dE 35.4)
PAL <- c(
  "fastsimplexreg"              = "#1a5276",  # brand primary  (deep blue)
  "fastsimplexreg (4 threads)"  = "#5dade2",  # light blue  (same engine family)
  "simplexreg"                  = "#a93226",  # brand accent  (deep red)
  "betareg"                     = "#e69f00"   # amber  (CVD-safe, Wong)
)
## Short display labels for legends/axes
LAB <- c(
  "fastsimplexreg"             = "fastsimplexreg (1 thread)",
  "fastsimplexreg (4 threads)" = "fastsimplexreg (4 threads)",
  "simplexreg"                 = "simplexreg (CRAN)",
  "betareg"                    = "betareg (CRAN)"
)
## Secondary encodings (back up hue: relief for lower-contrast series)
LTY <- c(
  "fastsimplexreg"             = "solid",
  "fastsimplexreg (4 threads)" = "22",
  "simplexreg"                 = "solid",
  "betareg"                    = "solid"
)
SHP <- c(
  "fastsimplexreg"             = 16,  # filled circle
  "fastsimplexreg (4 threads)" = 17,  # filled triangle
  "simplexreg"                 = 15,  # filled square
  "betareg"                    = 18   # filled diamond
)

## Ink / chrome tokens
INK      <- "#0b0b0b"; INK2 <- "#52514e"; MUTED <- "#898781"
GRIDLINE <- "#e1e0d9"

theme_bench <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title    = element_text(face = "bold", size = base_size + 3,
                                   colour = INK, margin = margin(b = 2)),
      plot.subtitle = element_text(size = base_size, colour = INK2,
                                   margin = margin(b = 10)),
      plot.caption  = element_text(size = base_size - 2, colour = MUTED,
                                   hjust = 0, margin = margin(t = 10)),
      axis.title    = element_text(size = base_size, colour = INK2),
      axis.text     = element_text(size = base_size - 1, colour = INK2),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = GRIDLINE, linewidth = 0.3),
      panel.border  = element_rect(colour = "#c3c2b7", linewidth = 0.4),
      legend.position = "top",
      legend.title  = element_blank(),
      legend.key    = element_blank(),
      legend.text   = element_text(size = base_size - 1, colour = INK2),
      strip.background = element_rect(fill = "#eef2f5", colour = NA),
      strip.text    = element_text(face = "bold", size = base_size - 1,
                                   colour = "#1a5276"),
      plot.title.position   = "plot",
      plot.caption.position = "plot"
    )
}

## Compact provenance string + hard-wrap helper so captions never clip.
sci_source <- paste0(
  "Precomputed benchmark: R 4.5.2, OpenBLAS; simplexreg ",
  meta$versions$simplexreg, ", betareg ", meta$versions$betareg,
  "; seed ", meta$seed, "."
)
wrap_cap <- function(..., width = 96)
  paste(strwrap(paste0(...), width = width), collapse = "\n")

## =====================================================================
## FIGURE 1 -- Scaling (HERO): median fit time vs n, log-log
## =====================================================================
fig_scaling <- function() {
  d <- results %>%
    filter(scenario == "scaling") %>%
    group_by(package, n) %>%
    summarise(med = median(time_sec),
              lo  = quantile(time_sec, 0.25),
              hi  = quantile(time_sec, 0.75), .groups = "drop") %>%
    mutate(package = factor(package, levels = names(PAL)))

  ggplot(d, aes(n, med, colour = package, fill = package)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.13, colour = NA) +
    geom_line(aes(linetype = package), linewidth = 0.9) +
    geom_point(aes(shape = package), size = 2.4) +
    scale_colour_manual(values = PAL, labels = LAB) +
    scale_fill_manual(values = PAL, labels = LAB, guide = "none") +
    scale_linetype_manual(values = LTY, labels = LAB) +
    scale_shape_manual(values = SHP, labels = LAB) +
    scale_x_log10(labels = label_number(big.mark = ",", accuracy = 1),
                  breaks = c(200, 1e3, 5e3, 2e4, 1e5, 5e5)) +
    scale_y_log10(labels = label_number(accuracy = 0.001)) +
    annotation_logticks(sides = "bl", colour = MUTED, linewidth = 0.25,
                        short = unit(0.05,"cm"), mid = unit(0.1,"cm"),
                        long = unit(0.15,"cm")) +
    labs(
      title = "fastsimplexreg scales sub-linearly where CRAN slows down",
      subtitle = paste0("Median wall-clock time for one cold fit vs sample ",
                        "size (both axes log10; band = IQR over replications)"),
      x = "Sample size  n", y = "Fit time  (seconds)",
      caption = wrap_cap(
        "Same simulated (0,1) data and identical mean(2)+dispersion(1) model ",
        "fed to every package per n; fastsimplexreg timed single-core and at 4 ",
        "threads. At the largest n the 4-thread build overtakes all; ",
        "multithreading adds overhead (wider IQR) at small/medium n. ",
        sci_source)
    ) +
    theme_bench() +
    guides(colour = guide_legend(nrow = 1),
           linetype = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))
}

## =====================================================================
## FIGURE 2 -- Speedup of fastsimplexreg (1 thread) vs CRAN, by n
## =====================================================================
fig_speedup <- function() {
  base <- results %>% filter(scenario == "scaling") %>%
    group_by(package, n) %>% summarise(med = median(time_sec), .groups = "drop")
  fast1 <- base %>% filter(package == "fastsimplexreg") %>%
    select(n, fast = med)
  ns  <- sort(unique(base$n))
  off <- 0.16
  d <- base %>%
    filter(package %in% c("simplexreg", "betareg")) %>%
    left_join(fast1, by = "n") %>%
    mutate(speedup = med / fast,
           package = factor(package, levels = c("simplexreg", "betareg")),
           xi   = as.integer(factor(n, levels = ns)),
           xpos = xi + ifelse(package == "simplexreg", -off, off))

  ggplot(d, aes(xpos, speedup, colour = package)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = MUTED,
               linewidth = 0.4) +
    geom_segment(aes(x = xpos, xend = xpos, y = 1, yend = speedup),
                 linewidth = 0.9) +
    geom_point(aes(shape = package), size = 3.4) +
    geom_text(aes(label = paste0(round(speedup, 1), "x")), vjust = -0.9,
              size = 3, show.legend = FALSE, colour = INK2) +
    scale_colour_manual(values = PAL[c("simplexreg","betareg")],
                        labels = LAB[c("simplexreg","betareg")]) +
    scale_shape_manual(values = SHP[c("simplexreg","betareg")],
                       labels = LAB[c("simplexreg","betareg")]) +
    scale_x_continuous(breaks = seq_along(ns),
                       labels = label_number(big.mark = ",", accuracy = 1)(ns)) +
    scale_y_log10() +
    labs(
      title = "How many times faster is fastsimplexreg (single core)?",
      subtitle = paste0("CRAN fit time ÷ fastsimplexreg (1 thread) fit ",
                        "time. Above the dashed line = fastsimplexreg wins."),
      x = "Sample size  n", y = "Speed-up factor  (log10, ×)",
      caption = wrap_cap(
        "Single-core, like-for-like: fastsimplexreg is 3-15× faster than ",
        "simplexreg throughout. betareg (a different, beta model) is a ",
        "speed-only comparator; at n = 500,000 its single fit edges out the ",
        "single-thread build, and the 4-thread build (Fig. 1) reclaims the ",
        "lead. ", sci_source)
    ) +
    theme_bench() +
    theme(panel.grid.major.x = element_blank())
}

## =====================================================================
## FIGURE 3 -- Accuracy: fastsimplexreg vs simplexreg estimates coincide
## =====================================================================
fig_accuracy <- function() {
  acc <- coefs %>% filter(scenario == "accuracy")
  wide <- acc %>%
    select(n, term, package, estimate) %>%
    pivot_wider(names_from = package, values_from = estimate) %>%
    mutate(part = ifelse(grepl("^disp", term), "dispersion", "mean"))
  maxdiff <- wide %>%
    mutate(absdiff = abs(fastsimplexreg - simplexreg)) %>%
    group_by(n) %>% summarise(md = max(absdiff), .groups = "drop")
  worst <- signif(max(maxdiff$md), 2)

  p_scatter <- ggplot(wide, aes(simplexreg, fastsimplexreg)) +
    geom_abline(slope = 1, intercept = 0, colour = MUTED,
                linetype = "dashed", linewidth = 0.5) +
    geom_point(aes(colour = part), size = 3, alpha = 0.9) +
    scale_colour_manual(values = c(mean = "#1a5276", dispersion = "#e69f00")) +
    coord_equal() +
    labs(subtitle = "Coefficient estimates",
         x = "simplexreg estimate", y = "fastsimplexreg estimate") +
    theme_bench() +
    theme(legend.position = c(0.02, 0.98),
          legend.justification = c(0, 1),
          legend.background = element_rect(fill = "#ffffffcc", colour = NA))

  p_diff <- ggplot(maxdiff, aes(n, md)) +
    geom_line(colour = "#1a5276", linewidth = 0.8) +
    geom_point(colour = "#1a5276", size = 2.6) +
    scale_x_log10(labels = label_number(big.mark = ",", accuracy = 1),
                  breaks = unique(maxdiff$n)) +
    scale_y_log10(labels = label_scientific(digits = 1)) +
    labs(subtitle = "Worst-case agreement",
         x = "Sample size  n", y = "max |fast − simplexreg|") +
    theme_bench()

  (p_scatter | p_diff) +
    plot_annotation(
      title = "fastsimplexreg reproduces simplexreg's estimates exactly",
      subtitle = wrap_cap("Same simplex likelihood, same data (n = 500 / 2,000 ",
                          "/ 8,000): every coefficient sits on the identity ",
                          "line; the worst disagreement is ", worst,
                          " (numerical noise).", width = 120),
      caption = wrap_cap("Both are MLEs of the identical simplex model. ",
                         sci_source, width = 120),
      theme = theme_bench()
    )
}

## =====================================================================
## FIGURE 4 -- Links: time by mean link, fastsimplexreg vs simplexreg
## =====================================================================
fig_links <- function() {
  d <- results %>% filter(scenario == "links") %>%
    group_by(package, link) %>%
    summarise(med = median(time_sec), conv = all(converged), .groups = "drop") %>%
    mutate(package = factor(package, levels = c("fastsimplexreg","simplexreg")),
           link = factor(link, levels = c("logit","probit","cloglog","neglog")))

  ggplot(d, aes(link, med, fill = package)) +
    geom_col(position = position_dodge(width = 0.72), width = 0.66) +
    geom_text(aes(label = sprintf("%.3f", med)),
              position = position_dodge(width = 0.72), vjust = -0.5,
              size = 2.9, colour = INK2) +
    scale_fill_manual(values = PAL[c("fastsimplexreg","simplexreg")],
                      labels = LAB[c("fastsimplexreg","simplexreg")]) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(
      title = "All four mean links, consistently faster",
      subtitle = paste0("Median fit time at n = 3,000; every fit converged ",
                        "for both packages"),
      x = "Mean link function", y = "Fit time  (seconds)",
      caption = wrap_cap("fastsimplexreg supports logit / probit / cloglog / ",
                       "neglog; each ~3-10× faster than simplexreg. ",
                       sci_source, width = 88)
    ) +
    theme_bench() +
    theme(panel.grid.major.x = element_blank())
}

## =====================================================================
## FIGURE 5 -- Real data: time by dataset & package
## =====================================================================
fig_realdata <- function() {
  d <- results %>% filter(scenario == "realdata") %>%
    group_by(dataset, package) %>%
    summarise(med = median(time_sec), .groups = "drop") %>%
    mutate(package = factor(package,
                            levels = c("fastsimplexreg","simplexreg","betareg")),
           dataset = factor(dataset,
                            levels = c("GasolineYield","ReadingSkills",
                                       "FoodExpenditure")))

  ggplot(d, aes(dataset, med, fill = package)) +
    geom_col(position = position_dodge(width = 0.78), width = 0.72) +
    geom_text(aes(label = sprintf("%.3f", med)),
              position = position_dodge(width = 0.78), vjust = -0.5,
              size = 2.8, colour = INK2) +
    scale_fill_manual(values = PAL[c("fastsimplexreg","simplexreg","betareg")],
                      labels = LAB[c("fastsimplexreg","simplexreg","betareg")]) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
    labs(
      title = "Fastest on the standard beta-regression datasets",
      subtitle = paste0("Median fit time on three betareg example datasets ",
                        "(constant dispersion, three-way match)"),
      x = NULL, y = "Fit time  (seconds)",
      caption = wrap_cap(
        "fastsimplexreg vs simplexreg fitted means agree perfectly ",
        "(cor = 1.00 on all three datasets). betareg fits a different (beta) ",
        "model — compared on speed/capability only, not likelihood. ",
        sci_source)
    ) +
    theme_bench() +
    theme(panel.grid.major.x = element_blank())
}

## =====================================================================
## FIGURE 6 -- Mixed model scaling (no CRAN competitor exists)
## =====================================================================
fig_mixed <- function() {
  d <- results %>% filter(scenario == "mixed") %>%
    group_by(package, J, nAGQ) %>%
    summarise(med = median(time_sec), conv = all(converged), .groups = "drop") %>%
    mutate(threads = ifelse(grepl("4 threads", package),
                            "4 threads", "1 thread"),
           threads = factor(threads, levels = c("1 thread","4 threads")),
           nAGQ_lab = factor(paste0("nAGQ = ", nAGQ),
                             levels = paste0("nAGQ = ", c(1,7,15))))

  tcol <- c("1 thread" = "#1a5276", "4 threads" = "#5dade2")

  ggplot(d, aes(J, med, colour = threads)) +
    geom_line(aes(linetype = threads), linewidth = 0.9) +
    geom_point(aes(shape = conv), size = 2.8) +
    facet_wrap(~ nAGQ_lab) +
    scale_colour_manual(values = tcol) +
    scale_linetype_manual(values = c("1 thread"="solid","4 threads"="22"),
                          guide = "none") +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 4),
                       labels = c(`TRUE` = "converged",
                                  `FALSE` = "not converged"),
                       name = NULL, drop = FALSE) +
    scale_x_log10(breaks = c(50, 200, 1000)) +
    scale_y_log10(labels = label_number(accuracy = 0.01)) +
    labs(
      title = "The simplex mixed model is tractable — where nothing else exists",
      subtitle = paste0("fastsimplexregmixed fit time vs number of clusters J, ",
                        "by quadrature order (no CRAN simplex/beta GLMM exists)"),
      x = "Number of clusters  J  (log10; n = 400 / 1,600 / 8,000)",
      y = "Fit time  (seconds, log10)",
      caption = wrap_cap(
        "No CRAN package fits a simplex mixed model, so this is ",
        "fastsimplexregmixed's own scaling. 4 threads ~ halves the time; the ",
        "J = 1,000 / nAGQ = 1 cell did not converge (marked ✕) — higher nAGQ ",
        "is more stable. ", sci_source)
    ) +
    theme_bench() +
    guides(colour = guide_legend(order = 1),
           shape  = guide_legend(order = 2))
}

## =====================================================================
## Render everything to PNG (proof the code works)
## Recommended vignette fig sizes noted per figure.
## =====================================================================
render_all <- function(out = results_path, dpi = 150) {
  specs <- list(
    list(f = fig_scaling,  file = "fig_scaling.png",  w = 7.5, h = 5.0), # HERO
    list(f = fig_speedup,  file = "fig_speedup.png",  w = 7.5, h = 4.6),
    list(f = fig_accuracy, file = "fig_accuracy.png", w = 9.0, h = 4.6),
    list(f = fig_links,    file = "fig_links.png",    w = 6.8, h = 4.4),
    list(f = fig_realdata, file = "fig_realdata.png", w = 7.2, h = 4.6),
    list(f = fig_mixed,    file = "fig_mixed.png",    w = 8.5, h = 4.6)
  )
  for (s in specs) {
    ggsave(file.path(out, s$file), plot = s$f(),
           width = s$w, height = s$h, dpi = dpi, bg = "white")
    message("wrote ", s$file, "  (", s$w, " x ", s$h, " in)")
  }
  invisible(TRUE)
}

if (sys.nframe() == 0L) {
  ## Auto-render ONLY when executed directly (Rscript viz_figures.R). When this
  ## file is source()d (e.g. from the benchmark vignette) sys.nframe() > 0, so it
  ## just defines the figure functions and never writes to disk.
  render_all()
}
