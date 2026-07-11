#' Diagnostic Plots for a Simplex Regression Fit
#'
#' @description
#' Produces model-diagnostic plots for a fitted `"simplex_fast"` object using
#' \pkg{ggplot2}. Up to four panels are available:
#'
#' \enumerate{
#'   \item Residuals against the fitted mean \eqn{\hat\mu}.
#'   \item Normal quantile-quantile plot of the residuals, with a reference line
#'     through the first and third quartiles.
#'   \item Scale-location plot (\eqn{\sqrt{|\text{residual}|}} against
#'     \eqn{\hat\mu}).
#'   \item Observed response against the fitted mean, with the identity line.
#' }
#'
#' @details
#' Each panel is a self-contained `ggplot` object. When more than one panel is
#' requested the panels are combined into a single figure with \pkg{patchwork}
#' if that package is installed; otherwise they are drawn one after another. The
#' (list of) `ggplot` object(s) is returned invisibly, so the plots can be
#' further customised or re-arranged by the caller.
#'
#' @param x An object of class `"simplex_fast"`.
#' @param which Integer vector selecting the panels to draw, a subset of
#'   `1:4`.
#' @param type Type of residual used in panels 1-3: one of `"deviance"`,
#'   `"pearson"` or `"response"`. See [residuals.simplex_fast()].
#' @param smooth Logical; if `TRUE` (the default) a LOESS smoother is added to
#'   the residual and scale-location panels when there are enough observations.
#' @param ... Additional arguments, currently ignored.
#'
#' @return Invisibly, a single `ggplot` object when one panel is requested, a
#'   combined \pkg{patchwork} object when several panels are requested and
#'   \pkg{patchwork} is available, or a named list of `ggplot` objects otherwise.
#'
#' @seealso [fastsimplexreg()], [residuals.simplex_fast()]
#'
#' @examples
#' set.seed(6)
#' n <- 400
#' dat <- data.frame(x1 = rnorm(n), z1 = rnorm(n))
#' mu <- simplex_linkinv(0.2 + 0.7 * dat$x1, link = "logit")
#' dat$y <- rsimplex(n, mu, exp(-0.5 + 0.4 * dat$z1))
#' fit <- fastsimplexreg(y ~ x1 | z1, data = dat, n_threads = 1L)
#' p <- plot(fit, which = 1:4)
#' # 'p' is a ggplot/patchwork object that can be further customised.
#'
#' @export
plot.simplex_fast <- function(x,
                              which = 1:4,
                              type = c("deviance", "pearson", "response"),
                              smooth = TRUE,
                              ...) {
  .simplex_diag_plot(x, which = which, type = match.arg(type), smooth = smooth)
}


# Internal shared ggplot2 diagnostic builder, used by the plot methods of both
# the fixed-effects and mixed-effects simplex fits. It relies only on the
# generic accessors residuals()/fitted() and on x$residuals, so it applies to
# any fitted object exposing them.
.simplex_diag_plot <- function(x, which = 1:4,
                               type = c("deviance", "pearson", "response"),
                               smooth = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for the diagnostic plots.", call. = FALSE)
  }
  type <- match.arg(type)
  which <- as.integer(which)
  if (anyNA(which) || !all(which %in% 1:4)) {
    stop("'which' must be a subset of 1:4.", call. = FALSE)
  }
  which <- sort(unique(which))

  mu <- x$fitted.values
  res <- residuals(x, type = type)
  y <- mu + x$residuals
  n <- length(res)
  res_label <- paste0(
    toupper(substring(type, 1L, 1L)), substring(type, 2L), " residuals"
  )

  df <- data.frame(
    fitted = mu,
    response = y,
    resid = res,
    root_abs = sqrt(abs(res))
  )

  add_smooth <- function(p) {
    if (isTRUE(smooth) && n >= 10L) {
      p <- p + ggplot2::geom_smooth(
        method = "loess", formula = y ~ x, se = FALSE,
        colour = "#2c7fb8", linewidth = 0.6
      )
    }
    p
  }

  base_theme <- ggplot2::theme_bw(base_size = 11)

  panels <- vector("list", 4L)

  # Panel 1: residuals vs fitted mean.
  panels[[1L]] <- add_smooth(
    ggplot2::ggplot(df, ggplot2::aes(x = .data$fitted, y = .data$resid)) +
      ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
      ggplot2::geom_point(alpha = 0.5)
  ) +
    ggplot2::labs(
      x = expression(hat(mu)), y = res_label, title = "Residuals vs fitted"
    ) +
    base_theme

  # Panel 2: normal Q-Q plot with a quartile reference line.
  qq <- data.frame(
    theoretical = stats::qnorm(stats::ppoints(n)),
    sample = sort(res)
  )
  qy <- stats::quantile(res, c(0.25, 0.75), names = FALSE, type = 7)
  qx <- stats::qnorm(c(0.25, 0.75))
  slope <- diff(qy) / diff(qx)
  intercept <- qy[1L] - slope * qx[1L]
  panels[[2L]] <-
    ggplot2::ggplot(qq, ggplot2::aes(x = .data$theoretical, y = .data$sample)) +
    ggplot2::geom_abline(
      slope = slope, intercept = intercept, colour = "grey40", linetype = 2
    ) +
    ggplot2::geom_point(alpha = 0.5) +
    ggplot2::labs(
      x = "Theoretical quantiles", y = "Sample quantiles", title = "Normal Q-Q"
    ) +
    base_theme

  # Panel 3: scale-location.
  panels[[3L]] <- add_smooth(
    ggplot2::ggplot(df, ggplot2::aes(x = .data$fitted, y = .data$root_abs)) +
      ggplot2::geom_point(alpha = 0.5)
  ) +
    ggplot2::labs(
      x = expression(hat(mu)),
      y = expression(sqrt(abs(residuals))),
      title = "Scale-Location"
    ) +
    base_theme

  # Panel 4: observed vs fitted mean.
  panels[[4L]] <-
    ggplot2::ggplot(df, ggplot2::aes(x = .data$fitted, y = .data$response)) +
    ggplot2::geom_abline(
      slope = 1, intercept = 0, colour = "grey40", linetype = 2
    ) +
    ggplot2::geom_point(alpha = 0.5) +
    ggplot2::labs(
      x = expression(hat(mu)), y = "Observed response",
      title = "Observed vs fitted"
    ) +
    base_theme

  selected <- panels[which]
  names(selected) <- as.character(which)

  if (length(selected) == 1L) {
    print(selected[[1L]])
    return(invisible(selected[[1L]]))
  }
  if (requireNamespace("patchwork", quietly = TRUE)) {
    combined <- patchwork::wrap_plots(
      selected, ncol = if (length(selected) > 1L) 2L else 1L
    )
    print(combined)
    return(invisible(combined))
  }
  for (p in selected) print(p)
  invisible(selected)
}
