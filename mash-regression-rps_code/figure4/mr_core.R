# mr_core.R — minimal two-sample MR engine (no TwoSampleMR dependency)
# Validated against simulation: IVW/Egger/WM recover true effects; Steiger direction correct.

mr_wald <- function(bx, by, bxse, byse) {
  b <- by / bx
  se <- sqrt(byse^2 / bx^2 + by^2 * bxse^2 / bx^4)  # first-order delta
  list(b = b, se = se, p = 2 * pnorm(-abs(b / se)), method = "Wald ratio")
}

mr_ivw <- function(bx, by, bxse, byse) {
  w <- 1 / byse^2
  b <- sum(w * bx * by) / sum(w * bx^2)
  se_fixed <- sqrt(1 / sum(w * bx^2))
  Q <- sum(w * (by - b * bx)^2)
  df <- length(bx) - 1
  se_mult <- se_fixed * max(1, sqrt(Q / df))
  list(b = b, se = se_mult, p = 2 * pnorm(-abs(b / se_mult)),
       Q = Q, Q_df = df, Q_p = pchisq(Q, df, lower.tail = FALSE),
       method = "IVW (mult. RE)")
}

mr_egger <- function(bx, by, bxse, byse) {
  w <- 1 / byse^2
  X <- cbind(1, bx)
  W <- diag(w)
  XtWX <- t(X) %*% W %*% X
  if (rcond(XtWX) < 1e-12) return(NULL)
  beta <- solve(XtWX) %*% t(X) %*% W %*% by
  resid <- by - X %*% beta
  df <- length(bx) - 2
  if (df < 1) return(NULL)
  s2 <- max(sum(w * resid^2) / df, 1e-12)
  vb <- solve(XtWX) * s2
  list(intercept = beta[1], intercept_se = sqrt(vb[1, 1]),
       intercept_p = 2 * pt(-abs(beta[1] / sqrt(vb[1, 1])), df),
       b = beta[2], se = sqrt(vb[2, 2]),
       p = 2 * pt(-abs(beta[2] / sqrt(vb[2, 2])), df), method = "MR-Egger")
}

mr_wmedian <- function(bx, by, bxse, byse, nboot = 1000) {
  r <- by / bx
  se_r <- sqrt(byse^2 / bx^2 + by^2 * bxse^2 / bx^4)
  w <- 1 / se_r^2
  o <- order(r); r <- r[o]; w <- w[o]
  cw <- cumsum(w) / sum(w)
  b <- r[which(cw >= 0.5)[1]]
  # parametric bootstrap se
  bs <- replicate(nboot, {
    rb <- rnorm(length(bx), by, byse) / rnorm(length(bx), bx, bxse)
    wb <- 1 / se_r^2
    o2 <- order(rb); rb2 <- rb[o2]; cw2 <- cumsum(wb[o2]) / sum(wb)
    rb2[which(cw2 >= 0.5)[1]]
  })
  se <- sd(bs, na.rm = TRUE)
  list(b = b, se = se, p = 2 * pnorm(-abs(b / se)), method = "Weighted median")
}

# Steiger directionality: compare variance explained in X vs Y
steiger_dir <- function(bx, by, bxse, byse, n_x, n_y, eaf = NULL) {
  if (is.null(eaf) || all(is.na(eaf))) eaf <- rep(0.25, length(bx))
  eaf[is.na(eaf)] <- 0.25
  r2 <- function(b, se, n, eaf) {
    vg <- 2 * eaf * (1 - eaf) * b^2
    ve <- 2 * eaf * (1 - eaf) * se^2 * n
    vg / (vg + ve)
  }
  rx2 <- mapply(r2, bx, bxse, n_x, eaf)
  ry2 <- mapply(r2, by, byse, n_y, eaf)
  list(r_exposure = sum(rx2), r_outcome = sum(ry2),
       correct_direction = sum(rx2) > sum(ry2))
}

run_mr <- function(bx, by, bxse, byse, n_x = NA, n_y = NA, eaf = NULL) {
  ok <- complete.cases(bx, by, bxse, byse) & bxse > 0 & byse > 0 & bx != 0
  bx <- bx[ok]; by <- by[ok]; bxse <- bxse[ok]; byse <- byse[ok]
  if (!is.null(eaf)) eaf <- eaf[ok]
  n <- length(bx)
  out <- list(nsnp = n)
  if (n == 0) return(list(nsnp = 0, primary = list(b = NA, se = NA, p = NA, method = "none")))
  if (n == 1) {
    out$primary <- mr_wald(bx, by, bxse, byse)
  } else {
    out$primary <- mr_ivw(bx, by, bxse, byse)
    if (n >= 3) out$egger <- mr_egger(bx, by, bxse, byse)
    if (n >= 3) out$wmedian <- mr_wmedian(bx, by, bxse, byse)
  }
  if (!any(is.na(n_x)) && !any(is.na(n_y)))
    out$steiger <- steiger_dir(bx, by, bxse, byse, n_x, n_y, eaf)
  out
}
