# Summaries and figures from results/sim_*.rds and results/pbc3_splits.rds.
#
#   Rscript analysis/summarise.R
#
# Writes results/summary_sim_*.csv, results/summary_pbc3_*.csv,
# results/fig_coverage_sim.png and results/fig_pbc3.png.

suppressMessages(library(ggplot2))
pal <- c(native = "#2a78d6", hpd = "#eb6834", cdf = "#1baf7a", abs = "#eda100",
         ipcw = "#e87ba4", conservative = "#008300")
q <- function(x, p) unname(stats::quantile(x, p, na.rm = TRUE))
summ <- function(d, by, vars) {
  do.call(rbind, lapply(split(d, d[by], drop = TRUE), function(g) {
    out <- g[1, by, drop = FALSE]
    for (v in vars) {
      out[[paste0(v, "_mean")]] <- mean(g[[v]], na.rm = TRUE)
      out[[paste0(v, "_p05")]] <- q(g[[v]], 0.05)
      out[[paste0(v, "_p95")]] <- q(g[[v]], 0.95)
    }
    out$n <- nrow(g)
    out
  }))
}
theme_c <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        legend.position = "none", strip.background = element_blank(),
        strip.text = element_text(face = "bold"))

# ---------------- simulation ----------------
sim_files <- list.files("results", "^sim_.*\\.rds$", full.names = TRUE)
if (length(sim_files) > 0) {
  long <- event <- risk <- cond <- NULL
  for (f in sim_files) {
    x <- readRDS(f)
    for (r in x$results) {
      long <- rbind(long, cbind(scenario = x$scenario, rep = r$rep, r$long))
      event <- rbind(event, cbind(scenario = x$scenario, rep = r$rep, r$event))
      risk <- rbind(risk, cbind(scenario = x$scenario, rep = r$rep, r$risk))
      cond <- rbind(cond, cbind(scenario = x$scenario, rep = r$rep, r$long_cond))
    }
  }
  sl <- summ(long, c("scenario", "method"), c("coverage", "width"))
  se <- summ(event, c("scenario", "method"), c("coverage", "median_bound"))
  sr <- summ(risk, c("scenario", "method"), c("citl", "ece", "brier", "width"))
  sc <- summ(cond, c("scenario", "method", "tertile"), "coverage")
  write.csv(sl, "results/summary_sim_biomarker.csv", row.names = FALSE)
  write.csv(se, "results/summary_sim_eventtime.csv", row.names = FALSE)
  write.csv(sr, "results/summary_sim_risk.csv", row.names = FALSE)
  write.csv(sc, "results/summary_sim_conditional.csv", row.names = FALSE)
  print(sl); print(se); print(sr); print(sc)

  pd <- rbind(
    data.frame(target = "Biomarker at s + h", scenario = long$scenario,
               method = long$method, coverage = long$coverage),
    data.frame(target = "Residual event time (lower bound)", scenario = event$scenario,
               method = event$method, coverage = event$coverage))
  pd$scenario <- factor(pd$scenario, c("mild", "severe"),
                        c("Mild misspecification", "Severe misspecification"))
  pd$method <- factor(pd$method, c("native", "hpd", "cdf", "abs", "ipcw", "conservative"))
  ps <- summ(pd, c("target", "scenario", "method"), "coverage")
  g <- ggplot(ps, aes(method, coverage_mean, colour = method)) +
    geom_hline(yintercept = 0.9, linetype = "dashed", colour = "grey45") +
    geom_errorbar(aes(ymin = coverage_p05, ymax = coverage_p95), width = 0, linewidth = 0.8) +
    geom_point(size = 2.8) +
    facet_grid(scenario ~ target, scales = "free_x", space = "free_x") +
    scale_colour_manual(values = pal) +
    labs(x = NULL, y = "Coverage (mean, 5th-95th pct over replicates)",
         title = "Nominal 90%: BJM native vs. conformal") + theme_c
  ggsave("results/fig_coverage_sim.png", g, width = 8, height = 5.5, dpi = 200)
}

# ---------------- pbc3 ----------------
if (file.exists("results/pbc3_splits.rds")) {
  x <- readRDS("results/pbc3_splits.rds")
  long <- do.call(rbind, lapply(x$results, function(r) cbind(split = r$split, r$long)))
  event <- do.call(rbind, lapply(x$results, function(r) cbind(split = r$split, r$event)))
  risk <- do.call(rbind, lapply(x$results, function(r) cbind(split = r$split, r$risk)))
  risk$diff <- risk$mean_pred - risk$observed_km
  # pooled coverage over splits (each split's test patients, weighted by count)
  pooled <- do.call(rbind, lapply(split(long, long[c("biomarker", "method")], drop = TRUE),
    function(g) data.frame(biomarker = g$biomarker[1], method = g$method[1],
                           coverage_pooled = sum(g$coverage * g$n_test) / sum(g$n_test))))
  sl <- merge(summ(long, c("biomarker", "method"), c("coverage", "width")), pooled)
  se <- summ(event, "method", c("coverage_ipcw", "median_bound"))
  sr <- summ(risk, "method", c("mean_pred", "observed_km", "diff", "width"))
  write.csv(sl, "results/summary_pbc3_biomarker.csv", row.names = FALSE)
  write.csv(se, "results/summary_pbc3_eventtime.csv", row.names = FALSE)
  write.csv(sr, "results/summary_pbc3_risk.csv", row.names = FALSE)
  print(sl); print(se); print(sr)

  pd <- rbind(
    data.frame(panel = long$biomarker, method = long$method, coverage = long$coverage),
    data.frame(panel = "Event time (IPCW)", method = event$method, coverage = event$coverage_ipcw))
  pd$method <- factor(pd$method, c("native", "hpd", "cdf", "ipcw", "conservative"))
  pd$panel <- factor(pd$panel, c("serBilir", "albumin", "Event time (IPCW)"))
  ps <- summ(pd, c("panel", "method"), "coverage")
  g <- ggplot(ps, aes(method, coverage_mean, colour = method)) +
    geom_hline(yintercept = 0.9, linetype = "dashed", colour = "grey45") +
    geom_errorbar(aes(ymin = coverage_p05, ymax = coverage_p95), width = 0, linewidth = 0.8) +
    geom_point(size = 2.8) +
    facet_wrap(~ panel, scales = "free_x") +
    scale_colour_manual(values = pal) +
    labs(x = NULL, y = "Coverage (mean, 5th-95th pct)",
         title = sprintf("pbc3, landmark %g, horizon %g: %d random splits",
                         x$s, x$h, length(x$results))) + theme_c
  ggsave("results/fig_pbc3.png", g, width = 8, height = 3.8, dpi = 200)
}
