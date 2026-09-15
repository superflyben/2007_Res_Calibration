# Script to compare the 2000 and 2007 calibration surfaces
#
# The question: the two calibrations were made on different sensors, and
# conductivity depends partly on electrode surface area, which varies between
# sensors through fabrication (GCA section 3.3.4). If that geometry difference is
# the whole story, then one multiplicative constant should reconcile the two
# surfaces everywhere, with no dependence on temperature or salt. That is exactly
# what licenses applying a single scale factor per deployment.
#
# So: fit ONE scale factor over the entire overlapping region and look at what is
# left. A residual that is flat and small supports the single-constant model. A
# residual with structure in T or S does not.
#
# The scale factor is fitted here only as a diagnostic. The factor the
# application applies to field data is a separate, user-supplied quantity -- see
# the note at the end of this script.

# Housekeeping ------------------------------------------------------------

library(plotly)
library(htmlwidgets)
library(tidyverse)

source("1_Read_Calibration_Data.R")        # 2007 surface, R interp::interp
source("2_Read_2000_Calibration_Data.R")   # 2000 surface, Matlab griddata

# Both surfaces are on the grid of record, laid out z[i, j] = f(T[i], S[j]).
stopifnot(all.equal(T_grid, T_grid_00), all.equal(S_grid, S_grid_00))
z07 <- R_surfdat$z

# Overlap -----------------------------------------------------------------

# Each surface is NA outside the convex hull of its own measurements, so the
# overlap falls out of the two NA masks without any extra masking. Nothing in
# what follows is fitted against extrapolated surface.
overlap <- is.finite(z07) & is.finite(z00)

idx <- which(overlap, arr.ind = TRUE)
message(sprintf("Overlap: %d of %d nodes (2007 has %d, 2000 has %d)",
                sum(overlap), length(overlap), sum(is.finite(z07)), sum(is.finite(z00))))
message(sprintf("  T %.1f-%.1f C, S %.0f-%.0f mmol/kg",
                min(T_grid[idx[, 1]]), max(T_grid[idx[, 1]]),
                min(S_grid[idx[, 2]]), max(S_grid[idx[, 2]])))

# Scale factor ------------------------------------------------------------

# 2007 is the reference: f is the factor that puts a 2000 conductivity into the
# 2007 frame, the same sense as Eq. 2 in the GCA paper, where a field reading is
# multiplied by f before the surface is inverted. Least squares over every
# overlapping node, one factor for the whole region.
f <- sum(z00[overlap] * z07[overlap]) / sum(z00[overlap]^2)

message(sprintf("\nSingle scale factor over the whole overlap:"))
message(sprintf("  f (2000 -> 2007) = %.4f      1/f = %.4f", f, 1 / f))

# For comparison, the node-by-node ratio, which f is a least-squares summary of
ratio <- z07[overlap] / z00[overlap]
message(sprintf("  node-by-node C_2007/C_2000: mean %.4f, median %.4f, sd %.4f (CV %.1f%%)",
                mean(ratio), median(ratio), sd(ratio), 100 * sd(ratio) / mean(ratio)))
message(sprintf("  range %.4f to %.4f", min(ratio), max(ratio)))

# Deviation surface -------------------------------------------------------

# What one constant fails to reconcile. Both absolute (in conductivity) and
# relative (as a percentage of the 2007 surface).
deviation <- f * z00 - z07
deviation[!overlap] <- NA
deviation_pct <- 100 * deviation / z07
deviation_pct[!overlap] <- NA

d <- deviation[overlap]
message(sprintf("\nAfter scaling, what is left:"))
message(sprintf("  residual RMS %.4f, max |residual| %.4f (in conductivity)",
                sqrt(mean(d^2)), max(abs(d))))
message(sprintf("  relative RMS %.1f%%, max %.1f%%",
                sqrt(mean(deviation_pct[overlap]^2)), max(abs(deviation_pct[overlap]))))

# Is the residual structured? ---------------------------------------------

# The single-constant model predicts a residual with no dependence on T or S.
# A one-directional trend in T is specifically expected if the 2000 surface is
# under-resolved through the 300-382 C plateau, so test for it rather than
# assuming the scatter is noise.
resid <- tibble(T = T_grid[idx[, 1]],
                S = S_grid[idx[, 2]],
                dev = d,
                dev_pct = deviation_pct[overlap])

message(sprintf("\nStructure in the residual:"))
message(sprintf("  correlation with T: %+.3f      with S: %+.3f",
                cor(resid$T, resid$dev), cor(resid$S, resid$dev)))

by_T <- resid %>%
    mutate(band = cut(T, breaks = c(0, 250, 300, 350, Inf),
                      labels = c("<250", "250-300", "300-350", ">350"))) %>%
    group_by(band) %>%
    summarise(n = n(), mean_pct = mean(dev_pct), rms_pct = sqrt(mean(dev_pct^2)),
              .groups = "drop")
message("  by temperature band (mean and RMS of the relative deviation):")
for (i in seq_len(nrow(by_T))) {
    message(sprintf("    %-8s n=%3d   mean %+6.1f%%   RMS %5.1f%%",
                    by_T$band[i], by_T$n[i], by_T$mean_pct[i], by_T$rms_pct[i]))
}

by_S <- resid %>%
    mutate(band = cut(S, breaks = c(0, 200, 350, 500, Inf),
                      labels = c("<200", "200-350", "350-500", ">500"))) %>%
    group_by(band) %>%
    summarise(n = n(), mean_pct = mean(dev_pct), rms_pct = sqrt(mean(dev_pct^2)),
              .groups = "drop")
message("  by salt band (mean and RMS of the relative deviation):")
for (i in seq_len(nrow(by_S))) {
    message(sprintf("    %-8s n=%3d   mean %+6.1f%%   RMS %5.1f%%",
                    by_S$band[i], by_S$n[i], by_S$mean_pct[i], by_S$rms_pct[i]))
}

# Plot --------------------------------------------------------------------

p_dev <- plot_ly(x = matrix(T_grid, nrow = length(T_grid), ncol = length(S_grid)),
                 y = matrix(S_grid, nrow = length(T_grid), ncol = length(S_grid),
                            byrow = TRUE),
                 z = deviation_pct,
                 type = "surface",
                 colorscale = "RdBu", cmid = 0,
                 colorbar = list(title = "deviation (%)")) %>%
    layout(title = sprintf("2000 vs 2007 after one scale factor (f = %.3f)", f),
           scene = list(xaxis = list(title = "Temperature (C)"),
                        yaxis = list(title = "NaCl (mmol/kg)"),
                        zaxis = list(title = "deviation (%)")))

saveWidget(as_widget(p_dev), "Deviation_Surface.html")

# Note on the applied scale factor ----------------------------------------

# f above is a diagnostic: it measures how well one constant reconciles the two
# calibrations. It is NOT the factor to apply to a field deployment. That one
# comes from pairing field readings with discrete fluid samples of known
# chloride (GCA Eq. 2, f = C_lab / C_field), which is a per-deployment quantity
# and needs sample data this comparison does not have. The application takes it
# as an input, defaulting to 1.
