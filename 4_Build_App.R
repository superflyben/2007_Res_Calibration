# Script to build the self-contained calibration application
#
# Writes app/Calibration_App.html: one static page carrying both surfaces, the
# deviation between them, and a conductivity-to-chloride converter. No server and
# no R at view time -- the surfaces are embedded as JSON and everything else is
# plotly from a CDN, so the page can be served from GitHub Pages as-is.
#
# The page template lives in app/app_template.html with a /*__DATA__*/ placeholder.
# Edit the template for anything about layout or behaviour; edit this script for
# anything about what data gets embedded.

# Housekeeping ------------------------------------------------------------

library(jsonlite)
library(interp)
library(tidyverse)

source("3_Compare_Surfaces.R")   # which sources scripts 1 and 2

# Display surfaces --------------------------------------------------------

# The grid of record is 20x20, which is right for differencing the two surfaces
# but coarse for display and for inverting conductivity back to salt. Re-grid
# each surface more finely over its own range, with the same fit each one already
# uses, purely for the page.
n_fine <- 80

# Span the measured points rather than the grid of record: the grid starts at
# T = 202 while the data starts at 202.1, and that sliver is outside the convex
# hull, so a grid-spanning surface has a dead edge where the converter silently
# returns nothing.
#
# Inset the range by a hair as well. The extreme nodes otherwise land exactly on
# the hull boundary -- the top salt edge runs straight between the 1500 mmol/kg
# measurements -- and interp drops boundary nodes, which cost 15 of the 80
# columns their S = 1500 value and made the converter refuse the corner
# calibration points themselves. 1e-9 of the range is enough to put them inside.
t_rng <- range(cal07_no_dup$Temperature_C)
s_rng <- range(cal07_no_dup$NaCl_mmol_kg)
t_eps <- diff(t_rng) * 1e-9
s_eps <- diff(s_rng) * 1e-9

fine07 <- interp(x = cal07_no_dup$Temperature_C,
                 y = cal07_no_dup$NaCl_mmol_kg,
                 z = cal07_no_dup$C,
                 xo = seq(t_rng[1] + t_eps, t_rng[2] - t_eps, length.out = n_fine),
                 yo = seq(s_rng[1] + s_eps, s_rng[2] - s_eps, length.out = n_fine),
                 extrap = FALSE)

# The 2000 surface is already dense from Matlab (200x200 over its own range);
# subsample rather than refit, so the page shows the surface Matlab produced.
i_sub <- round(seq(1, length(T_nat_00), length.out = n_fine))
j_sub <- round(seq(1, length(S_nat_00), length.out = n_fine))

# Serialisation helpers ---------------------------------------------------

# plotly wants z[row][col] at (x[col], y[row]), i.e. rows indexed by salt and
# columns by temperature -- the transpose of the z[T, S] layout the R scripts
# use. NaN has no JSON representation, so it goes out as null and the page's
# interpolation treats null as "no data here".
as_z <- function(z) {
    m <- t(z)                       # -> [salt, temperature]
    out <- lapply(seq_len(nrow(m)), function(i) {
        row <- m[i, ]
        lapply(row, function(v) if (is.finite(v)) v else NULL)
    })
    out
}

surface_payload <- function(tvec, svec, z) {
    list(t = tvec, s = svec, z = as_z(z))
}

# Band statistics ---------------------------------------------------------

# Recomputed here rather than passed along, so the page and the console report of
# script 3 cannot drift apart.
band_stats <- resid %>%
    mutate(z07 = z07[overlap], z00 = z00[overlap],
           band = cut(S, breaks = c(0, 200, 350, 500, Inf),
                      labels = c("<200", "200-350", "350-500", ">500"))) %>%
    group_by(band) %>%
    summarise(n = n(),
              f = sum(z00 * z07) / sum(z00^2),
              mean = mean(dev_pct),
              .groups = "drop") %>%
    mutate(rel = f / .env$f)

# Assemble ----------------------------------------------------------------

payload <- list(
    f = f,
    s07 = surface_payload(fine07$x, fine07$y, fine07$z),
    s00 = surface_payload(T_nat_00[i_sub], S_nat_00[j_sub], z00_nat[i_sub, j_sub]),
    dev = surface_payload(T_grid, S_grid, deviation_pct),
    p07 = list(t = cal07_no_dup$Temperature_C,
               s = cal07_no_dup$NaCl_mmol_kg,
               c = cal07_no_dup$C),
    p00 = list(t = cal00$Temperature_C,
               s = cal00$NaCl_mmol_kg,
               c = cal00$C),
    stats = list(
        n = sum(overlap),
        tmin = min(T_grid[idx[, 1]]), tmax = max(T_grid[idx[, 1]]),
        smin = min(S_grid[idx[, 2]]), smax = max(S_grid[idx[, 2]]),
        rms = sqrt(mean(deviation_pct[overlap]^2)),
        corT = cor(resid$T, resid$dev),
        corS = cor(resid$S, resid$dev),
        # the single worst node, which the plot makes prominent: it sits where
        # the 2007 surface is at its lowest over the overlap, so a middling
        # absolute error reads as a large percentage
        worstT = resid$T[which.max(abs(resid$dev_pct))],
        worstS = resid$S[which.max(abs(resid$dev_pct))],
        worstPct = max(abs(resid$dev_pct)),
        worstC = z07[overlap][which.max(abs(resid$dev_pct))],
        rmsExcl = sqrt(mean(resid$dev_pct[-which.max(abs(resid$dev_pct))]^2)),
        # a data frame serialises as an array of objects, which is what the
        # page's table loop expects -- splitting it into one-row frames instead
        # yields an array of arrays and the loop reads undefined
        bands = band_stats %>%
            transmute(band = as.character(band), n, f, rel, mean)
    )
)

json <- toJSON(payload, auto_unbox = TRUE, digits = 6, null = "null", na = "null")

template <- readLines("app/app_template.html", warn = FALSE) %>% paste(collapse = "\n")
stopifnot(grepl("/*__DATA__*/", template, fixed = TRUE))
page <- sub("/*__DATA__*/", json, template, fixed = TRUE)

writeLines(page, "app/Calibration_App.html", useBytes = TRUE)

message(sprintf("\nWrote app/Calibration_App.html (%.0f kB)",
                file.size("app/Calibration_App.html") / 1024))
message(sprintf("  2007 display surface %dx%d, %d finite",
                length(fine07$x), length(fine07$y), sum(is.finite(fine07$z))))
message(sprintf("  2000 display surface %dx%d, %d finite",
                n_fine, n_fine, sum(is.finite(z00_nat[i_sub, j_sub]))))
message(sprintf("  deviation surface %dx%d, %d finite",
                length(T_grid), length(S_grid), sum(overlap)))
