# Script to read in the 2000 Juan de Fuca calibration data and build its surface
#
# This is the second of the two calibrations being compared. It is an independent
# experiment on a different sensor: 5 isohalines (constant salt, varying
# temperature) against the 2007 run's 6 isotherms, 20 points against 34, and a
# much narrower salt range. See data/2000_JdF/SOURCE.md for provenance.
#
# The surface is built here with the same algorithm and on the same grid as the
# 2007 surface in 1_Read_Calibration_Data.R, so that any difference between the
# two reflects the data rather than the fitting method. In particular this does
# NOT reuse the archived TSfit3 surface, which was built with a spline pre-pass
# and gridfit.

# Housekeeping ------------------------------------------------------------

library(interp)
library(R.matlab)
library(tidyverse)

# Read in the data --------------------------------------------------------

# The five isohaline files, unlabelled, one per NaCl concentration.
# Columns: NaCl (mol/kg), Temperature (deg C), Conductivity (1/V)
iso_files <- list.files("data/2000_JdF", pattern = "^mm[0-9]+\\.csv$",
                        full.names = TRUE)

cal00 <- iso_files %>%
    set_names(basename(.)) %>%
    map_dfr(read_csv,
            col_names = c("NaCl_mol_kg", "Temperature_C", "C"),
            col_types = "ddd",
            .id = "source_file")

# Salt is stored in mol/kg here and mmol/kg in the 2007 data. Convert once, at
# the point of ingest, so everything downstream is in mmol/kg.
cal00 <- cal00 %>%
    mutate(NaCl_mmol_kg = NaCl_mol_kg * 1000) %>%
    arrange(NaCl_mmol_kg, Temperature_C)

message(sprintf("2000 JdF: %d points, %d isohalines at %s mmol/kg",
                nrow(cal00),
                n_distinct(cal00$NaCl_mmol_kg),
                paste(sort(unique(cal00$NaCl_mmol_kg)), collapse = ", ")))
message(sprintf("  T %.1f-%.1f C, C %.3f-%.3f",
                min(cal00$Temperature_C), max(cal00$Temperature_C),
                min(cal00$C), max(cal00$C)))

# Output grid -------------------------------------------------------------

# The grid of record, defined in 1_Read_Calibration_Data.R: the 20x20 nodes the
# Matlab surface uses, read from surfdat.mat rather than restated, so that the
# two surfaces are on identical nodes by construction and can be differenced.
surfdat <- readMat("../Calibration/2007/surfdat.mat")
T_grid <- sort(unique(as.vector(surfdat$Tdat)))
S_grid <- sort(unique(as.vector(surfdat$Sdat)))

# Interpolation -----------------------------------------------------------

# Same call as the 2007 surface: linear, no extrapolation. Nodes outside the
# convex hull of the 20 data points come back NA, which is what restricts this
# surface to the overlap region.
R_surfdat_00 <- interp(x = cal00$Temperature_C,
                       y = cal00$NaCl_mmol_kg,
                       z = cal00$C,
                       xo = T_grid,
                       yo = S_grid,
                       extrap = FALSE)

# How much of the grid of record does the 2000 surface actually cover?
# The grid spans T 202-395.5 and S 100-1500, while this experiment only reaches
# T 104-382 and S 54-634, so most of it is necessarily empty.
n_finite <- sum(is.finite(R_surfdat_00$z))
message(sprintf("2000 surface on the grid of record: %d of %d nodes finite (%.0f%%)",
                n_finite, length(R_surfdat_00$z),
                100 * n_finite / length(R_surfdat_00$z)))

covered <- which(is.finite(R_surfdat_00$z), arr.ind = TRUE)
if (nrow(covered) > 0) {
    message(sprintf("  covering T %.1f-%.1f C, S %.0f-%.0f mmol/kg",
                    min(T_grid[covered[, 1]]), max(T_grid[covered[, 1]]),
                    min(S_grid[covered[, 2]]), max(S_grid[covered[, 2]])))
}

# Matrix format for plotting, matching 1_Read_Calibration_Data.R:
# interp returns z[i, j] = f(x[i], y[j]), temperature down rows, salt across
# columns -- the transpose of the Matlab layout.
R_surfdat_00$T_grid <- matrix(data = R_surfdat_00$x,
                              nrow = length(R_surfdat_00$x),
                              ncol = length(R_surfdat_00$y),
                              byrow = FALSE)
R_surfdat_00$S_grid <- matrix(data = R_surfdat_00$y,
                              nrow = length(R_surfdat_00$x),
                              ncol = length(R_surfdat_00$y),
                              byrow = TRUE)

# Caveats -----------------------------------------------------------------

# Two things to keep in mind before reading anything into this surface:
#
# 1. Four temperature points per isohaline, and C(T) is strongly non-monotonic --
#    it rises to a maximum near 300 C and falls away by 382 C. Linear
#    triangulation across 20 points represents that curvature far more crudely
#    than TSfit3's spline pre-pass did. This is the price of holding the
#    algorithm constant across the two datasets, and it is paid on the 2000 side
#    because that dataset is the sparser of the two.
# 2. Pressure is not constant: 220 bar except the 382 C block at 300 bar, where
#    the 2007 run sat near 256 kgf/cm2 throughout. The two calibrations are
#    therefore not at matched pressure at the top of the temperature range.
