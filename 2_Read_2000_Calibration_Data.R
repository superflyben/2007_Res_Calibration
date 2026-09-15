# Script to read in the 2000 Juan de Fuca calibration data and its surface
#
# This is the second of the two calibrations being compared. It is an independent
# experiment on a different sensor: 5 isohalines (constant salt, varying
# temperature) against the 2007 run's 6 isotherms, 20 points against 34, and a
# much narrower salt range. See data/2000_JdF/SOURCE.md for provenance.
#
# The two surfaces are deliberately NOT built with the same algorithm. Each is
# fitted the way its own data warrants:
#
#   2007 -> R, interp::interp on the measured points (1_Read_Calibration_Data.R)
#   2000 -> Matlab, matlab/Build_2000_Surface.m
#
# The 2000 data has four measured temperatures per isohaline spanning a C(T)
# curve that plateaus near 300 C and falls away sharply above ~350 C. Forcing it
# through plain interpolation on the raw points does not just add noise, it adds
# a one-directional error: straight-lining the 300->382 gap underestimates C by
# 12-28% through the plateau (measured by removing the 2007 350 C point and
# reconstructing it the same way). The original fit handled this with a spline
# pre-pass along temperature, and that machinery is reused rather than replaced.
#
# So this script reads the Matlab surface rather than building one. Regenerate it
# with, from the matlab/ directory:
#
#   matlab -batch "Build_2000_Surface;"

# Housekeeping ------------------------------------------------------------

library(R.matlab)
library(tidyverse)

# Read in the measured points ---------------------------------------------

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

# Read in the surface -----------------------------------------------------

# Built by matlab/Build_2000_Surface.m: interpmmlarson spline along temperature
# within each isohaline, then griddata (cubic) onto the grid of record. griddata
# returns NA outside the convex hull of the data, so the surface is only defined
# where there are measurements behind it -- which is what lets the comparison in
# script 3 be masked honestly.
surf00 <- readMat("matlab/surfdat2000.mat")

# Grid of record, and the surface on it. Matlab's meshgrid puts salt down rows
# and temperature across columns; transpose to the layout the R scripts use,
# z[i, j] = f(T[i], S[j]).
T_grid_00 <- sort(unique(as.vector(surf00$Tdat)))
S_grid_00 <- sort(unique(as.vector(surf00$Sdat)))
z00 <- t(surf00$Cdat)

n_finite <- sum(is.finite(z00))
message(sprintf("2000 surface (griddata '%s') on the grid of record: %d of %d nodes finite (%.0f%%)",
                as.character(surf00$method), n_finite, length(z00),
                100 * n_finite / length(z00)))

covered <- which(is.finite(z00), arr.ind = TRUE)
message(sprintf("  covering T %.1f-%.1f C, S %.0f-%.0f mmol/kg",
                min(T_grid_00[covered[, 1]]), max(T_grid_00[covered[, 1]]),
                min(S_grid_00[covered[, 2]]), max(S_grid_00[covered[, 2]])))

# Check that the points carried through Matlab are the ones read here, so a
# stale surfdat2000.mat cannot silently disagree with the CSVs in data/.
stopifnot(nrow(cal00) == length(surf00$Craw))
raw_matlab <- tibble(NaCl_mmol_kg = as.vector(surf00$Sraw),
                     Temperature_C = as.vector(surf00$Traw),
                     C = as.vector(surf00$Craw)) %>%
    arrange(NaCl_mmol_kg, Temperature_C)
raw_here <- cal00 %>% arrange(NaCl_mmol_kg, Temperature_C)
if (max(abs(raw_here$C - raw_matlab$C)) > 1e-9 ||
    max(abs(raw_here$Temperature_C - raw_matlab$Temperature_C)) > 1e-9) {
    stop("surfdat2000.mat was built from different data than data/2000_JdF -- rebuild it")
}

# The dense surface over the data's own range, for plotting the 2000 surface on
# its own rather than on the (much larger) grid of record.
T_nat_00 <- sort(unique(as.vector(surf00$Tnat)))
S_nat_00 <- sort(unique(as.vector(surf00$Snat)))
z00_nat <- t(surf00$Cnat)

# Caveat ------------------------------------------------------------------

# Pressure is not constant across this experiment: 220 bar except the 382 C
# block, which ran at 300 bar, where the 2007 run sat near 256 kgf/cm2
# throughout. The two calibrations are therefore not at matched pressure at the
# top of the temperature range.
