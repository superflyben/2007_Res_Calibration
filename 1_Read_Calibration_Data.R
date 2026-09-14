# Script to read in the experimental data 

# Housekeeping ------------------------------------------------------------

library(readxl)
library(interp)
library(R.matlab)
library(plotly)
library(htmlwidgets)
library(tidyverse)

# Documentation -----------------------------------------------------------

#   - May be more involved than I thought. 
#   - What do we really want to do here?
#   - 3 things
#       1. ability to plug in T and S and get C from calibration data (Pressure is essentially constant)
#       2. ability to plot the surface and show the data points (ideally interactive)
#       3. ability to provide T and C and get S
#       
#   - Packages to explore
#       - akima
#       - interp
#        - glr
#        - pserp
#       

# Read in the data --------------------------------------------------------

# Read in calibration data
# NOTE: Reading in the processed experimental data, so relying on the decisions
#       made at the time of the experiment that went into converting raw data
#       to processed data 
# --> eventually put the raw data file in the git project folder so and repoint
#     to this, so it will be available on github
cal07 <- read_xls(path = "/Users/ben/Dropbox/9North/2004_Alvin/GGG/excel/calibration_data.xls", 
                  sheet = "extrapolate",
                  range = "A1:E46",
                  col_names = TRUE,
                  .name_repair = 'universal')

# Read in Matlab benchmark
# /Users/ben/Dropbox/9North/Calibration/2007/TSCsurface.m
surfdat <- readMat("../Calibration/2007/surfdat.mat")

# R version ---------------------------------------------------------------

# Remove consecutive and trailing periods introduced by name repair
names(cal07) <- names(cal07) %>% 
    gsub("\\.+$", "", .) %>% 
    gsub("\\.+", "_", .)

# Remove duplicate data (keep only corrected data points at the lower pressure)
# NOTE: 700 mmol/kg @ 202.1 C is duplicated at the *same* pressure (256.2), with
#       C = 2.7726 and 2.9852, so pressure alone cannot break that tie. The
#       Matlab pipeline (mmall.csv) uses 2.9852, the later of the two rows, so
#       ties are broken by taking the last matching row in sheet order. This is
#       the only node in the sheet where the tie-break matters; the genuine
#       two-pressure duplicates (379.7 and 395.5 blocks, 256.2 vs 307.4) are
#       resolved by the pressure filter as before.
cal07_no_dup <- cal07 %>% 
    group_by(NaCl_mmol_kg, Temperature_C) %>% 
    filter(Pressure_kgf_cm2 == min(Pressure_kgf_cm2)) %>% 
    slice_tail(n = 1) %>% 
    ungroup()

# Replace the bad 1000 mmol/kg @ 350.4 C measurement
# The spreadsheet retains the measured C = 4.0717, but that point was rejected
# when the Matlab surface was built: TSCsurface.m drops it and reconstructs it by
# linear interpolation across the rest of the 350 C isotherm (C = 3.59725, read
# from c350.csv). Neighbouring isotherms put 1000 mmol/kg near 3.6-3.7, so 4.07
# is anomalous. Reconstruct here the same way, but interpolate the isotherm as it
# appears in the .xls rather than in c350.csv, which is a 3-decimal rounded copy
# -- hence 3.597536 here versus Matlab's 3.597250, a difference of 3e-4 that is
# entirely that rounding.
bad_pt <- with(cal07_no_dup, NaCl_mmol_kg == 1000 & Temperature_C == 350.4)
iso350 <- cal07_no_dup %>% 
    filter(Temperature_C == 350.4, !bad_pt) %>% 
    arrange(NaCl_mmol_kg)
C_recon <- approx(x = iso350$NaCl_mmol_kg, y = iso350$C, xout = 1000)$y
message(sprintf(
    "1000 mmol/kg @ 350.4 C: replacing measured C = %.6f with reconstructed %.6f (linear interpolation across the 350.4 C isotherm in calibration_data.xls)",
    cal07_no_dup$C[bad_pt], C_recon))
cal07_no_dup$C[bad_pt] <- C_recon
cal07_no_dup$Rx1[bad_pt] <- NA_real_

# Output grid --------------------------------------------------------------

# Every surface in this project must be evaluated on an explicitly specified
# grid. interp()'s default is 40x40 spanning whatever data it was handed, which
# is neither comparable to the Matlab benchmark's 20x20 nor shared with the 2000
# Juan de Fuca surface -- and two surfaces cannot be subtracted node for node
# unless they sit on the same grid.
#
# The grid of record is the one the Matlab surface already uses: 20x20 over
# T 202-395.5 C and S 100-1500 mmol/kg. Taking it from surfdat.mat rather than
# restating the numbers keeps the benchmark comparison exact by construction.
T_grid <- sort(unique(as.vector(surfdat$Tdat)))
S_grid <- sort(unique(as.vector(surfdat$Sdat)))

# Interpolation
R_surfdat <- interp(x = cal07_no_dup$Temperature_C,
                    y = cal07_no_dup$NaCl_mmol_kg, 
                    z = cal07_no_dup$C,
                    xo = T_grid,
                    yo = S_grid,
                    extrap = FALSE)

# Convert to matrix format for creating surface plot
# NOTE: interp returns z[i, j] = f(x[i], y[j]), so temperature varies down rows
#       and salt across columns -- the transpose of what's read in from Matlab.
R_surfdat$T_grid <- matrix(data = R_surfdat$x, 
                           nrow = length(R_surfdat$x), ncol = length(R_surfdat$y), 
                           byrow = FALSE)
R_surfdat$S_grid <- matrix(data = R_surfdat$y, 
                           nrow = length(R_surfdat$x), ncol = length(R_surfdat$y), 
                           byrow = TRUE)

# Compare against the Matlab benchmark, node for node
# Both surfaces are now on the same grid, so this is a like-for-like difference.
# NOTE: interp() defaults to linear interpolation while TSCsurface.m uses
#       griddata(...,'cubic'), so a residual here is expected and is a *method*
#       difference, not a data one -- the inputs now agree to mmall.csv's rounding.
R_vs_matlab <- R_surfdat$z - t(surfdat$Cdat)
message(sprintf(
    "R vs Matlab on the shared %dx%d grid: %d of %d nodes finite in both, max |dC| = %.4f, RMS = %.4f",
    length(T_grid), length(S_grid),
    sum(is.finite(R_vs_matlab)), length(R_vs_matlab),
    max(abs(R_vs_matlab), na.rm = TRUE),
    sqrt(mean(R_vs_matlab^2, na.rm = TRUE))))

# Visualize the surface
p <- plot_ly(x = R_surfdat$T_grid, 
        y = R_surfdat$S_grid, 
        z = R_surfdat$z, 
        type = "surface") %>% 
    add_trace(x = cal07_no_dup$Temperature_C,
              y = cal07_no_dup$NaCl_mmol_kg,
              z = cal07_no_dup$C,
              mode = "markers",
              type = "scatter3d", 
              marker = list(size = 5, color = "red", symbol = 104))

htmlwidgets::saveWidget(as_widget(p), "Calibration_Surface_2007.html")

# Matlab Benchmark --------------------------------------------------------

# Visualize the surface
plot_ly(x = surfdat$Tdat, 
        y = surfdat$Sdat, 
        z = surfdat$Cdat, 
        type = "surface") %>% 
    add_trace(x = cal07_no_dup$Temperature_C,
              y = cal07_no_dup$NaCl_mmol_kg,
              z = cal07_no_dup$C,
              mode = "markers",
              type = "scatter3d", 
              marker = list(size = 5, color = "red", symbol = 104))


# Next Steps --------------------------------------------------------------

# Very next step:
# - Find the raw data and Matlab processing script, and clean up just enough
#   to include in the repo.

# --> pick data section to use in calculating scale factor
#       - refer to  /Users/ben/Dropbox/9North/Calibration/2007/scalefactor.m
#       - ok, but what data are we using?
#       - simplify by using cleaned output from Matlab Script
#       - Clean up Matlab code just enough to produce a file that can be ingested by R
#       - Main point of the R project is to demo the calibration process
#       - NOT the rando data cleaning or outdataed phase separation analysis
#       - Could potentially include one figure to illustrate utility
#       - This does not need to be a Nature paper, just functional.

# --> Apply scale factors and convert conductivity to salt concentration 
#       - /Users/ben/Dropbox/9North/Calibration/2007/adjust.m
#       - /Users/ben/Dropbox/9North/Calibration/2007/convert_CtoS.m

# --> Visualize surface
#       - show calibration points
#       - show position on surface that gives measured chloride

# --> Interactive visualization 
#       - one tab with calibration demo including data stream and data points on surface
#       - one tab with date entry and highlighted points on surface and on data

# Some thoughts
# - not sure it's worth going all the way back to Bio9 Prime data to replicate
#   the steps from raw data to calibrated data because some additional complexity
#   in correcting the temperature reading and electronic gain which are not
#   things that anybody using the OOI data would have access to. 
# - Probably better to use the 2015 EGU project as the template, b/c this used
#   observatory data with the 2007 calibration surface, which is the use case
#   we're trying to replicate
# - Also doing it this way would mean not having to revisit any of the data wrangling
#   to the clean bio9 prime data
# - So, the example to include with the Github Repo would be the data pulled for 
#   the 2015 EGU paper. 
# - Could separately include the process used to get the Bio9 Prime data into 
#   shape, or just include the cleaned data without re-doing the processing for
#   the guthub repo.

# So what does this project look like?
# - Motivation at this point is to put a pin in it
# - Perfect world: calibration process on Git hub wih 2 datasets
# - each datset is published and can be linked to a paper
# - EGU is done
# - 9N Should be bare f'ing minimum to get into a journal
# - This is NOT a modeling paper. 
# - Would anybody ever use this or care about it?
# - Or is it just for me and Marv?
# - How much work would be involved
# - In the beginng: a lot!
# - Simplify at every available opportunity 
# - Funded grant would make this so much more doable. 
# - Maybe don't worry about publishing
# - Just finish it off as a Github repo wit two examples and move on
# - Or self publish one or two graphics and be done.
# - Still involves recreating the steps in the Matlab script
# - Or just use the cleaned output and cite the script. 
# - Because nobody will bother to inspect that work
# - And just because it could be re-done in R doesn't mean it should be
# - Also the repo can contain Matlab scripts for some of the analysis
# - And keep in mind, nobody will look at this.