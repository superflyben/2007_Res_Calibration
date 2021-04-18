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
cal07_no_dup <- cal07 %>% 
    group_by(NaCl_mmol_kg, Temperature_C) %>% 
    arrange(Pressure_kgf_cm2) %>% 
    summarise_all(first) %>% 
    ungroup()

# Use Matlab benchmark to define output grid
Tdat <- unique(surfdat$Tdat)
Sdat <- unique(surfdat$Sdat)

# Interpolation
# --> Fix to include custom grid 
R_surfdat <- interp(x = cal07_no_dup$Temperature_C,
                    y = cal07_no_dup$NaCl_mmol_kg, 
                    z = cal07_no_dup$C)

# Convert to matrix format for creating surface plot
# NOTE: works with temperature constant across rows and salt constant down columns
#       which is the transpose of what's read in from Matlab. Not sure why 
R_surfdat$T_grid <- matrix(data = R_surfdat$x, nrow = 40, ncol = 40, byrow = FALSE)
R_surfdat$S_grid <- matrix(data = R_surfdat$y, nrow = 40, ncol = 40, byrow = TRUE)

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

# --> pick data section to use in calculating scale factor
#       - refer to  /Users/ben/Dropbox/9North/Calibration/2007/scalefactor.m

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
#
#   