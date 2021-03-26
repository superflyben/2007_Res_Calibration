# Script to read in the experimental data 
library(readxl)
library(tidyverse)

# Read in calibration data
# NOTE: Reading in the processed experimental data, so relying on the decisions
#       made at the time of the experiment that went into converting raw data
#       to processed data 
cal07 <- read_xls(path = "/Users/ben/Dropbox/9North/2004_Alvin/GGG/excel/calibration_data.xls", 
                  sheet = "extrapolate",
                  range = "A1:E46",
                  col_names = TRUE, )

# Fix the names
cal07 <- cal07 %>% 
    rename()

# --> grid or otherwise interpolate to create an interpolation surface
#       - refer to /Users/ben/Dropbox/9North/Calibration/2007/TSCsurface.m

# --> pick data section to use in calculating scale factor
#       - refer to  /Users/ben/Dropbox/9North/Calibration/2007/scalefactor.m

# --> Apply scale factors anf convert conductivity to salt concentration 
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
#   the steps from raw data to calibrated data becasue some additional complexity
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