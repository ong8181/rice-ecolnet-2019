####
#### Data compile
#### 2021.1.12 Ushio
####

# Load library
library(tidyverse); packageVersion("tidyverse") # 1.3.0, 2021.1.12
library(lubridate); packageVersion("lubridate") # 1.7.9.2, 2021.1.12

# Load data
d_gr0 <- read.csv("data/data_original/2019_CER_MonitoringSheet.csv")
d_hd0 <- read.csv("data/data_original/2019_CER_RiceHeading.csv")
d_yd0 <- read.csv("data/data_original/2019_CER_YieldData.csv")

# Compile rice growth data
dim(d_gr0); colnames(d_gr0)
d_gr1 <- data.frame(date = ymd(d_gr0$Date),
                    date_time = ymd_hm(paste(d_gr0$Date, d_gr0$Time)),
                    plot = d_gr0$Plot.No.,
                    treatment = d_gr0$Treatment,
                    spad1 = d_gr0$SPAD1,
                    spad2 = d_gr0$SPAD2,
                    spad3 = d_gr0$SPAD3,
                    height1 = d_gr0$Height1,
                    height2 = d_gr0$Height2,
                    height3 = d_gr0$Height3)
d_hd1 <- data.frame(date = ymd(d_gr0$Date),
                    date_time = ymd_hm(paste(d_gr0$Date, d_gr0$Time)),
                    plot = d_gr0$Plot.No.,
                    treatment = d_gr0$Treatment,
                    n_stem1 = d_gr0$N.Stem1,
                    n_stem2 = d_gr0$N.Stem2,
                    n_stem3 = d_gr0$N.Stem3,
                    head1 = d_gr0$Head1,
                    head2 = d_gr0$Head2,
                    head3 = d_gr0$Head3) %>% na.omit()

# Calculate growth rate (cm/24hours)
d_gr1$gr1 <- (d_gr1$height1 - dplyr::lag(d_gr1$height1, n = 10))/as.numeric((d_gr1$date_time - dplyr::lag(d_gr1$date_time, n = 10))/24)
d_gr1$gr2 <- (d_gr1$height2 - dplyr::lag(d_gr1$height2, n = 10))/as.numeric((d_gr1$date_time - dplyr::lag(d_gr1$date_time, n = 10))/24)
d_gr1$gr3 <- (d_gr1$height3 - dplyr::lag(d_gr1$height3, n = 10))/as.numeric((d_gr1$date_time - dplyr::lag(d_gr1$date_time, n = 10))/24)
d_gr1$spad_r1 <- (d_gr1$spad1 - dplyr::lag(d_gr1$spad1, n = 10))/as.numeric((d_gr1$date_time - dplyr::lag(d_gr1$date_time, n = 10))/24)
d_gr1$spad_r2 <- (d_gr1$spad2 - dplyr::lag(d_gr1$spad2, n = 10))/as.numeric((d_gr1$date_time - dplyr::lag(d_gr1$date_time, n = 10))/24)
d_gr1$spad_r3 <- (d_gr1$spad3 - dplyr::lag(d_gr1$spad3, n = 10))/as.numeric((d_gr1$date_time - dplyr::lag(d_gr1$date_time, n = 10))/24)

## Extract and save compiled dataset
d_gr2 <- tibble(d_gr1[d_gr1$plot != "NC",] %>% .[10:nrow(.),])
write.csv(d_gr2, "data/data_rice_growth.csv", row.names = F)
saveRDS(d_gr2, "data/data_rice_growth.obj")

## Save heading data
write.csv(d_hd1, "data/data_rice_heading.csv", row.names = F)
saveRDS(d_hd1, "data/data_rice_heading.obj")

# Compile yield data
write.csv(d_yd0, "data/data_rice_yield.csv", row.names = F)
saveRDS(d_yd0, "data/data_rice_yield.obj")

