library(here)
library(dplyr)

source(here("src", "helper_funs.R"))

df1 <- readRDS(here("data", "station", "cleaned", "zambia_1979_update_qc.RDS"))
df2 <- readRDS(here("data","station", "cleaned", "zm_gridded_updated.RDS"))

zambia <- left_join(df1, df2, by = c("station", "date"))

enacts_old <- readRDS(here("data", "station", "cleaned", "zambia_gridded.RDS")) 

enacts_old <- enacts_old %>%
  select(station, date, enacts_rain) %>%
  rename(enacts0_rain = "enacts_rain") %>%
  filter(year(date) >= 1983)

zambia_combined <- library(dplyr)

df_combined <- full_join(zambia, enacts_old, by = c("station", "date"))

saveRDS(df_combined, here("data", "station", "cleaned", "zambia_gridded_merged.RDS"))
