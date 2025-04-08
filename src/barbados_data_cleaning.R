# Setup -------------------------------------------------------------------

library(dplyr)
library(tibble)
library(lubridate)
library(rio)
library(reshape2)
library(forcats)
library(readxl)
library(here)
library(stringr)

data_dir <- "data"
prim_station_dir <- paste(data_dir, "station/primary", sep = "/")
proc_station_dir <- paste(data_dir, "station/processed", sep = "/")

husbands <- read_excel(here(prim_station_dir, "barbados", "husbands_rainfall.xlsx"),
                       skip = 3, na = "Tr")

husbands <- husbands %>% 
  mutate(date = as.Date(paste(Year, Month, Day), format = "%Y %m %d"))

husbands <- husbands %>%
  mutate(country = "Barbados", station = "Husbands") %>%
  select(country, station, date, rain = `RR (mm)`)

station_metadata <- tibble(
  country = "Barbados",
  station = "Husbands",
  latitude = 13.148,
  longitude = -59.624)

# Save processed station data ---------------------------------------------

saveRDS(husbands, here(proc_station_dir, "barbados_station.RDS"))

saveRDS(station_metadata, here(proc_station_dir, "barbados_station_metadata.RDS"))
