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

# Zambia ------------------------------------------------------------------

moorings <- readRDS(here(prim_station_dir, "zambia", "moorings.RDS"))
moorings <- moorings$get_data_frame("Moorings")
moorings <- moorings %>% 
  mutate(country = "Zambia", station = "Moorings") %>%
  dplyr::select(country, station, date = Date, rain = rain)

station_metadata <- tibble(
                     country = "Zambia",
                     station = "Moorings",
                     latitude = -16.19,
                     longitude = 27.53)

zam_files <- list.files(here(prim_station_dir, "zambia"), 
                        pattern = "\\.DLY$", full.names = TRUE)
zam_data <- import_list(zam_files, format = "txt", na.strings = "-99.0")
zam_data <- bind_rows(zam_data, .id = "element")
zam_data$element <- sub('.*\\_', '', zam_data$element)
names(zam_data) <- c("element", "station", "year", "month", paste0("d", 1:31))
zam_data <- melt(zam_data, measure.vars = paste0("d", 1:31), variable.name = "day")
zam_data <- dcast(zam_data, station + year + month + day ~ element, value.var = "value")
zam_data <- zam_data %>% 
  mutate(day = substr(day, 2, 4), 
         date = as.Date(paste(year, month, day, sep = "/")),
         station = trimws(station, whitespace = "\'"),
         station = dplyr::recode(station, 
                                 CHOMA001 = "Choma", 
                                 LIVING01 = "Livingstone",
                                 MANSA001 = "Mansa",
                                 KASAMA01 = "Kasama",
                                 MAGOYE01 = "Magoye",
                                 MPIKA001 = "Mpika"),
         country = "Zambia") %>%
  filter(!is.na(date)) %>%
  dplyr::select(country, station, date, rain, tmin, tmax)

zam_data <- bind_rows(zam_data, moorings)

station_metadata <- station_metadata %>% 
  add_row(tibble_row(country = "Zambia",
                     station = "Choma",
                     latitude = -16.84,
                     longitude = 27.07))

station_metadata <- station_metadata %>% 
  add_row(tibble_row(country = "Zambia",
                     station = "Livingstone",
                     latitude = -17.82,
                     longitude = 25.82))

station_metadata <- station_metadata %>% 
  add_row(tibble_row(country = "Zambia",
                     station = "Mansa",
                     latitude = -11.14,
                     longitude = 28.87))

station_metadata <- station_metadata %>% 
  add_row(tibble_row(country = "Zambia",
                     station = "Kasama",
                     latitude = -10.22,
                     longitude = 31.14))

station_metadata <- station_metadata %>% 
  add_row(tibble_row(country = "Zambia",
                     station = "Magoye",
                     latitude = -16.00,
                     longitude = 27.62))

station_metadata <- station_metadata %>% 
  add_row(tibble_row(country = "Zambia",
                     station = "Mpika",
                     latitude = -11.90,
                     longitude = 31.43))

# Save processed station data ---------------------------------------------

zam_data <- zam_data %>%
  mutate(country = stringr::str_to_title(country),
         station = stringr::str_to_title(station))

saveRDS(zam_data, here(proc_station_dir, "zambia_stations.RDS"))

# Summarise Metadata ------------------------------------------------------

station_metadata <- station_metadata %>% 
  mutate(country = stringr::str_to_title(country),
         station = stringr::str_to_title(station))

station_summary <- zam_data %>%
  mutate(year = year(date)) %>%
  group_by(country, station) %>%
  summarise(min_year = min(year),
            max_year = max(year),
            years = max_year - min_year + 1,
            rain_complete = mean(!is.na(rain)),
            tmin_complete = mean(!is.na(tmin)),
            tmax_complete = mean(!is.na(tmax)))

station_metadata <- station_metadata %>% 
  full_join(station_summary, by = c("country", "station"))

saveRDS(station_metadata, here(proc_station_dir, "zambia_stations_metadata.RDS"))
