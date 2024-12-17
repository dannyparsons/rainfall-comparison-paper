library(here)
library(dplyr)
library(ggplot2)
library(lubridate)

source(here("src", "helper_funs.R"))

zambia <- readRDS(here("data", "station", "processed", "zambia_stations_updated.RDS"))
data_dir <- "data"
proc_station_dir <- paste(data_dir, "station/processed", sep = "/")

zambia <- zambia %>% 
  mutate(month = factor(lubridate::month(date)),
         year = lubridate::year(date)) %>%
  filter(country == "Zambia" & year >= 1983) %>%
  dplyr::select(-country)

zambia_year <- zambia %>%
  group_by(station, year) %>%
  summarise(n_rain = sum(rain > 0.85),
            t_rain = sum(rain))

zambia_month <- zambia %>%
  group_by(station, year, month) %>%
  summarise(n_rain = sum(rain > 0.85),
            t_rain = sum(rain))

# Rain amounts
ggplot(zambia %>% filter(rain > 0), aes(x = month, y = rain)) +
  geom_boxplot() +
  facet_wrap(~station)

# Number rain days
ggplot(zambia_month, aes(x = month, y = n_rain)) +
  geom_boxplot() +
  facet_wrap(~station)

# Inventory plot
ggplot(zambia, aes(x = date, y = station, fill = !is.na(rain))) +
  geom_tile() +
  geom_hline(yintercept = seq(0.5, by = 1, length.out = length(unique(zambia$station)) + 1))

# Fill Date gaps
dates_list <- list()
for(s in unique(zambia$station)) {
  
  dates <- seq(min((zambia %>% filter(station == s))$date), 
               max((zambia %>% filter(station == s))$date),
               by = 1)
  dd <- data.frame(station = s, date = dates)
  dates_list[[length(dates_list) + 1]] <- dd
}
date_df <- bind_rows(dates_list)

nr <- nrow(zambia)
zambia <- full_join(date_df, zambia, by = c("station", "date"))
print(paste("Filled", nrow(zambia) - nr, "rows"))

zambia <- zambia %>%
  mutate(year = year(date), month = factor(month(date)))

# Inventory plot again
ggplot(zambia, aes(x = date, y = station, fill = !is.na(rain))) +
  geom_tile() +
  geom_hline(yintercept = seq(0.5, by = 1, length.out = length(unique(zambia$station)) + 1))

# Large or negative values check
large_check <- zambia %>% 
  filter(rain < 0 | rain > 200)
if(nrow(large_check) > 0) View(large_check)

# remove one very large value for now until confirmed
zambia$rain <- ifelse(zambia$rain > 200, NA, zambia$rain)

# Consecutive non-zero values check
consec_check <- zambia %>% 
  group_by(station) %>%
  mutate(same = rep(rle(as.numeric(rain))$lengths, rle(as.numeric(rain))$lengths)) %>%
  filter(rain > 1.5 & same >= 2)
if(nrow(consec_check) > 0) View(consec_check)

# Consecutive rain days check
raindays_check <- zambia %>%
  group_by(station) %>%
  mutate(raindays = cumsum(rain > 0) - cummax(((rain > 0) == 0) * cumsum(rain > 0))) %>%
  filter(raindays > 15)
if(nrow(raindays_check) > 0) View(raindays_check)

# Dry months check - strict
drymonths_check <- zambia %>%
  filter(!month %in% 4:10) %>%
  group_by(station, year, month) %>%
  summarise(t_rain = sum(rain)) %>%
  filter(t_rain == 0)
if(nrow(drymonths_check) > 0) View(drymonths_check)

display_daily(zambia %>% filter(station == "Livingstone" & year == 2019), 
              Stations = "Livingstone", Years = 2019, Variables = "rain")

# Data ends in October 2019. 
# Suggest making 2019 missing until data is updated
zambia <- zambia %>% mutate(rain = replace(rain, 
                                           station == "Livingstone" & 
                                             year == 2019, 
                                           NA))

display_daily(zambia %>% filter(station == "Mansa" & year == 2015), 
              Stations = "Mansa", Years = 2015, Variables = "rain")

#zambia <- zambia %>% mutate(rain = replace(rain, 
#                                           station == "Mansa" & 
#                                             year == 2015, 
#                                           NA))

display_daily(zambia %>% filter(station == "Mansa" & year == 2016), 
              Stations = "Mansa", Years = 2016, Variables = "rain")

# Data ends in 2016
# Suggest making 2016 missing until data is updated
zambia <- zambia %>% mutate(rain = replace(rain, 
                                           station == "Mansa" & 
                                             year == 2016, 
                                           NA))

display_daily(zambia %>% filter(station == "Magoye" & year == 2014), 
              Stations = "Magoye", Years = 2014, Variables = "rain")

# Lots of missing around this period
# Suggest making Nov 2014 missing until data is updated
#zambia <- zambia %>% mutate(rain = replace(rain, 
#                                           station == "Magoye" & 
#                                             year == 2014 &
#                                             month == 11, 
#                                           NA))

# Dry months check - April/October
drymonths_check <- zambia %>%
  filter(month %in% c(4, 10)) %>%
  group_by(station, year, month) %>%
  summarise(t_rain = sum(rain)) %>%
  filter(t_rain == 0)
if(nrow(drymonths_check) > 0) View(drymonths_check)

 #daily graph
 #for(s in unique(zambia$station)) {
#   g <- ggplot(zambia %>% filter(station == s), aes(x = date, y = rain)) +
#     geom_col(colour = "blue") +
#     geom_rug(data = filter(zambia, station == s & is.na(rain)), colour = "red") +
#     scale_y_continuous(limits = c(0, 100))
#   print(g)
# }

saveRDS(zambia, here("data", "station", "cleaned", "zambia_1979_update_qc.RDS")) 

# Summarise Metadata ------------------------------------------------------
zambia$country <- "Zambia"
station_metadata <- readRDS(here(proc_station_dir, "zambia_stations_metadata.RDS"))

station_metadata <- station_metadata %>% 
  mutate(country = stringr::str_to_title(country),
         station = stringr::str_to_title(station))

station_summary <- zambia %>%
  mutate(year = year(date)) %>%
  group_by(country, station) %>%
  summarise(min_year = min(year),
            max_year = max(year),
            years = max_year - min_year + 1,
            rain_complete = mean(!is.na(rain)))

station_metadata <- station_metadata %>% 
  full_join(station_summary, by = c("country", "station"))

saveRDS(station_metadata, here(proc_station_dir, "zambia_stations_metadata_updated.RDS"))
