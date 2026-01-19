library(here)
library(dplyr)
library(ggplot2)
library(lubridate)

source(here("src", "helper_funs.R"))

zim_five_stations <- read.csv("C:/Users/lclem/Downloads/zim_five_stations.csv")

zim_five_stations$station <- factor(zim_five_stations$station)
zim_five_stations$date <- as.Date(zim_five_stations$date)

zim_five_stations <- zim_five_stations %>%
  dplyr::filter(date >= as.Date("1979-01-01") & date <= as.Date("2023-06-30") )

zim_five_stations <- zim_five_stations %>% 
  mutate(month = factor(lubridate::month(date)),
         year = lubridate::year(date))

# Yearly rainfall - rainy days (n_rain) and total rainfall (t_rain)
zim_five_stations_year <- zim_five_stations %>%
  group_by(station, year) %>%
  summarise(n_rain = sum(rain > 0.85),
            t_rain = sum(rain))

# Monthly-yearly rainfall - rain days and total rainfall
zim_five_stations_month <- zim_five_stations %>%
  group_by(station, year, month) %>%
  summarise(n_rain = sum(rain > 0.85),
            t_rain = sum(rain))

# Rain amounts
ggplot(zim_five_stations %>% filter(rain > 0), aes(x = month, y = rain)) +
  geom_boxplot() +
  facet_wrap(~station)

# Number rain days
ggplot(zim_five_stations_month, aes(x = month, y = n_rain)) +
  geom_boxplot() +
  facet_wrap(~station)

# Inventory plot
ggplot(zim_five_stations, aes(x = date, y = station, fill = !is.na(rain))) +
  geom_tile() +
  geom_hline(yintercept = seq(0.5, by = 1, length.out = length(unique(zim_five_stations$station)) + 1))

# Fill Date gaps
dates_list <- list()
for(s in unique(zim_five_stations$station)) {
  
  dates <- seq(min((zim_five_stations %>% filter(station == s))$date), 
               max((zim_five_stations %>% filter(station == s))$date),
               by = 1)
  dd <- data.frame(station = s, date = dates)
  dates_list[[length(dates_list) + 1]] <- dd
}
date_df <- bind_rows(dates_list)

nr <- nrow(zim_five_stations)
zim_five_stations <- full_join(date_df, zim_five_stations, by = c("station", "date"))
print(paste("Filled", nrow(zim_five_stations) - nr, "rows"))

zim_five_stations <- zim_five_stations %>%
  mutate(year = year(date), month = factor(month(date)))

# Inventory plot again
ggplot(zim_five_stations, aes(x = date, y = station, fill = !is.na(rain))) +
  geom_tile() +
  geom_hline(yintercept = seq(0.5, by = 1, length.out = length(unique(zim_five_stations$station)) + 1))

# Large or negative values check
large_check <- zim_five_stations %>% 
  filter(rain < 0 | rain > 200)
if(nrow(large_check) > 0) View(large_check)

# remove one very large value for now until confirmed
#zim_five_stations$rain <- ifelse(zim_five_stations$rain > 200, NA, zim_five_stations$rain)

# Consecutive non-zero values check
consec_check <- zim_five_stations %>% 
  group_by(station) %>%
  mutate(same = rep(rle(as.numeric(rain))$lengths, rle(as.numeric(rain))$lengths)) %>%
  filter(rain > 1.5 & same >= 2)
if(nrow(consec_check) > 0) View(consec_check)

# Consecutive rain days check
raindays_check <- zim_five_stations %>%
  group_by(station) %>%
  mutate(raindays = cumsum(rain > 0) - cummax(((rain > 0) == 0) * cumsum(rain > 0))) %>%
  filter(raindays > 15)
if(nrow(raindays_check) > 0) View(raindays_check)

# Dry months check - strict
drymonths_check <- zim_five_stations %>%
  filter(!month %in% 4:10) %>%
  group_by(station, year, month) %>%
  summarise(t_rain = sum(rain)) %>%
  filter(t_rain == 0)
if(nrow(drymonths_check) > 0) View(drymonths_check)

display_daily(zim_five_stations %>% filter(station == "Chisumbanje" & year %in% c(2002, 2003)), 
              Stations = "Chisumbanje", Years = 2002:2003, Variables = "rain")

# FROM Previous script: 
# # Data ends in October 2019. 

# 2005, 2008, 2009 MARCH

# Suggest making 2019 missing until data is updated
display_daily(zim_five_stations %>% filter(station == "Buffalo_Range" & year == 2021),
              Stations = "Buffalo_Range", Years = 2021, Variables = "rain")

display_daily(zim_five_stations %>% filter(station == "Chisumbanje" & year == 2016),
              Stations = "Chisumbanje", Years = 2016, Variables = "rain")

display_daily(zim_five_stations %>% filter(station == "Masvingo" & year == 2004),
              Stations = "Masvingo", Years = 2004, Variables = "rain")

display_daily(zim_five_stations %>% filter(station == "Mt_Darwin" & year == 2020),
              Stations = "Mt_Darwin", Years = 2020, Variables = "rain")

display_daily(zim_five_stations %>% filter(station == "Plumtree" & year == 2019),
              Stations = "Plumtree", Years = 2019, Variables = "rain")


# 
# #zim_five_stations <- zim_five_stations %>% mutate(rain = replace(rain, 
# #                                           station == "Mansa" & 
# #                                             year == 2015, 
# #                                           NA))
# 
# display_daily(zim_five_stations %>% filter(station == "Mansa" & year == 2016), 
#               Stations = "Mansa", Years = 2016, Variables = "rain")
# 
# # Data ends in 2016
# # Suggest making 2016 missing until data is updated
# zim_five_stations <- zim_five_stations %>% mutate(rain = replace(rain, 
#                                            station == "Mansa" & 
#                                              year == 2016, 
#                                            NA))
# 
# display_daily(zim_five_stations %>% filter(station == "Magoye" & year == 2014), 
#               Stations = "Magoye", Years = 2014, Variables = "rain")
# 
# Lots of missing around this period
# Suggest making Nov 2014 missing until data is updated
# zim_five_stations <- zim_five_stations %>% mutate(rain = replace(rain,
#                                           station == "Magoye" &
#                                             year == 2014 &
#                                             month == 11,
#                                           NA))

# Dry months check - April/October
drymonths_check <- zim_five_stations %>%
  filter(month %in% c(4, 10)) %>%
  group_by(station, year, month) %>%
  summarise(t_rain = sum(rain)) %>%
  filter(t_rain == 0)
if(nrow(drymonths_check) > 0) View(drymonths_check)

#daily graph
#for(s in unique(zim_five_stations$station)) {
#   g <- ggplot(zim_five_stations %>% filter(station == s), aes(x = date, y = rain)) +
#     geom_col(colour = "blue") +
#     geom_rug(data = filter(zim_five_stations, station == s & is.na(rain)), colour = "red") +
#     scale_y_continuous(limits = c(0, 100))
#   print(g)
# }

#saveRDS(zim_five_stations, here("data", "station", "cleaned", "zim_five_stations_1979_update_qc.RDS")) 

# # Summarise Metadata ------------------------------------------------------
# zim_five_stations$country <- "zim_five_stations"
# station_metadata <- readRDS(here(proc_station_dir, "zim_five_stations_stations_metadata.RDS"))
# 
# station_metadata <- station_metadata %>% 
#   mutate(country = stringr::str_to_title(country),
#          station = stringr::str_to_title(station))
# 
# station_summary <- zim_five_stations %>%
#   mutate(year = year(date)) %>%
#   group_by(country, station) %>%
#   summarise(min_year = min(year),
#             max_year = max(year),
#             years = max_year - min_year + 1,
#             rain_complete = mean(!is.na(rain)))
# 
# station_metadata <- station_metadata %>% 
#   full_join(station_summary, by = c("country", "station"))
# 
# saveRDS(station_metadata, here(proc_station_dir, "zim_five_stations_stations_metadata_updated.RDS"))
