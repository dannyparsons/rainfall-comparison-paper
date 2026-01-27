library(here)
library(dplyr)
library(ggplot2)
library(lubridate)

source(here("src", "helper_funs.R"))

# Read in the data (change file directory)
zim_five_stations <- read.csv("C:/Users/lclem/Downloads/zim_five_stations.csv")

# Set up for QC cod
zim_five_stations$station <- factor(zim_five_stations$station)
zim_five_stations$date <- as.Date(zim_five_stations$date)

zim_five_stations <- zim_five_stations %>%
  dplyr::filter(date >= as.Date("1979-01-01") & date <= as.Date("2023-06-30") )

zim_five_stations <- zim_five_stations %>% 
  mutate(month = factor(lubridate::month(date)),
         year = lubridate::year(date))

# Rainfall amounts =======================================================================
# Yearly Rainfall
# rainy days (n_rain) and total rainfall (t_rain)
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

# Inventory plot =======================================================================
# There are 46 missing values.
# This includes a long period of NAs in December 1992 in Buffalo Range. Otherwise some random missing values throughout for Mt. Darwin (6 values), Buffalo Range, and Plumtree (3 values)
# The long period in Buffalo Range will probably mean that year will be excluded from annual summaries. 
ggplot(zim_five_stations, aes(x = date, y = station, fill = !is.na(rain))) +
  geom_tile() +
  geom_hline(yintercept = seq(0.5, by = 1, length.out = length(unique(zim_five_stations$station)) + 1))

# Fill Date gaps =======================================================================
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

# Inventory plot again =======================================================================
ggplot(zim_five_stations, aes(x = date, y = station, fill = !is.na(rain))) +
  geom_tile() +
  geom_hline(yintercept = seq(0.5, by = 1, length.out = length(unique(zim_five_stations$station)) + 1))

# Large or negative values check =======================================================================
# Three instances, all look plausible.

# 1. Buffalo Range February 2019: Looks plausible.
# Has a big values on the next day also, 68.
# Nearby stations have high-ish rainfall (Chisumbanje has 24.7mm that day; Masvingo has 20.1mm that day; 25.9mm the next). (The other two are further away and have rainfall ~5mm).
# Cyclone Dineo hit Mozambique over that day: "Widespread flooding took place in Zimbabwe, with Mutare, Chiredzi, and Beitbridge particularly hard-hit.". Chiredzi is ~18km from Buffalo Range

# 2. Chisumbanje, January 2005: Looks plausible.
# In the middle of a 4 day rain spell
# Nearby stations have high-ish rainfall (Buffalo Range 42.9mm; Masvingo 69.0mm)
# Potentially links with Storm Chedza - "heavy rainfall over Mozambique occurred between 11 and 13 January with some stations recording about 80% of their total January 2015 rainfall as resulting from this event."

# 3. Masvingo, March 2003: Looks plausible.
# In the middle of a very big 6 day rain spell
# No major values in Buffalo Range or Chisumbanje around then. Mt. Darwin has 62.7mm rainfall.
# Other data sources show that Rupike and Zaka had rainfall >100mm that day, and they're within 100km of Masvingo
# Cyclone Japhet hit over that time, which seems to have been pretty devastating and was reported to be in the South of Zimbabwe.

large_check <- zim_five_stations %>% 
  filter(rain < 0 | rain > 200)
if(nrow(large_check) > 0) View(large_check)

# Consecutive non-zero values check =======================================================================
# Plausible. Only 2 days cases, which isn't impossible so we can leave these in.
consec_check <- zim_five_stations %>% 
  group_by(station) %>%
  mutate(same = rep(rle(as.numeric(rain))$lengths, rle(as.numeric(rain))$lengths)) %>%
  filter(rain > 1.5 & same >= 2)
if(nrow(consec_check) > 0) View(consec_check)

# Consecutive rain days check =======================================================================
# Plausible
# Consecutive rainy days > 15 days in a row are all in the middle of the rainy season so this seems pretty normal - especially as it's never more than 20.
raindays_check <- zim_five_stations %>%
  group_by(station) %>%
  mutate(raindays = cumsum(rain > 0) - cummax(((rain > 0) == 0) * cumsum(rain > 0))) %>%
  filter(raindays > 15)
if(nrow(raindays_check) > 0) View(raindays_check)

# Dry months check - strict =======================================================================
# Some instances of 0mm rain in some years for the stations in the rainy period. 
# Some suspicious. E.g., Chisumbanje in the 2002-2003 period and 2009-2010 period. 
# Mostly occurring in Nov and Mar - the beginning and end of the rainy seasons. Hard to tell if this is ordinary or not.
# I might want to check with Roger for his view since he knows the region better.

drymonths_check <- zim_five_stations %>%
  filter(!month %in% 4:10) %>%
  group_by(station, year, month) %>%
  summarise(t_rain = sum(rain)) %>%
  filter(t_rain == 0)
if(nrow(drymonths_check) > 0) View(drymonths_check)

display_daily(zim_five_stations %>% filter(station == "Chisumbanje" & year %in% c(2002, 2003)), 
              Stations = "Chisumbanje", Years = 2002:2003, Variables = "rain")

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


# 1. Buffalo Range, 1992: Plausible.
# No rain in February at all.
# No rain from 24th Jan - 12th March (> 0.85mm), this looks plausible. 

# 2. Buffalo Range, 2004 and 2021: Both had no rain in March
# 2004: 0mm, 3mm, 2mm, 0mm rain on February 26th-29th; 13.2mm on April 1st, then <0.85mm until April 11th. Maybe coincidence, but slightly suspicious.
# 2021: 0mm 19th-26th Feb; 15mm 27th Feb; 0.4mm 28th Feb. 0mm for all of April, until May 23rd. Unsure.

* For 2021, no rain until May 23rd that year, so, they were in the dry season?
Might have been a short rainy season. 

# 3. Chisumbanje: No rain in March.
# a) 1982: There's also no rain in April, May, June. May be plausible and not unusual. 
# b) 2005: No rain in March, until November. Slightly unusual but not impossible.
# c) 2008 and 2011: No rain in March, until October. Slightly unusual but not impossible.
# d) 2009: No rain in March for the rest of the year. Suspicious: There's then no rainy season in this year (2009-2010 are all 0 in the rainy season, see 4.)  

# 4. Chisumbanje, 2002-2003; 2009-2010
# There's no rainfall in the 5 month period we are looking at.
# Almost certainly incorrectly entered missing values. We'll probably change this to NA. 
# Chisumbanje is close to Buffalo Range, so I’d have expected to see similar 0mm patterns in Buffalo Range if these values were real.
# From what I can find online, severe drought episodes were recorded in 1991-92, 1994-95, 2002-03, 2015-16, and 2018-19, primarily affecting Matabeleland North and South rather than Manicaland (where Chisumbanje is located). That matches our 2002-03 drought period, but doesn't explain the apparent dryness in 2009-2010.

# Chisumbanje is in Manicaland, and I can see that Manicaland is reported to have experienced a drought in 2010, but it’s unclear whether this refers to the 2009–10 season or 2010–11.

### TODO ###########################################################
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
