library(here)
library(dplyr)
library(ggplot2)
library(lubridate)

source(here("src", "helper_funs.R"))

daily <- readRDS(here("data", "station", "processed", "stations_all.RDS"))

westafrica <- daily %>% 
  mutate(month = factor(lubridate::month(date)),
         year = lubridate::year(date)) %>%
  filter(country %in% c("Niger", "Ghana") & year >= 1979) %>%
  select(-region)

westafrica <- westafrica %>% arrange(country, station, date)

metadata_station <- readRDS(here("data", "station", "processed", "stations_metadata.RDS"))
metadata_wa <- metadata_station %>%
  filter(country %in% c("Niger", "Ghana"))
rm(metadata_station)

sf_wa <- ne_countries(returnclass = "sf")
ggplot(sf_wa) + 
  geom_sf() +
  geom_point(data = metadata_wa, aes(x = longitude, y = latitude)) +
  geom_text_repel(data = metadata_wa, aes(x = longitude, y = latitude, label = station)) +
  coord_sf(xlim = c(-8, 15), ylim = c(0, 25), expand = FALSE)

rm(sf_wa)

ggsave(here("results", "westafrica", "stations.png"), width = 12, height = 6)

westafrica <- westafrica %>% mutate(rain = na_if(rain, -99.9))

neg_rain <- westafrica %>% filter(rain < 0)
if(nrow(neg_rain) > 0) View(neg_rain)
if(nrow(neg_rain) > 0) {
  westafrica <- westafrica %>% 
    mutate(rain = ifelse(rain < 0, NA, rain))
  }

by_year <- westafrica %>%
  group_by(station, year) %>%
  summarise(n_rain = sum(rain > 0.85),
            t_rain = sum(rain))

by_month <- westafrica %>%
  group_by(station, year, month) %>%
  summarise(n_rain = sum(rain > 0.85),
            t_rain = sum(rain))

# Rain amounts
ggplot(westafrica %>% filter(rain > 0), aes(x = month, y = rain)) +
  geom_boxplot() +
  facet_wrap(~station)

# Number rain days
ggplot(by_month, aes(x = month, y = n_rain)) +
  geom_boxplot() +
  facet_wrap(~station)

# Inventory plot
ggplot(westafrica, aes(x = date, y = station, fill = !is.na(rain))) +
  geom_tile() +
  geom_hline(yintercept = seq(0.5, by = 1, length.out = length(unique(westafrica$station)) + 1)) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y")

ggsave(here("results", "westafrica", "inventory.png"), width = 12, height = 6)

# Fill Date gaps
dates_list <- list()
for(s in unique(westafrica$station)) {
  
  dates <- seq(min((westafrica %>% filter(station == s))$date), 
               max((westafrica %>% filter(station == s))$date),
               by = 1)
  dd <- data.frame(station = s, date = dates)
  dates_list[[length(dates_list) + 1]] <- dd
}
date_df <- bind_rows(dates_list)

nr <- nrow(date_df)
if(nrow(westafrica) < nr) {
  print(paste("Filling", nrow(westafrica) - nr, "rows"))
  westafrica <- full_join(date_df, westafrica, by = c("station", "date"))
  westafrica <- westafrica %>%
    mutate(year = year(date), month = factor(month(date)))
  # Inventory plot again
  ggplot(westafrica, aes(x = date, y = station, fill = !is.na(rain))) +
    geom_tile() +
    geom_hline(yintercept = seq(0.5, by = 1, length.out = length(unique(westafrica$station)) + 1))
}

# Large or negative values check
large_check <- westafrica %>% 
  filter(rain < 0 | rain > 200)
if(nrow(large_check) > 0) View(large_check)
# High rain amount on 2005/06/05 appear to be correct from reports

# Consecutive non-zero values check
consec_check <- westafrica %>% 
  group_by(station) %>%
  mutate(same = rep(rle(as.numeric(rain))$lengths, rle(as.numeric(rain))$lengths)) %>%
  filter(rain > 1.5 & same >= 2)
if(nrow(consec_check) > 0) View(consec_check)

# Consecutive rain days check
raindays_check <- westafrica %>%
  group_by(station) %>%
  mutate(raindays = cumsum(rain > 0) - cummax(((rain > 0) == 0) * cumsum(rain > 0))) %>%
  filter(raindays > 15)
if(nrow(raindays_check) > 0) View(raindays_check)

# Dry months check - strict
drymonths_check <- westafrica %>%
  filter(month %in% 6:9) %>%
  group_by(station, year, month) %>%
  summarise(t_rain = sum(rain, na.rm = TRUE)) %>%
  filter(t_rain == 0)
if(nrow(drymonths_check) > 0) View(drymonths_check)

display_daily(westafrica %>% filter(station == "Tamale" & year == 2015), 
              Stations = "Tamale", Years = 2015, Variables = "rain")

# remove Sadore 2015
westafrica <- westafrica %>% filter(!(year == 2015 & station == "Sadore"))

# remove most suspicious years until confirmed by Met
# westafrica <- westafrica %>% mutate(rain = replace(rain, 
#                                            station == "Mhlume" & 
#                                              year == 2012, 
#                                            NA))

# daily graph
for(s in unique(westafrica$station)) {
  g <- ggplot(westafrica, aes(x = date, y = rain)) +
    geom_col(colour = "blue") +
    geom_rug(data = filter(westafrica, is.na(rain)), colour = "red") +
    scale_y_continuous(limits = c(0, 100))
  print(g)
}

saveRDS(westafrica, here("data", "station", "cleaned", "westafrica_1979_qc.RDS"))
