library(here)
library(dplyr)
library(ggplot2)
library(lubridate)

source(here("src", "helper_funs.R"))

daily <- readRDS(here("data", "station", "processed", "barbados_station.RDS"))

husbands <- daily %>% 
  mutate(month = factor(lubridate::month(date)),
         year = lubridate::year(date)) %>%
  filter(station == "Husbands" & year >= 1979) %>%
  dplyr::select(-c(country))

neg_rain <- husbands %>% filter(rain < 0)
if(nrow(neg_rain) > 0) View(neg_rain)
if(nrow(neg_rain) > 0) {
  husbands <- husbands %>% 
    mutate(rain = ifelse(rain < 0, NA, rain))
}

by_year <- husbands %>%
  group_by(station, year) %>%
  summarise(n_rain = sum(rain > 0.85),
            t_rain = sum(rain))

by_month <- husbands %>%
  group_by(station, year, month) %>%
  summarise(n_rain = sum(rain > 0.85),
            t_rain = sum(rain))

# Rain amounts
ggplot(husbands %>% filter(rain > 0), aes(x = month, y = rain)) +
  geom_boxplot()

# Number rain days
ggplot(by_month, aes(x = month, y = n_rain)) +
  geom_boxplot()

# Total rain time series
ggplot(by_month, aes(x = year, y = t_rain)) +
  geom_line() +
  facet_wrap(~month)

# Number rainy days time series
ggplot(by_month, aes(x = year, y = n_rain)) +
  geom_line() +
  facet_wrap(~month)

# Inventory plot
ggplot(husbands, aes(x = date, y = 1, fill = !is.na(rain))) +
  geom_tile()

# Fill Date gaps
dates <- seq(min(husbands$date), max(husbands$date), by = 1)
dd <- data.frame(date = dates)

nr <- nrow(dd)
if(nrow(husbands) < nr) {
  husbands <- full_join(date_df, husbands, by = c("date"))
  print(paste("Filled", nrow(husbands) - nr, "rows"))
  husbands <- husbands %>%
    mutate(year = year(date), month = factor(month(date)))
  # Inventory plot again
  ggplot(husbands, aes(x = date, y = 1, fill = !is.na(rain))) +
    geom_tile()
}

# Large or negative values check
large_check <- husbands %>% 
  filter(rain < 0 | rain > 200)
if(nrow(large_check) > 0) View(large_check)

# Consecutive non-zero values check
consec_check <- husbands %>% 
  mutate(same = rep(rle(as.numeric(rain))$lengths, rle(as.numeric(rain))$lengths)) %>%
  filter(rain > 1.5 & same >= 2)
if(nrow(consec_check) > 0) View(consec_check)

# Consecutive rain days check
raindays_check <- husbands %>%
  mutate(raindays = cumsum(rain > 0) - cummax(((rain > 0) == 0) * cumsum(rain > 0))) %>%
  filter(raindays > 15)
if(nrow(raindays_check) > 0) View(raindays_check)

# Dry months check - all
drymonths_check <- husbands %>%
  group_by(year, month) %>%
  summarise(t_rain = sum(rain)) %>%
  filter(t_rain < 5)
if(nrow(drymonths_check) > 0) View(drymonths_check)

# daily graph
ggplot(husbands %>%
         mutate(year_groups <- cut(year, breaks = seq(1980, 2020, 5))),
  aes(x = date, y = rain)) +
  geom_col(colour = "blue") +
  geom_rug(data = filter(husbands, is.na(rain)), colour = "red") +
  scale_x_date(breaks = "4 month", date_labels = "%b %y") +
  scale_y_continuous(limits = c(0, 100)) + 
  facet_wrap(~year_groups, scales = "free_x", ncol = 1, strip.position = "right")

saveRDS(husbands, here("data", "station", "cleaned", "husbands_qc.RDS"))
