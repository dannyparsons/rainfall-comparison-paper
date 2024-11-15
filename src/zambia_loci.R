library(here)
library(ggplot2)
library(lubridate)
library(reshape2)
library(viridis)
library(RColorBrewer)
library(tidyr)
library(hydroGOF)
library(stringr)
library(knitr)
library(kableExtra)
library(rnaturalearthdata)
library(rnaturalearth)
library(ggrepel)
library(sp)
library(tibble)
library(verification)
library(purrr)
library(dplyr)

source(here("src", "helper_funs.R"))


# Data Import -------------------------------------------------------------

zm <- readRDS(here("data", "station", "cleaned", "zambia_gridded.RDS")) %>%
  dplyr::select(-(arc2_rain))

s_doy_start <- 214
zm <- zm %>% mutate(doy = yday_366(date),
                    s_doy = (doy - s_doy_start + 1) %% 366,
                    s_doy = ifelse(s_doy == 0, 366, s_doy),
                    syear = year,
                    syear = ifelse(s_doy > (366 - s_doy_start + 1), syear - 1, syear),
                    month = factor(month, levels = c(8:12, 1:7)),
                    month_abb = factor(month, labels = month.abb[c(8:12, 1:7)]),
                    rain = ifelse(rain < 0, 0, rain),
                    era5_rain = ifelse(era5_rain < 0, 0, era5_rain),
                    season = ifelse(as.character(month) %in% 5:10, "dry", as.character(month)),
                    season = factor(season, levels = c("dry", 11:12, 1:4))) %>%
  filter(syear >= 1983)

zm_long_st <- zm %>% 
  melt(id.vars = c("station", "date", "year", "syear", "month", "month_abb", 
                   "doy", "s_doy", "season", "rain"),
       measure.vars = names(zm)[endsWith(names(zm), "rain")][-1],
       variable.name = "product", value.name = "pr_rain")

zm_long <- zm %>% 
  melt(id.vars = c("station", "date", "year", "syear", "month", "month_abb", "doy", "s_doy", "season"),
       measure.vars = names(zm)[endsWith(names(zm), "rain")],
       variable.name = "product", value.name = "rain") %>%
  mutate(rainday = rain > 1)

zm_long$product <- recode(zm_long$product, rain = "station")
products <- levels(zm_long$product)
products <- products[-1]
names(products) <- substr(products, 1, nchar(products) - 5)


# Thresholds --------------------------------------------------------------

thresh_station_season <- zm_long_st %>%
  mutate(st_wd0p85 = rain > 0.85) %>%
  filter(!is.na(rain) & !is.na(pr_rain)) %>%
  group_by(product, station, season) %>%
  summarise(prob = mean(st_wd0p85), 
            thres = quantile(pr_rain, 1 - prob)) %>%
  mutate(thres = ifelse(thres < 0.85, 0.85, thres))

thresh_season <- zm_long_st %>%
  mutate(st_wd0p85 = rain > 0.85) %>%
  filter(!is.na(rain) & !is.na(pr_rain)) %>%
  group_by(product, season) %>%
  summarise(prob = mean(st_wd0p85), 
            thres = quantile(pr_rain, 1 - prob)) %>%
  mutate(thres = ifelse(thres < 0.85, 0.85, thres))

ggplot(thresh_station_season, aes(x = season, y = thres, group = station, colour = station)) +
  geom_line() +
  geom_line(aes(x = season, y = thres, group = 1), data = thresh_season, inherit.aes = FALSE, 
            colour = "black", size = 1) +
  facet_wrap(~product)

zm_long_st <- full_join(zm_long_st, thresh_season, by = c("product", "season"))

scale_season1 <- zm_long_st %>%
  filter(!is.na(rain) & !is.na(pr_rain)) %>%
  filter(rain > 0.85) %>%
  group_by(product, season) %>%
  summarise(s_rain = mean(rain)) %>%
  mutate(s_rain = s_rain - 0.85)

scale_season2 <- zm_long_st %>%
  filter(!is.na(rain) & !is.na(pr_rain)) %>%
  filter(pr_rain > thres) %>%
  group_by(product, season) %>%
  summarise(s_prod = mean(pr_rain),
            thres = dplyr::first(thres)) %>%
  mutate(s_prod = s_prod - thres)

scale_season <- full_join(scale_season1, scale_season2, by = c("product", "season"))
scale_season$s <- scale_season$s_rain/scale_season$s_prod
scale_season$s_rain <- NULL
scale_season$s_prod <- NULL

zm_long_st$prob <- NULL
zm_long_st$thres <- NULL

zm_long_st <- full_join(zm_long_st, scale_season, by = c("product", "season"))

zm_long_st <- zm_long_st %>%
  mutate(pr_rain_loci = ifelse(pr_rain <= thres, 0, pr_rain - thres),
         pr_rain_loci = pr_rain_loci * s)


# Checks ------------------------------------------------------------------

by_syear1 <- zm_long_st %>%
  group_by(station, syear) %>%
  filter(product == "chirps_rain") %>%
  summarise(total_rain = sum(naif_nmin(rain, 355)), 
            n_rain = sum(naif_nmin(rain, 355) >= 0.85), 
            mean_rain = total_rain/n_rain
  )

by_syear2 <- zm_long_st %>%
  group_by(station, syear, product) %>%
  summarise(total_rain = sum(naif_nmin(pr_rain, 355)), 
            n_rain = sum(naif_nmin(pr_rain, 355) >= 0.85), 
            mean_rain = total_rain/n_rain
  )

by_syear3 <- zm_long_st %>%
  group_by(station, syear, product) %>%
  summarise(total_rain = sum(naif_nmin(pr_rain_loci, 355)), 
            n_rain = sum(naif_nmin(pr_rain_loci, 355) >= 0.85), 
            mean_rain = total_rain/n_rain
  )

ggplot(by_syear1, aes(x = syear, y = total_rain)) +
  geom_line(aes(colour = "station"), size = 1) +
  geom_line(aes(colour = "product"), by_syear2 %>% filter(product == "era5_rain")) +
  geom_line(aes(colour = "loci"), by_syear3 %>% filter(product == "era5_rain")) +
  facet_wrap(~station)

ggplot(by_syear1, aes(x = syear, y = n_rain)) +
  geom_line(aes(colour = "station"), size = 1) +
  geom_line(aes(colour = "product"), by_syear2 %>% filter(product == "era5_rain")) +
  geom_line(aes(colour = "loci"), by_syear3 %>% filter(product == "era5_rain")) +
  facet_wrap(~station)

ggplot(by_syear1, aes(x = syear, y = mean_rain)) +
  geom_line(aes(colour = "station"), size = 1) +
  geom_line(aes(colour = "product"), by_syear2 %>% filter(product == "chirps_rain")) +
  geom_line(aes(colour = "loci"), by_syear3 %>% filter(product == "chirps_rain")) +
  facet_wrap(~station)


# Start of rains ----------------------------------------------------------

zm_start_station <- zm_long_st %>%
  filter(product == "chirps_rain") %>%
  group_by(station) %>%
  mutate(roll_sum_rain = RcppRoll::roll_sumr(x = rain, n = 3, fill = NA, na.rm = FALSE),
         dry_spell = .spells(rain <= 0.85),
         roll_max_dry_spell = dplyr::lead(RcppRoll::roll_maxl(dry_spell, n = 30, fill = NA))) %>%
  filter((rain >= 0.85 & roll_sum_rain > 20 & roll_max_dry_spell <= 10) | 
           is.na(rain) | is.na(roll_sum_rain) | is.na(roll_max_dry_spell)) %>%
  group_by(station, syear) %>%
  filter(s_doy >= 107 & s_doy <= 213) %>%
  summarise(start_rain_doy_dry = ifelse(is.na(dplyr::first(rain)) | is.na(dplyr::first(roll_sum_rain)) | is.na(dplyr::first(roll_max_dry_spell)),
                                        NA, dplyr::first(s_doy))
  )

zm_year_pr_rain_dry <- zm_long_st %>%
  group_by(product, station) %>%
  mutate(roll_sum_pr_rain = RcppRoll::roll_sumr(x = pr_rain, n = 3, fill = NA, na.rm = FALSE),
         dry_spell = .spells(pr_rain <= 0.85),
         roll_max_dry_spell = dplyr::lead(RcppRoll::roll_maxl(dry_spell, n = 30, fill = NA))) %>%
  filter((pr_rain >= 0.85 & roll_sum_pr_rain > 20 & roll_max_dry_spell <= 10) | 
           is.na(pr_rain) | is.na(roll_sum_pr_rain) | is.na(roll_max_dry_spell)) %>%
  group_by(product, station, syear) %>%
  filter(s_doy >= 107 & s_doy <= 213) %>%
  summarise(start_pr_rain_doy_dry = ifelse(is.na(dplyr::first(pr_rain)) | is.na(dplyr::first(roll_sum_pr_rain)) | is.na(dplyr::first(roll_max_dry_spell)),
                                           NA, dplyr::first(s_doy))
  )

zm_year_pr_rain_loci_dry <- zm_long_st %>%
  group_by(product, station) %>%
  mutate(roll_sum_pr_rain = RcppRoll::roll_sumr(x = pr_rain_loci, n = 3, fill = NA, na.rm = FALSE),
         dry_spell = .spells(pr_rain_loci <= 0.85),
         roll_max_dry_spell = dplyr::lead(RcppRoll::roll_maxl(dry_spell, n = 30, fill = NA))) %>%
  filter((pr_rain_loci >= 0.85 & roll_sum_pr_rain > 20 & roll_max_dry_spell <= 10) | 
           is.na(pr_rain_loci) | is.na(roll_sum_pr_rain) | is.na(roll_max_dry_spell)) %>%
  group_by(product, station, syear) %>%
  filter(s_doy >= 107 & s_doy <= 213) %>%
  summarise(start_pr_rain_loci_doy_dry = ifelse(is.na(dplyr::first(pr_rain_loci)) | is.na(dplyr::first(roll_sum_pr_rain)) | is.na(dplyr::first(roll_max_dry_spell)),
                                           NA, dplyr::first(s_doy))
  )

zm_start <- full_join(zm_year_pr_rain_dry, zm_year_pr_rain_loci_dry, by = c("product", "station", "syear"))
zm_start_era5 <- zm_start %>% filter(product == "era5_rain")
names(zm_start_era5)[4:5] <- c("ERA5", "ERA5 LOCI")
zm_start_era5$product <- NULL

zm_start_era5 <- full_join(zm_start_station, zm_start_era5, by = c("station", "syear"))
names(zm_start_era5)[3] <- "Gauge"
zm_start_era5 <- pivot_longer(zm_start_era5, Gauge:`ERA5 LOCI`, names_to = "Source",
                              names_ptypes = list(Source = factor(levels = c("Gauge", "ERA5", "ERA5 LOCI"))))

g <- ggplot(zm_start_era5, aes(x = syear, y = value, colour = Source, size = Source)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = c("black", c25[1:2])) +
  scale_x_continuous("Year", breaks = seq(1980, 2020, 5)) +
  scale_y_continuous("Start of rains") +
  scale_size_manual(values = c(1, 0.5, 0.5)) +
  facet_wrap(~station)

##! Figure zm_start (Section: 4.3.4)
ggsave(here("results", "zambia_start_rains_era5_loci.png"),
       g, width = 12, height = 6)

zm_start_era5 <- zm_start_station
names(zm_start_era5)[3] <- "start_station"
zm_start_era5 <- full_join(zm_start_era5, zm_start %>% filter(product == "era5_rain"), by = c("station", "syear"))
names(zm_start_era5)[5:6] <- c("start_era5", "start_era5_loci")

zm_start_era5 <- pivot_longer(zm_start_era5, cols = c(start_era5, start_era5_loci),
                              names_to = "type", values_to = "start_prod")

zm_start_metrics <- zm_start_era5 %>%
  group_by(type, station) %>%
  summarise(r = cor(start_prod, start_station, use = "na.or.complete"),
            me = hydroGOF::me(start_prod, start_station),
            mae = hydroGOF::mae(start_prod, start_station),
            rSD = hydroGOF::rSD(start_prod, start_station)
  )

df <- zm_start_metrics %>%
  mutate(r = as.character(round(r, 2)),
         me = as.character(round(me, 1)),
         mae = as.character(round(mae, 1)),
         rSD = as.character(round(rSD, 2))) %>%
  pivot_longer(cols = c("r", "me", "mae", "rSD"), names_to = "metric") %>%
  mutate(metric = factor(metric, levels = c("r", "me", "mae", "rSD"))) %>%
  dplyr::select(station, type, metric, value) %>%
  mutate(value = format(value, digits = 2)) %>%
  pivot_wider(names_from = "type", values_from = "value") %>%
  arrange(metric, station)

##! Table zm_start_metrics (Section: 4.3.4)
df %>% 
  kable() %>%
  collapse_rows(columns = 1) %>%
  skable()
