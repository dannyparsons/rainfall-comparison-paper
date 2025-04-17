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

zm <- readRDS(here("data", "station", "cleaned", "husbands_gridded.RDS"))
zm <- zm %>% 
  dplyr::select(station, date, rain, tmax, tmin, month, year, year_groups, chirps_rain, era5_rain)
# 1 Jan = 1
s_doy_start <- 1
zm <- zm %>% mutate(doy = yday_366(date),
                    s_doy = (doy - s_doy_start + 1) %% 366,
                    s_doy = ifelse(s_doy == 0, 366, s_doy),
                    syear = year,
                    syear = ifelse(s_doy > (366 - s_doy_start + 1), syear - 1, syear),
                    month = factor(month, levels = c(1:12)),
                    month_abb = factor(month, labels = month.abb)
                    )

zm_long_st <- zm %>% 
  melt(id.vars = c("station", "date", "year", "syear", "month", "month_abb", 
                   "doy", "s_doy", "rain"),
       measure.vars = names(zm)[endsWith(names(zm), "rain")][-1],
       variable.name = "product", value.name = "pr_rain")

zm_long <- zm %>% 
  melt(id.vars = c("station", "date", "year", "syear", "month", "month_abb", "doy", "s_doy"),
       measure.vars = names(zm)[endsWith(names(zm), "rain")],
       variable.name = "product", value.name = "rain") %>%
  mutate(rainday = rain > 1)

zm_long$product <- recode(zm_long$product, rain = "station")
stations <- c("Husbands")
products <- levels(zm_long$product)
products <- products[-1]
names(products) <- substr(products, 1, nchar(products) - 5)

metadata_station <- readRDS(here("data", "station", "processed", "stations_metadata.RDS"))
metadata_zm <- metadata_station %>%
  filter(country == "Barbados")
rm(metadata_station)

metadata_zm$station <- factor(metadata_zm$station, levels = stations)
zm_long$station <- factor(zm_long$station, levels = stations)

by_station <- zm_long %>%
  group_by(station) %>%
  filter(!is.na(rain) & product == "station") %>%
  summarise(first_date = first(date),
            last_date = last(date))

metadata_zm <- left_join(metadata_zm, by_station, by = "station")

rm(by_station)

zm_long <- left_join(zm_long, metadata_zm, by = "station")

zm_long <- zm_long %>% 
  filter(date >= first_date & date <= last_date)

skable <- function(kable_input) {
  kable_input %>% kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                                full_width = FALSE)
}

## ----prop_complete_days-------------------------------------------------------------
zm %>% 
  group_by(station) %>% 
  summarise(prop_n = 100 * (1 - naflex::na_prop(rain))) %>%
  kable(digits = 1) %>%
  skable()

## ----na_funs------------------------------------------------------------------------
naif_nmin <- function(x, n_min) {
  if(length(na.omit(x)) > 0 && sum(!is.na(x)) >= n_min) {
    na.omit(x)
  } else NA
}

## ----gof_fun------------------------------------------------------------------------
dgof <- function(df, sim, obs, na.rm = TRUE) {
  g <- hydroGOF::gof(sim = df[[sim]], obs = df[[obs]], na.rm = na.rm)
  as.list(g[ ,1])
}

comp_stats <- c("r", "ME", "PBIAS %", "MAE", "NSE", "rSD")
names(comp_stats) <- c("Correlation coefficient (1 = Perfect)",
                       "Mean bias (same units)",
                       "Percentage bias (%)",
                       "Mean absolute bias (same units)",
                       "Nash-Sutcliffe efficiency (1 = Perfect)",
                       "Ratio of standard deviations (< 1 less variable, > 1 more variable)")

comp_stats_digits <- c(2, 0, 0, 0, 2, 3)

## ----yearly_calcs-------------------------------------------------------------------
by_syear <- zm_long %>%
  group_by(station, syear, product) %>%
  summarise(total_rain = sum(naif_nmin(rain, 355)), 
            n_rain0p1 = sum(naif_nmin(rain, 355) > 0.1),
            n_rain0p85 = sum(naif_nmin(rain, 355) > 0.85),
            max_rain = max(naif_nmin(rain, 350)), 
            mean_rain = total_rain/n_rain0p1, 
            n_na = sum(is.na(rain))
            )

by_syear_st <- by_syear %>%
  pivot_wider(id_cols = c(station, syear), 
              names_from = product, values_from = total_rain:n_na, names_sep = "__") %>%
  pivot_longer(cols = -c(station, syear, ends_with("station")), 
               names_to = c(".value", "product"), names_sep = "__")

gof_syear <- by_syear_st %>%
  group_by(station, product) %>%
  nest() %>%
  mutate(n = purrr::map_int(data, 
                        ~sum(!is.na(.$total_rain__station) & !is.na(.$total_rain))),
         gof__total_rain = purrr::map(data, dgof, "total_rain", "total_rain__station", 
                                      na.rm = TRUE),
         gof__n_rain0p85 = purrr::map(data, dgof, "n_rain0p85", "n_rain0p85__station", 
                                  na.rm = TRUE),
         gof__mean_rain = purrr::map(data, dgof, "mean_rain", "mean_rain__station", 
                                     na.rm = TRUE)
         )

gof_pr <- gof_syear %>% 
  unnest(cols = data) %>% 
  group_by(product) %>% 
  nest()

## ----yearly_plots_fun---------------------------------------------------------------
yearly_plots <- function(df, gof_col, stat_pr, stat_st, product_name) {
  max_y <- max(c(df[[stat_pr]], df[[stat_st]]), na.rm = TRUE)
  vals_relace <- c(product_name, "station")
  names(vals_relace) <- c(stat_pr, stat_st)
  dat <- df %>% 
    filter(!is.na(.data[[stat_pr]]) & !is.na(.data[[stat_st]])) %>%
    pivot_longer(cols = c(stat_pr, stat_st), names_to = "product", values_to = stat_pr) %>%
    mutate(product = recode(product, !!!vals_relace),
           ME = purrr::map_dbl(.data[[gof_col]], "ME"),
           r = purrr::map_dbl(.data[[gof_col]], "r"),
           rSD = purrr::map_dbl(.data[[gof_col]], "rSD")
           )
  mean_df <- dat %>% 
    group_by(station, product) %>% 
    summarise(m = mean(.data[[stat_pr]], na.rm = TRUE))
  g <- ggplot(dat, aes(x = syear, y = .data[[stat_pr]], colour = product)) +
    geom_line() +
    geom_hline(data = mean_df, aes(yintercept = m, colour = product)) +
    scale_x_continuous(limits = c(1979, 2020)) +
    # n
    geom_text(data = dat, aes(label = paste("n", n)), size = 4,
              x = 1979, y = max_y, na.rm = TRUE, 
              inherit.aes = FALSE) +
    # bias
    geom_text(data = dat, aes(label = paste("bias", signif(ME, 2))), 
              size = 4, x = 1979 + 6, y = max_y, na.rm = TRUE, 
              inherit.aes = FALSE) +
    # cor
    geom_text(aes(label = paste("cor", round(r, 2))), 
              size = 4, x = 1979 + 12, y = max_y, na.rm = TRUE, 
              inherit.aes = FALSE) +
    # rSD
    geom_text(aes(label = paste("rSD", round(rSD, 2))), 
              size = 4, x = 1979 + 18, y = max_y, na.rm = TRUE, 
              inherit.aes = FALSE) +
    ggtitle(paste(stat_pr, ":", product_name, "vs", "station")) +
    facet_wrap(~station)
  g
}

## ----stats_tables_fun---------------------------------------------------------------
stats_tables <- function(df, obj_col, obj_stats = comp_stats) {
  for(i in seq_along(obj_stats)) {
    dat <- df %>% 
      ungroup() %>%
      mutate(purrr::map_dbl(df[[obj_col]], obj_stats[i]))
    names(dat)[ncol(dat)] <- obj_stats[i]
    dat <- dat %>% 
      mutate(station = as.character(station)) %>%
      pivot_wider(id_cols = "station", names_from = "product", 
                  values_from = obj_stats[i]) %>%
      data.frame()
  dat[nrow(dat) + 1, ] <- c(list("(All)"), as.list(as.numeric(colMeans(dat[ , -1]))))
  dat %>%
    kable(digits = comp_stats_digits[i], caption = names(obj_stats[i]), 
          format.args = list(big.mark = ",")) %>%
    skable() %>%
    row_spec(nrow(dat), bold = TRUE) %>%
    print()
  }
}

## ----stats_tables_month_fun---------------------------------------------------------
stats_tables_month <- function(df, obj_col, obj_stats = comp_stats) {
  for(i in seq_along(obj_stats)) {
    dat <- df %>% 
      ungroup() %>%
      mutate(purrr::map_dbl(df[[obj_col]], obj_stats[i]))
    names(dat)[ncol(dat)] <- obj_stats[i]
    dat <- dat %>% 
      mutate(station = as.character(station)) %>%
      pivot_wider(id_cols = c("station", "month_abb"), names_from = "product", 
                  values_from = obj_stats[i]) %>%
      data.frame()
  dat[nrow(dat) + 1, ] <- c(list("(All)", "(All)"), 
                            as.list(as.numeric(colMeans(dat[ , -c(1, 2)], na.rm = TRUE))))
  dat %>%
    kable(digits = comp_stats_digits[i], caption = names(obj_stats[i]), 
          format.args = list(big.mark = ",")) %>%
    skable() %>%
    row_spec(nrow(dat), bold = TRUE) %>%
    print()
  }
}

# TODO Remove IMERG from graphs and tables as no longer using in

## ----yearly_total_plots-------------------------------------------------------------
hus_year <- gof_pr %>%
  unnest(cols = c(data)) %>%
  dplyr::select(product, station, syear, total_rain__station, total_rain) %>%  
  mutate(total_rain__station = ifelse(is.na(total_rain), NA, total_rain__station),
         total_rain = ifelse(is.na(total_rain__station), NA, total_rain)) %>%
  pivot_longer(cols = c("total_rain__station", "total_rain"), 
               names_to = "Source", values_to = "total_rain") %>%
  ungroup() %>%
  mutate(Source = forcats::fct_recode(Source, 
                                      Gauge = "total_rain__station", Estimate = "total_rain"),
         Source = factor(Source, levels = c("Gauge", "Estimate")),
         product = toupper(substr(product, 1, nchar(product) - 5)),
         product = factor(product, levels = c("CHIRPS", "ERA5")))

mean_df <- hus_year %>%
  group_by(product, Source) %>%
  summarise(m = mean(total_rain, na.rm = TRUE))
g <- ggplot(hus_year, aes(x = syear, y = total_rain, colour = Source)) +
  geom_line() +
  geom_point() +
  geom_hline(data = mean_df, aes(yintercept = m, colour = Source)) +
  scale_x_continuous(breaks = seq(1980, 2020, 5), limits = c(1978, 2022), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = "", y = "Total Rainfall (mm/year)") +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 12),          
    legend.title = element_text(size = 14),         
    legend.key.size = unit(2, "lines"),
    strip.text.y = element_text(margin = margin(r = 1, l = 1)), 
    panel.spacing = unit(0.3, "lines"),
    #axis.text.x = element_text(angle = 45),
    axis.text = element_text(face = "bold", size = 12, family = "Helvetica"), 
    text = element_text(face = "bold", size = 12, family = "Helvetica")
  ) +
  facet_wrap(~product) +
  scale_color_manual(
    values = c("Gauge" = "black", "Estimate" = "dodgerblue2"),  
    name = "Source" 
  ) +
  #ggtitle("Gauge vs estimate Annual total rainfall at Husbands") +
  facet_wrap(~product, ncol = 2)
g
##! Figure hus_total (Section 4.3.2.1)
ggsave(here("results", paste0("Fig5.jpeg")),
       width = 12, height = 6)

## ----yearly_n_obs---------------------------------------------------
# dat <- gof_syear %>% 
#   select(station, product, n) %>%
#   pivot_wider(id_cols = c(station), names_from = product, values_from = n)
# 
# names(dat)[endsWith(names(dat), "rain")] <- 
#   substr(names(dat)[endsWith(names(dat), "rain")], 
#          1, nchar(names(dat)[endsWith(names(dat), "rain")]) - 5)
# 
# dat %>%
#   kable(caption = "Number of years compared") %>%
#   skable()

## ----yearly_total_rain_tables---------------------------------------
df <- gof_syear %>%
  ungroup() %>%
  mutate(r = purrr::map_dbl(gof__total_rain, "r"),
         pbias = purrr::map_dbl(gof__total_rain, "PBIAS %"),
         MAE = purrr::map_dbl(gof__total_rain, "MAE")) %>%
  pivot_longer(cols = c("r", "pbias", "MAE"), names_to = "metric") %>%
  mutate(metric = factor(metric, levels = c("r", "pbias", "MAE"))) %>%
  dplyr::select(station, product, metric, value) %>%
  mutate(value = ifelse(metric == "r", format(value, digits = 2), 
                        round(value))) %>%
  pivot_wider(names_from = "product", values_from = "value") %>%
  arrange(metric, station)
names(df)[3:ncol(df)] <- toupper(substr(names(df)[3:ncol(df)], 1, 
                                        nchar(names(df)[3:ncol(df)]) - 5))

##! Table hus_total_metrics (Section 4.3.2.1)
df %>% 
  dplyr::select(metric, CHIRPS, ERA5) %>%
  kable() %>%
  collapse_rows(columns = 1) %>%
  skable()

## ----monthly_calcs------------------------------------------------------------------
by_month <- zm_long %>%
  group_by(station, syear, month_abb, product) %>%
  summarise(total_rain = sum(naif_nmin(rain, 25)), 
            n_rain0p85 = sum(naif_nmin(rain, 25) > 0.85), 
            max_rain = max(naif_nmin(rain, 20)), 
            mean_rain0p85 = total_rain/n_rain0p85, 
            n_na = sum(is.na(rain))
            )

by_month_st <- by_month %>%
  pivot_wider(id_cols = c(station, syear, month_abb), 
              names_from = product, values_from = total_rain:n_na, names_sep = "__") %>%
  pivot_longer(cols = -c(station, syear, month_abb, ends_with("station")), 
               names_to = c(".value", "product"), names_sep = "__")

gof_each_month <- by_month_st %>%
  group_by(station, product, month_abb) %>%
  nest() %>%
  mutate(n = purrr::map_int(data, 
                        ~sum(!is.na(.$total_rain__station) & !is.na(.$total_rain))),
         gof__total_rain = purrr::map(data, dgof, "total_rain", "total_rain__station", 
                                      na.rm = TRUE),
         gof__n_rain0p85 = purrr::map(data, dgof, "n_rain0p85", "n_rain0p85__station", 
                                      na.rm = TRUE)
         )

gof_month <- by_month_st %>%
  group_by(station, product) %>%
  nest() %>%
  mutate(n = purrr::map_int(data, 
                        ~sum(!is.na(.$total_rain__station) & !is.na(.$total_rain))),
         gof__total_rain = purrr::map(data, dgof, "total_rain", "total_rain__station", 
                                      na.rm = TRUE),
         gof__n_rain0p85 = purrr::map(data, dgof, "n_rain0p85", "n_rain0p85__station", 
                                      na.rm = TRUE),
         # gof__mean_rain = purrr::map(data, dgof, "mean_rain", "mean_rain__station", 
         #                             na.rm = TRUE)
         )

gof_pr_month <- gof_month %>%
  unnest(cols = data) %>% 
  group_by(product) %>% 
  nest()

## ----monthly_plots_fun--------------------------------------------------------------
monthly_plots <- function(df, stat_pr, stat_st, product_name) {
  vals_relace <- c(product_name, "station")
  names(vals_relace) <- c(stat_pr, stat_st)
  dat <- df %>% 
    filter(!is.na(.data[[stat_pr]]) & !is.na(.data[[stat_st]]) & 
             month_abb %in% month.abb) %>%
    pivot_longer(cols = c(stat_pr, stat_st), names_to = "product", values_to = stat_pr) %>%
    mutate(product = recode(product, !!!vals_relace))

  mean_df <- dat %>% 
    group_by(station, product, month_abb) %>% 
    summarise(m = mean(.data[[stat_pr]], na.rm = TRUE))

  g <- ggplot(dat, aes(x = syear, y = .data[[stat_pr]], colour = product)) +
    geom_line() +
    geom_hline(data = mean_df, aes(yintercept = m, colour = product)) +
    scale_x_continuous(limits = c(1979, 2020)) +
    ggtitle(paste(stat_pr, ":", product_name, "vs", "station")) +
    facet_wrap(~month_abb)
  g
}

## ----monthly_plots_total_rain_chirps------------------------------------------------
hus_month <- gof_pr_month %>%
  unnest(cols = c(data)) %>%
  dplyr::select(product, station, syear, month_abb, total_rain__station, total_rain) %>%  
  mutate(total_rain__station = ifelse(is.na(total_rain), NA, total_rain__station),
         total_rain = ifelse(is.na(total_rain__station), NA, total_rain)) %>%
  pivot_longer(cols = c("total_rain__station", "total_rain"), 
               names_to = "Source", values_to = "total_rain") %>%
  ungroup() %>%
  mutate(Source = forcats::fct_recode(Source, 
                                      Gauge = "total_rain__station", Estimate = "total_rain"),
         Source = factor(Source, levels = c("Gauge", "Estimate")),
         product = toupper(substr(product, 1, nchar(product) - 5)),
         product = factor(product, levels = c("CHIRPS", "ERA5", "IMERG")))

mean_df <- hus_month %>%
  group_by(product, Source, month_abb) %>%
  summarise(m = mean(total_rain, na.rm = TRUE))

g <- ggplot(hus_month %>% filter(product == "CHIRPS"), 
            aes(x = syear, y = total_rain, colour = Source)) +
  geom_line() +
  geom_point() +
  geom_hline(data = mean_df %>% filter(product == "CHIRPS"), aes(yintercept = m, colour = Source)) +
  scale_x_continuous(breaks = seq(1980, 2020, 5), limits = c(1978, 2022), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, NA)) +
  scale_colour_discrete(labels = c("Gauge", "CHIRPS")) +
  labs(x = "Year", y = "Total Rainfall (mm/month)") +
  #ggtitle("Gauge vs CHIRPS Monthly total rainfall at Husbands") +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 12),          
    legend.title = element_text(size = 14),         
    legend.key.size = unit(2, "lines"),
    strip.text.y = element_text(margin = margin(r = 1, l = 1)), 
    panel.spacing = unit(0.3, "lines"),
    axis.text.x = element_text(angle = 45),
    axis.text = element_text(face = "bold", size = 12, family = "Helvetica"), 
    text = element_text(face = "bold", size = 12, family = "Helvetica")
  ) +
  facet_wrap(~ month_abb)+
  scale_color_manual(
    values = c("Gauge" = "black", "Estimate" = "dodgerblue2"),  
    name = "Source" 
  )
g


##! Figure hus_total_month (Section 4.3.2.1)
ggsave(here("results", paste0("Fig6.jpeg")),
       width = 12, height = 6)
