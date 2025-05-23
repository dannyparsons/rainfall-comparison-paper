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

zm <- readRDS(here("data", "station", "cleaned", "westafrica_gridded.RDS")) %>%
  select(station, date, rain, month, year, chirps_rain, era5_rain, tamsat_rain) %>%
  filter(station == "Sadore")
ghana_df <- readRDS("data", "station", "cleaned", "ghana_df.RDS")

zm <- bind_rows(zm, ghana_df) %>%
  filter(year >= 1983 & year <= 2022) %>%
  mutate(chirps_rain = ifelse(chirps_rain < 0, 0, chirps_rain),
         era5_rain = ifelse(era5_rain < 0, 0, era5_rain),
         tamsat_rain = ifelse(tamsat_rain < 0, 0, tamsat_rain)
         )
# 1 Jan = 1
s_doy_start <- 1
zm <- zm %>% 
  mutate(doy = yday_366(date),
         s_doy = (doy - s_doy_start + 1) %% 366,
         s_doy = ifelse(s_doy == 0, 366, s_doy),
         syear = year,
         syear = ifelse(s_doy > (366 - s_doy_start + 1), syear - 1, syear),
         month = factor(month),
         month_abb = month(date, label = TRUE),
         
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
stations <- c("Sadore", "Wa", "Tamale", "Saltpond")
products <- levels(zm_long$product)
products <- products[-1]
names(products) <- substr(products, 1, nchar(products) - 5)

metadata_station <- readRDS(here("data", "station", "processed", "stations_metadata.RDS"))
metadata_zm <- metadata_station %>%
  filter(country %in% c("Niger")) %>% 
  select(country, station, latitude, longitude)
rm(metadata_station) 

gh_metadata <- read.csv(here("data", "station", "processed", "ghana_metadata.csv")) %>% 
  mutate(country = "Ghana") %>%
  select(country, station, latitude, longitude)

metadata_zm <- bind_rows(metadata_zm, gh_metadata)
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

min_max_year <- zm_long %>%
  select(station, first_date, last_date) %>%
  mutate(min_year = year(first_date),
         max_year = year(last_date)) %>%
  unique() %>%
  select(station, min_year, max_year)

skable <- function(kable_input) {
  kable_input %>% kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                                full_width = FALSE)
}

metadata_zm %>% 
  left_join(min_max_year, by="station") %>%
  select(station, latitude, longitude, min_year, max_year) %>%
  arrange(-latitude) %>%
  kable() %>%
  skable()

zm %>% 
  group_by(station) %>% 
  summarise(prop_n = 100* (1 - naflex::na_prop(rain))) %>%
  kable(digits = 1) %>%
  skable()

sf_wa <- ne_countries(returnclass = "sf")
ggplot(sf_wa) + 
  geom_sf() +
  geom_point(data = metadata_zm, aes(x = longitude, y = latitude)) +
  geom_text_repel(data = metadata_zm, aes(x = longitude, y = latitude, label = station)) +
  coord_sf(xlim = c(-8, 15), ylim = c(0, 25), expand = FALSE)

rm(sf_wa)


ggplot(zm, aes(x = date, y = station, fill = !is.na(rain))) +
  geom_tile() +
  geom_hline(yintercept = seq(0.5, by = 1, 
                              length.out = length(unique(zm$station)) + 1)) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y")


naif_nmin <- function(x, n_min) {
  if(length(na.omit(x)) > 0 && sum(!is.na(x)) >= n_min) {
    na.omit(x)
  } else NA
}

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

by_syear <- zm_long %>%
  group_by(station, syear, product) %>%
  summarise(total_rain = sum(naif_nmin(rain, 355)), 
            n_rain = sum(naif_nmin(rain, 355) > 0.1), 
            max_rain = max(naif_nmin(rain, 350)), 
            mean_rain = total_rain/n_rain, 
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
         gof__n_rain = purrr::map(data, dgof, "n_rain", "n_rain__station", na.rm = TRUE),
         gof__mean_rain = purrr::map(data, dgof, "mean_rain", "mean_rain__station",
                                     na.rm = TRUE)
         )

gof_pr <- gof_syear %>% 
  unnest(cols = data) %>% 
  group_by(product) %>% 
  nest()

yearly_plots <- function(df, gof_col, stat_pr, stat_st, product_name) {
  max_y <- max(c(df[[stat_pr]], df[[stat_st]]), na.rm = TRUE)
  vals_relace <- c(product_name, "station")
  names(vals_relace) <- c(stat_pr, stat_st)
  dat <- df %>% 
    #filter(!is.na(.data[[stat_pr]]) & !is.na(.data[[stat_st]])) %>%
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
    geom_point() +
    #geom_hline(data = mean_df, aes(yintercept = m, colour = product)) +
    scale_x_continuous(limits = c(1979, 2020)) +
    # n
    geom_text(data = dat, aes(label = paste("n", n)), size = 3,
              x = 1979, y = max_y, na.rm = TRUE, 
              inherit.aes = FALSE) +
    # bias
    geom_text(data = dat, aes(label = paste("bias", signif(ME, 2))), 
              size = 3, x = 1979 + 5, y = max_y, na.rm = TRUE, 
              inherit.aes = FALSE) +
    # cor
    geom_text(aes(label = paste("cor", round(r, 2))), 
              size = 3, x = 1979 + 12, y = max_y, na.rm = TRUE, 
              inherit.aes = FALSE) +
    # rSD
    geom_text(aes(label = paste("rSD", round(rSD, 2))), 
              size = 3, x = 1979 + 19, y = max_y, na.rm = TRUE, 
              inherit.aes = FALSE) +
    ggtitle(paste(stat_pr, ":", product_name, "vs", "station")) +
    facet_wrap(~station)
  g
}

stats_tables <- function(df, obj_col, obj_stats = comp_stats) {
  for(i in seq_along(obj_stats)) {
    dat <- df %>% 
      ungroup() %>%
      mutate(purrr::map_dbl(df[[obj_col]], obj_stats[i]))
    names(dat)[ncol(dat)] <- obj_stats[i]
    dat <- dat %>% 
      mutate(station = as.character(station))
    dat <- dat %>%
      pivot_wider(id_cols = "station", names_from = "product", 
                  values_from = tidyselect::all_of(obj_stats[[i]])) %>%
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

p_syear_total <- gof_pr %>%
  mutate(p = purrr::map(data, yearly_plots, gof_col = "gof__total_rain", 
                        stat_pr = "total_rain", stat_st = "total_rain__station", 
                        product_name = product),
         paths = here("results", "westafrica", 
                      paste0("westafrica_", "syear_", "total_rain", "_", product, ".png")))

walk2(p_syear_total$paths, p_syear_total$p, ggsave, width = 12, height = 6)
walk(p_syear_total$p, print)

dat <- gof_syear %>% 
  select(station, product, n) %>%
  pivot_wider(id_cols = c(station), names_from = product, values_from = n)

names(dat)[endsWith(names(dat), "rain")] <- 
  substr(names(dat)[endsWith(names(dat), "rain")], 
         1, nchar(names(dat)[endsWith(names(dat), "rain")]) - 5)
dat %>%
  kable(caption = "Number of years compared") %>%
  skable()

stats_tables(gof_syear, "gof__total_rain")

p_syear_n <- gof_pr %>%
  mutate(p = purrr::map(data, yearly_plots, gof_col = "gof__n_rain", 
                        stat_pr = "n_rain", stat_st = "n_rain__station", 
                        product_name = product),
         paths = here("results", "westafrica", 
                      paste0("westafrica_", "syear_", "n_rain", "_", product, ".png")))

walk2(p_syear_n$paths, p_syear_n$p, ggsave, width = 12, height = 6)
walk(p_syear_n$p, print)

stats_tables(gof_syear, "gof__n_rain")

by_month <- zm_long %>%
  group_by(station, syear, month_abb, product) %>%
  filter(month %in% c(3:10)) %>%
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

monthly_plots <- function(df, stat_pr, stat_st, product_name) {
  vals_relace <- c(product_name, "station")
  names(vals_relace) <- c(stat_pr, stat_st)
  dat <- df %>% 
    filter(month_abb %in% month.abb[c(3:10)]) %>%
    pivot_longer(cols = c(stat_pr, stat_st), names_to = "product", values_to = stat_pr) %>%
    mutate(product = recode(product, !!!vals_relace))

  mean_df <- dat %>% 
    group_by(station, product, month_abb) %>% 
    summarise(m = mean(.data[[stat_pr]], na.rm = TRUE))

  g <- ggplot(dat, aes(x = syear, y = .data[[stat_pr]], colour = product)) +
    geom_line() +
    geom_point(size = 1) +
    #geom_hline(data = mean_df, aes(yintercept = m, colour = product)) +
    scale_x_continuous(limits = c(1979, 2020)) +
    ggtitle(paste(stat_pr, ":", product_name, "vs", "station")) +
    facet_grid(station~month_abb, scales = "free_y") +
    theme(strip.text.y = element_text(size = 6))
  g
}

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

p_month_total <- gof_pr_month %>%
  mutate(p = purrr::map(data, monthly_plots, stat_pr = "total_rain", 
                        stat_st = "total_rain__station", product_name = product),
         paths = here("results", "westafrica", 
                      paste0("westafrica_", "month_", "total_rain", "_", product, ".png")))

walk2(p_month_total$paths, p_month_total$p, ggsave, width = 12, height = 6)
walk(p_month_total$p, print)

stats_tables(gof_month, "gof__total_rain")

p_month_total <- gof_pr_month %>%
  mutate(p = purrr::map(data, monthly_plots, stat_pr = "n_rain0p85", 
                        stat_st = "n_rain0p85__station", product_name = product),
         paths = here("results", "westafrica", 
                      paste0("westafrica_", "month_", "n_rain0p85", "_", product, ".png")))

walk2(p_month_total$paths, p_month_total$p, ggsave, width = 12, height = 6)
walk(p_month_total$p, print)

stats_tables(gof_month, "gof__n_rain0p85")

zambia_markov <- zm_long_st %>% 
  filter(!is.na(rain) & !is.na(pr_rain)) %>%
  mutate(rainday1 = rain > 0.85,
         pr_rainday1 = pr_rain > 0.85,
         pr_rainday3 = pr_rain > 3,
         pr_rainday4 = pr_rain > 4,
         pr_rainday5 = pr_rain > 5,
         pr_rainday2max = pr_rain > 1 & pr_rain <= 2,
         pr_rainday3max = pr_rain > 1 & pr_rain <= 3,
         pr_rainday4max = pr_rain > 1 & pr_rain <= 4,
         pr_rainday5max = pr_rain > 1 & pr_rain <= 5)

f_zero_order_station <- rainday1 ~ (cos(s_doy * 1 * 2 * pi/366) +
                                    sin(s_doy * 1 * 2 * pi/366) + 
                                    cos(s_doy * 2 * 2 * pi/366) + 
                                    sin(s_doy * 2 * 2 * pi/366) +
                                    cos(s_doy * 3 * 2 * pi/366) +
                                    sin(s_doy * 3 * 2 * pi/366))
f_zero_order_product <- update.formula(f_zero_order_station, pr_rainday1 ~ .)

predict_stack_lst <- list()
for(s in seq_along(stations)) {
  predict_df <- data.frame(station = stations[s], s_doy = 1:366,
                           s_doy_date = as.Date(1:366, origin = as.Date("1999/12/31")))
  dat <- zambia_markov %>%
    filter(station == stations[s])
  for(i in seq_along(products)) {
    dat_prod <- dat %>%
      filter(product %in% c("station", products[i])) %>% 
      filter(!is.na(rain) & !is.na(pr_rain))
    zero_order_station <- glm(f_zero_order_station, data = dat_prod, family = binomial)
    zero_order_product <- glm(f_zero_order_product, data = dat_prod, family = binomial)
    #print(anova(zero_order_station, test="Chisq"))
    predict_df[[paste0("station", "_", products[i])]] <- predict(zero_order_station, 
                                                                 newdata = predict_df, 
                                                                 type = "response")
    predict_df[[products[i]]] <- predict(zero_order_product, newdata = predict_df, 
                                         type = "response")
    
    f_zero_order_product_3thres <- update.formula(f_zero_order_station, pr_rainday3 ~ .)
    f_zero_order_product_4thres <- update.formula(f_zero_order_station, pr_rainday4 ~ .)
    f_zero_order_product_5thres <- update.formula(f_zero_order_station, pr_rainday5 ~ .)
    f_zero_order_product_2max <- update.formula(f_zero_order_station, pr_rainday2max ~ .)
    f_zero_order_product_3max <- update.formula(f_zero_order_station, pr_rainday3max ~ .)
    f_zero_order_product_4max <- update.formula(f_zero_order_station, pr_rainday4max ~ .)
    f_zero_order_product_5max <- update.formula(f_zero_order_station, pr_rainday5max ~ .)
    fms <- list(f_zero_order_product_2max, f_zero_order_product_3max, 
                f_zero_order_product_4max, f_zero_order_product_5max)
    fms_thres <- list(f_zero_order_product_3thres, f_zero_order_product_4thres,
                      f_zero_order_product_5thres)
    mds <- list()
    for(j in 2:(2 + length(fms) - 1)) {
      zero_order <- glm(fms[[j - 1]], data = dat_prod, family = binomial)
      predict_df[[paste0(products[i], "_", j, "max")]] <- predict(zero_order, 
                                                                  newdata = predict_df, 
                                                                  type = "response")
    }
    for(j in 3:(3 + length(fms_thres) - 1)) {
      zero_order <- glm(fms_thres[[j - 2]], data = dat_prod, family = binomial)
      predict_df[[paste0(products[i], "_", j, "thres")]] <- predict(zero_order, 
                                                                    newdata = predict_df, 
                                                                    type = "response")
    }
  }
  
  predict_stack <- predict_df %>% melt(id.vars = c("station", "s_doy", "s_doy_date"), 
                                       variable.name = "product", value.name = "prob")

  predict_stack$product <- as.character(predict_stack$product)
  predict_stack$product2 <- predict_stack$product
  predict_stack$type <- "product1"
  
  predict_stack <- predict_stack %>% 
    mutate(type = ifelse(startsWith(product, "station"), "station1", type),
           type = ifelse(endsWith(product, "max"), 
                         substr(product, nchar(product) - 3, nchar(product)), type),
           type = ifelse(endsWith(product, "thres"), 
                         substr(product, nchar(product) - 5, nchar(product)), type),
           product2 = ifelse(startsWith(product, "station"), 
                             substr(product, 9, nchar(product)), product2),
           product2 = ifelse(endsWith(product, "max"), 
                             substr(product, 1, nchar(product) - 5), product2),
           product2 = ifelse(endsWith(product, "thres"), 
                             substr(product, 1, nchar(product) - 7), product2)
           )
  
  predict_stack$type <- factor(predict_stack$type, levels = unique(predict_stack$type))
  predict_stack_lst[[length(predict_stack_lst) + 1]] <- predict_stack
  # Plot small amounts
  # g <- ggplot(predict_stack, aes(x = s_doy, y = prob, colour = type)) +
  #   geom_line() +
  #   facet_wrap(~product2) +
  #   scale_color_manual(values = c("black", c25[1:7])) +
  #   ggtitle(paste("Chance of rain:", stations[s]))
  # ggsave(here("results", "zambia", paste0("zambia_", "markov_zero", stations[s], ".png")), 
  #        plot = g, width = 12, height = 6)
}
predict_stack_all <- bind_rows(predict_stack_lst)

types <- c("Gauge > 0.85mm", " > 0.85mm",
           " > 2mm", " > 3mm",
           " > 4mm", " > 5mm")
for(s in seq_along(products)) {
  curve_labs <- paste0(c("", rep(toupper(names(products)[s]), 5)), types)
  dat <- predict_stack_all %>% 
    filter(product2 == products[s] & !grepl("max", type))
  g <- ggplot(dat, aes(x = s_doy_date, y = prob, 
                       colour = type, size = type, linetype = type)) +
    geom_line() +
    facet_wrap(~station) +
    scale_color_manual(values = c("black", viridis(4, end = 0.8)), name = "Curve", labels = curve_labs) +
    scale_size_manual(values = c(1.3, rep(0.8, 4)), name = "Curve", labels = curve_labs) +
    scale_linetype_manual(values = c("solid", rep("longdash", 4)), name = "Curve", labels = curve_labs) +
    scale_x_date(date_breaks = "2 months", date_labels = "%b", name = "Date") +
    scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1), name = "Rain day frequency (rain days/day)") +
    theme(aspect.ratio = 0.5) +
    ggtitle(paste0("Gauge and ", toupper(names(products)[s]), " ", 
                  "estimated rain day frequency (rain days/day) at various thresholds"))
    ggtitle(paste("Chance of rain:", products[s]))
  print(g)
  ggsave(here("results", "westafrica", 
              paste0("westafrica_", "markov_zero", "_product_", names(products)[s], ".png")),
         plot = g, width = 12, height = 6)
}

zambia_markov <- zm_long_st %>% 
  filter(!is.na(rain) & !is.na(pr_rain)) %>%
  mutate(rainday1 = rain > 0.85,
         pr_rainday1 = pr_rain > 0.85,
         pr_rainday2 = pr_rain > 2,
         pr_rainday3 = pr_rain > 3,
         pr_rainday4 = pr_rain > 4,
         pr_rainday5 = pr_rain > 5,
         lag1 = factor(ifelse(dplyr::lag(rainday1, 1), "w", "d")),
         pr_lag1 = factor(ifelse(dplyr::lag(pr_rainday1, 1), "w", "d")),
         pr_lag2 = factor(ifelse(dplyr::lag(pr_rainday2, 1), "w", "d")),
         pr_lag3 = factor(ifelse(dplyr::lag(pr_rainday3, 1), "w", "d")),
         pr_lag4 = factor(ifelse(dplyr::lag(pr_rainday4, 1), "w", "d")),
         pr_lag5 = factor(ifelse(dplyr::lag(pr_rainday5, 1), "w", "d"))
         )

f_first_order_station <- rainday1 ~ ((cos(s_doy * 1 * 2 * pi/366) +
                                       sin(s_doy * 1 * 2 * pi/366) +
                                       cos(s_doy * 2 * 2 * pi/366) +
                                       sin(s_doy * 2 * 2 * pi/366) +
                                       cos(s_doy * 3 * 2 * pi/366) +
                                       sin(s_doy * 3 * 2 * pi/366)) * lag1)

f_first_order_product <- pr_rainday1 ~ ((cos(s_doy * 1 * 2 * pi/366) +
                                       sin(s_doy * 1 * 2 * pi/366) +
                                       cos(s_doy * 2 * 2 * pi/366) +
                                       sin(s_doy * 2 * 2 * pi/366) +
                                       cos(s_doy * 3 * 2 * pi/366) +
                                       sin(s_doy * 3 * 2 * pi/366)) * pr_lag1)

predict_stack_lst <- list()
for(s in seq_along(stations)) {
  predict_df <- expand.grid(station = stations[s], s_doy = 1:366, wd = c("d", "w"))
  predict_df$s_doy_date = as.Date(predict_df$s_doy, origin = as.Date("1999/12/31"))
  dat <- zambia_markov %>%
    filter(station == stations[s])
  for(i in seq_along(products)) {
    dat_prod <- dat %>%
      filter(product %in% c("station", products[i])) %>% 
      filter(!is.na(rain) & !is.na(pr_rain))
    first_order_station <- glm(f_first_order_station, data = dat_prod, family = binomial)
    first_order_product <- glm(f_first_order_product, data = dat_prod, family = binomial)
    # print(anova(first_order_station, test = "Chisq"))
    names(predict_df)[3] <- "lag1"
    predict_df[[paste0("station", "_", products[i])]] <- predict(first_order_station, 
                                                                 newdata = predict_df, 
                                                                 type = "response")
    names(predict_df)[3] <- "pr_lag1"
    predict_df[[products[i]]] <- predict(first_order_product, newdata = predict_df, 
                                         type = "response")

    f_first_order_product_2thres <- pr_rainday2 ~ ((cos(s_doy * 1 * 2 * pi/366) +
                                       sin(s_doy * 1 * 2 * pi/366) +
                                       cos(s_doy * 2 * 2 * pi/366) +
                                       sin(s_doy * 2 * 2 * pi/366) +
                                       cos(s_doy * 3 * 2 * pi/366) +
                                       sin(s_doy * 3 * 2 * pi/366)) * pr_lag2)
    f_first_order_product_3thres <- pr_rainday3 ~ ((cos(s_doy * 1 * 2 * pi/366) +
                                       sin(s_doy * 1 * 2 * pi/366) +
                                       cos(s_doy * 2 * 2 * pi/366) +
                                       sin(s_doy * 2 * 2 * pi/366) +
                                       cos(s_doy * 3 * 2 * pi/366) +
                                       sin(s_doy * 3 * 2 * pi/366) * pr_lag3))
    f_first_order_product_4thres <- pr_rainday4 ~ ((cos(s_doy * 1 * 2 * pi/366) +
                                       sin(s_doy * 1 * 2 * pi/366) +
                                       cos(s_doy * 2 * 2 * pi/366) +
                                       sin(s_doy * 2 * 2 * pi/366) +
                                       cos(s_doy * 3 * 2 * pi/366) +
                                       sin(s_doy * 3 * 2 * pi/366) * pr_lag4))
    f_first_order_product_5thres <- pr_rainday5 ~ ((cos(s_doy * 1 * 2 * pi/366) +
                                       sin(s_doy * 1 * 2 * pi/366) +
                                       cos(s_doy * 2 * 2 * pi/366) +
                                       sin(s_doy * 2 * 2 * pi/366) +
                                       cos(s_doy * 3 * 2 * pi/366) +
                                       sin(s_doy * 3 * 2 * pi/366) * pr_lag5))
    fms_thres <- list(f_first_order_product_2thres, f_first_order_product_3thres,
                      f_first_order_product_4thres, f_first_order_product_5thres)
    mds <- list()
    for(j in seq_along(fms_thres)) {
      first_order <- glm(fms_thres[[j]], data = dat_prod, family = binomial)
      names(predict_df)[3] <- paste0("pr_lag", j + 1)
      predict_df[[paste0(products[i], "_", j + 1, "thres")]] <- predict(first_order,
                                                                    newdata = predict_df,
                                                                    type = "response")
    }
  }
  names(predict_df)[3] <- "lag1"
  predict_stack <- predict_df %>% melt(id.vars = c("station", "s_doy", "s_doy_date", "lag1"), 
                                       variable.name = "product", value.name = "prob")

  predict_stack$product <- as.character(predict_stack$product)
  predict_stack$product2 <- predict_stack$product
  predict_stack$type <- "product1"

  predict_stack <- predict_stack %>%
    mutate(type = ifelse(startsWith(product, "station"), "station1", type),
           type = ifelse(endsWith(product, "max"),
                         substr(product, nchar(product) - 3, nchar(product)), type),
           type = ifelse(endsWith(product, "thres"),
                         substr(product, nchar(product) - 5, nchar(product)), type),
           product2 = ifelse(startsWith(product, "station"),
                             substr(product, 9, nchar(product)), product2),
           product2 = ifelse(endsWith(product, "max"),
                             substr(product, 1, nchar(product) - 5), product2),
           product2 = ifelse(endsWith(product, "thres"),
                             substr(product, 1, nchar(product) - 7), product2)
           )

  predict_stack$type <- factor(predict_stack$type, levels = unique(predict_stack$type))
  predict_stack_lst[[length(predict_stack_lst) + 1]] <- predict_stack
}
predict_stack_all <- bind_rows(predict_stack_lst)
predict_stack_all$station <- factor(predict_stack_all$station, levels = stations)

thres <- c(0.85, 2, 3, 4, 5)
names_thres <- c("product1", paste0(2:5, "thres"))
for (i in seq_along(thres)) {
  for(s in seq_along(products)) {
    curve_labs <- paste0(c("", rep(toupper(names(products)[s]), 5)), types)
    dat <- predict_stack_all %>%
      filter(product2 == products[s] & type %in% c("station1", names_thres[i]))
    g <- ggplot(dat, aes(x = s_doy_date, y = prob,
                         colour = lag1, linetype = type)) +
      geom_line() +
      facet_wrap(~station) +
      scale_x_date(date_breaks = "2 months", date_labels = "%b", name = "Date") +
      scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1), name = "Rain day frequency (rain days/day)") +
      scale_linetype_manual(values = c("solid", "longdash")) +
      theme(aspect.ratio = 1) +
      ggtitle(paste0("Gauge and ", toupper(names(products)[s]), " ",
                     "estimated rain day frequency (rain days/day) at various thresholds"))
    ggsave(here("results", "westafrica",
                paste0("westafrica_", "markov_first", "_product_", thres[i], "_", names(products)[s], ".png")),
           plot = g, width = 12, height = 6)
  }
}
