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
  select(product, station, syear, total_rain__station, total_rain) %>%  
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

mean_df <- hus_year %>%
  group_by(product, Source) %>%
  summarise(m = mean(total_rain, na.rm = TRUE))

g <- ggplot(hus_year, aes(x = syear, y = total_rain, colour = Source)) +
  geom_line() +
  geom_point() +
  geom_hline(data = mean_df, aes(yintercept = m, colour = Source)) +
  scale_x_continuous(breaks = seq(1980, 2020, 5), limits = c(1978, 2020), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = "Year", y = "Total Rainfall (mm/year)") +
  ggtitle("Gauge vs estimate Annual total rainfall at Husbands") +
  facet_wrap(~product, ncol = 2)
g
##! Figure hus_total (Section 4.3.2.1)
ggsave(here("results", paste0("barbados_", "syear_", "total_rain", "_", "all", ".png")),
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
  select(station, product, metric, value) %>%
  mutate(value = ifelse(metric == "r", format(value, digits = 2), 
                        round(value))) %>%
  pivot_wider(names_from = "product", values_from = "value") %>%
  arrange(metric, station)
names(df)[3:ncol(df)] <- toupper(substr(names(df)[3:ncol(df)], 1, 
                                        nchar(names(df)[3:ncol(df)]) - 5))

##! Table hus_total_metrics (Section 4.3.2.1)
df %>% 
  dplyr::select(metric, CHIRPS, ERA5, IMERG) %>%
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
  select(product, station, syear, month_abb, total_rain__station, total_rain) %>%  
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
  scale_x_continuous(breaks = seq(1980, 2020, 5), limits = c(1978, 2020), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, NA)) +
  scale_colour_discrete(labels = c("Gauge", "CHIRPS")) +
  labs(x = "Year", y = "Total Rainfall (mm/month)") +
  ggtitle("Gauge vs CHIRPS Monthly total rainfall at Husbands") +
  facet_wrap(~ month_abb)

##! Figure hus_total_month (Section 4.3.2.1)
ggsave(here("results", paste0("barbados_", "monthly_", "total_rain", "_", "chirps", ".png")),
       width = 12, height = 6)

## Next

## ----all_months_tables_total_rain, results="asis"-----------------------------------
stats_tables(gof_month, "gof__total_rain")

#' 
## ----months_tables_total_rain, results="asis"---------------------------------------
stats_tables_month(gof_each_month, "gof__total_rain")

#' 
#' ### Comparison statistics for number of monthly raindays
#' - Similar to yearly values, good but relatively lower correlation and same over/under estimation
#' - ERA5 has a better correlation and fairly consistent each month
#' - CHIRPS has 0 correlation Dec-Feb, IMERG in Mar-Apr
#' - ERA5 and IMERG have a reasonably similar variation
#' 
## ----monthly_plots_n_rain-----------------------------------------------------------
p_month_total <- gof_pr_month %>%
  mutate(p = purrr::map(data, monthly_plots, stat_pr = "n_rain0p85", 
                        stat_st = "n_rain0p85__station", product_name = product),
         paths = here("results", "barbados", 
                      paste0("barbados_", "month_", "n_rain0p85", "_", product, ".png")))

walk2(p_month_total$paths, p_month_total$p, ggsave, width = 12, height = 6)
walk(p_month_total$p, print)

#' 
## ----all_month_tables_n_rain0p85, results="asis"------------------------------------
stats_tables(gof_month, "gof__n_rain0p85")

#' 
## ----months_tables_n_rain0p85, results="asis"---------------------------------------
stats_tables_month(gof_each_month, "gof__n_rain0p85")

#' 
#' ## Markov Chain models of the chance of rain
#' 
#' - ERA5 appears to model the shape of chance of rain closest to the station
#' - A simple shift of thresholds doesn't look sensible, in contrast to Zambia results
#' - e.g. a 3mm threshold for ERA5 matches for the most rainy months, whereas between 1 and 1.5 matches best in the drier months
#' - IMERG over and under estimates in different parts of the season
#' - CHIRPS significantly underestimates throughout so a threshold increase doesn't help
#' 
## ----markov_chain_setup-------------------------------------------------------------
zambia_markov <- zm_long_st %>% 
  filter(!is.na(rain) & !is.na(pr_rain)) %>%
  mutate(rainday1 = rain > 0.85,
         pr_rainday1 = pr_rain > 0.85,
         pr_rainday1p5 = pr_rain > 1.5,
         pr_rainday2 = pr_rain > 2,
         pr_rainday3 = pr_rain > 3,
         pr_rainday4 = pr_rain > 4,
         pr_rainday5 = pr_rain > 5,
         pr_rainday2max = pr_rain > 1 & pr_rain <= 2,
         pr_rainday3max = pr_rain > 1 & pr_rain <= 3,
         pr_rainday4max = pr_rain > 1 & pr_rain <= 4,
         pr_rainday5max = pr_rain > 1 & pr_rain <= 5)

threshs <- c(1.5, 2, 3, 4, 5)

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
    
    f_zero_order_product_1p5thres <- update.formula(f_zero_order_station, pr_rainday1p5 ~ .)
    f_zero_order_product_2thres <- update.formula(f_zero_order_station, pr_rainday2 ~ .)
    f_zero_order_product_3thres <- update.formula(f_zero_order_station, pr_rainday3 ~ .)
    f_zero_order_product_4thres <- update.formula(f_zero_order_station, pr_rainday4 ~ .)
    f_zero_order_product_5thres <- update.formula(f_zero_order_station, pr_rainday5 ~ .)
    f_zero_order_product_2max <- update.formula(f_zero_order_station, pr_rainday2max ~ .)
    f_zero_order_product_3max <- update.formula(f_zero_order_station, pr_rainday3max ~ .)
    f_zero_order_product_4max <- update.formula(f_zero_order_station, pr_rainday4max ~ .)
    f_zero_order_product_5max <- update.formula(f_zero_order_station, pr_rainday5max ~ .)
    fms <- list(f_zero_order_product_2max, f_zero_order_product_3max, 
                f_zero_order_product_4max, f_zero_order_product_5max)
    fms_thres <- list(f_zero_order_product_1p5thres, f_zero_order_product_2thres,
                      f_zero_order_product_3thres)
    mds <- list()
    for(j in 2:(2 + length(fms) - 1)) {
      zero_order <- glm(fms[[j - 1]], data = dat_prod, family = binomial)
      predict_df[[paste0(products[i], "_", j, "max")]] <- predict(zero_order, 
                                                                  newdata = predict_df, 
                                                                  type = "response")
    }
    for(j in seq_along(fms_thres)) {
      zero_order <- glm(fms_thres[[j]], data = dat_prod, family = binomial)
      predict_df[[paste0(products[i], "_", threshs[j], "thres")]] <- predict(zero_order, 
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
           type = ifelse(endsWith(product, "1.5thres"), 
                         substr(product, nchar(product) - 7, nchar(product)),
                         ifelse(endsWith(product, "thres"), 
                                substr(product, nchar(product) - 5, nchar(product)), type)),
           product2 = ifelse(startsWith(product, "station"), 
                             substr(product, 9, nchar(product)), product2),
           product2 = ifelse(endsWith(product, "max"), 
                             substr(product, 1, nchar(product) - 5), product2),
           product2 = ifelse(endsWith(product, "1.5thres"), 
                         substr(product, 1, nchar(product) - 9),
                         ifelse(endsWith(product, "thres"), 
                                substr(product, 1, nchar(product) - 7), product2)),
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

#' 
## ----markov_chain_plots_stations----------------------------------------------------
for(s in seq_along(stations)) {
  dat <- predict_stack_all %>% filter(station == stations[s] & !grepl("max", type))
  g <- ggplot(dat, aes(x = s_doy_date, y = prob, colour = type, size = type)) +
    geom_line() +
    facet_wrap(~product2) +
    scale_size_manual(values = c(0.8, rep(0.6, 6))) +
    scale_color_manual(values = c("black", c25[1:6])) +
    scale_x_date(date_breaks = "2 months", date_labels = "%b") +
    ggtitle(paste("Chance of rain:", stations[s]))
  print(g)
  ggsave(here("results", "barbados", 
              paste0("barbados_", "markov_zero", "_station_", stations[s], ".png")),
         plot = g, width = 12, height = 6)
}

#' 
## ----markov_chain_plots_products----------------------------------------------------
for(s in seq_along(products)) {
  dat <- predict_stack_all %>% 
    filter(product2 == products[s] & !grepl("max", type))
  g <- ggplot(dat, aes(x = s_doy_date, y = prob, 
                       colour = type, size = type, linetype = type)) +
    geom_line() +
    facet_wrap(~station) +
    scale_color_manual(values = c("black", viridis(4, end = 0.8))) +
    scale_size_manual(values = c(1.3, rep(0.8, 4))) +
    scale_x_date(date_breaks = "2 months", date_labels = "%b") +
    scale_linetype_manual(values = c("solid", rep("longdash", 4))) +
    ggtitle(paste("Chance of rain:", products[s]))
  print(g)
  ggsave(here("results", "barbados", 
              paste0("barbados_", "markov_zero", "_product_", names(products)[s], ".png")),
         plot = g, width = 12, height = 6)
}

#' 
#' ### Distribution of rainfall amounts
#' 
#' - Look into small rainfall amounts and very large amounts
## -----------------------------------------------------------------------------------
# for(i in c(200, 50, 10, 5)) {
# g <- ggplot(zm_long %>% filter(rain > 0), aes(x = rain, colour = product)) +
#   stat_ecdf(aes(y = 1 - ..y.., size = product), pad = FALSE, geom = "step") +
#   coord_cartesian(xlim = c(NA, i)) +
#   scale_color_manual(values = c("black", c25[1:6])) +
#   scale_size_manual(values = c(1, rep(0.5, 6))) +
#   facet_wrap(~station)
# print(g)
# }

#' 
## -----------------------------------------------------------------------------------
# for(i in c(200, 50, 10, 5)) {
#   g <- ggplot(zm_long %>% filter(rain > 0), 
#               aes(x = rain, y = after_stat(density), colour = product)) +
#     geom_freqpoly(aes(size = product), binwidth = 0.2, center = 0.1) +
#     scale_color_manual(values = c("black", c25[1:6])) +
#     scale_size_manual(values = c(1, rep(0.5, 6))) +
#     coord_cartesian(xlim = c(NA, i)) +
#     facet_wrap(~station)
#   print(g)
# }

#' 
#' ### Detection of rainfall on the same days
#' 
#' ## Rain intensity categories
## ----rain_cats----------------------------------------------------------------------
cat_labs <- c("Dry", "Light Rain", 
              "Moderate Rain", "Heavy Rain", 
              "Violent Rain", "Extreme")
zm_cats <- zm_long_st %>%
  mutate(rain_cats = cut(rain, c(0, 0.85, 5, 20, 40, 100, Inf), include.lowest = TRUE,
                         right = FALSE, labels = cat_labs),
         pr_rain_cats = cut(pr_rain, c(0, 0.85, 5, 20, 40, 100, Inf), include.lowest = TRUE,
                         right = FALSE, labels = cat_labs))

#' 
## ----rain_cats_overall--------------------------------------------------------------
zm_acc_all <- zm_cats %>%
  group_by(station, product) %>%
  summarise(accuracy = mean(rain_cats == pr_rain_cats, na.rm = TRUE))

zm_acc_all %>% 
  mutate(product = toupper(substr(product, 1, nchar(as.character(product)) - 5))) %>%
  pivot_wider(id_cols = "station", names_from = "product", values_from = "accuracy") %>%
  kable(digits = 2) %>%
  skable()

#' 
## ----rain_cats_each-----------------------------------------------------------------
zm_acc_each <- zm_cats %>%
  group_by(station, product, rain_cats) %>%
  filter(!is.na(rain_cats)) %>%
  summarise(n = n(),
            accuracy = mean(rain_cats == pr_rain_cats, na.rm = TRUE))

ggplot(zm_acc_each, aes(x = rain_cats, y = accuracy, fill = product, group = product)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c25[1:5], 
                    labels = c("CHIRPS", "ERA5", "IMERG")) +
  labs(x = "Rainfall Intensity categories (mm/day)", y = "Probability of detection") +
  theme(axis.text.x = element_text(size = 8)) +
  labs(fill = "Product") +
  facet_wrap(~station)

ggsave(here("results", "picsa-paper", 
            paste0("husbands_", "rain_cats.png")),
       width = 12, height = 6)

#' 
#' ![](cont_table.jpg)
#' 
#' - Reasonably good accuracy for all products ~70%, a bit less than in Zambia (75%)
#' - CHIRPS and IMERG miss most of the rain days, only 20%-30% detection, ERA5 get 70%
#' - All also estimate a large number of rain days that did not rain
#' - ERA5 performs best on detecting rain days, even adjusted for higher number of rainy day, but relatively low compared to Zambia
#' - ERA5 misses fewest rain days in station
#' - ERA5 performs best on overall scores, but relatively low compared to Zambia
#' 
## -----------------------------------------------------------------------------------
raindays <- zm_long_st %>% 
  mutate(station0p85 = rain > 0.85,
         pr_rainday0p85 = pr_rain > 0.85,
         pr_rainday3 = pr_rain > 3,
         pr_rainday5 = pr_rain > 5)

sverify <- function(df, obs, pred) {
  v <- verification::verify(obs = df[[obs]], pred = df[[pred]], 
                       obs.type = "binary", frcst.type = "binary")
  class(v) <- "list"
  v
}

by_st_pr <- raindays %>%
  group_by(station, product) %>%
  nest() %>%
  mutate(v = purrr::map(data, sverify, obs = "station0p85", pred = "pr_rainday0p85"),
         bias = purrr::map_dbl(v, "BIAS"))

#' 
## ----cont_tables, results="asis"----------------------------------------------------
cont_tables <- function(df, name) {
  df %>%
    kable(caption = name, digits = 2, format.args = list(big.mark = ",")) %>%
    skable() %>%
    row_spec(nrow(df), bold = TRUE)
}

rain_levs <- c("rain", "no_rain")
ver_zm <- raindays %>%
  group_by(product, station) %>%
  filter(month %in% c(11:12, 1:4)) %>%
  filter(!is.na(station0p85) & !is.na(pr_rainday0p85)) %>%
  summarise(hit = sum(station0p85 & pr_rainday0p85),
            fa = sum(!station0p85 & pr_rainday0p85),
            miss = sum(station0p85 & !pr_rainday0p85),
            cneg = sum(!station0p85 & !pr_rainday0p85),
            n = n(),
            accuracy = (hit + cneg)/n,
            bias = (hit + fa)/(hit + miss),
            hit_rate = hit/(hit + miss),
            far = fa/(hit + fa),
            ts = hit/(hit + miss + fa),
            ets = verification::verify(obs = station0p85, pred = pr_rainday0p85, 
                                       obs.type = "binary", frcst.type = "binary")$ETS,
            miss_frac = miss/n,
            hk = (hit/(hit + miss) - (fa/(fa + cneg))),
            hss = verification::verify(obs = station0p85, pred = pr_rainday0p85, 
                                       obs.type = "binary", frcst.type = "binary")$HSS
            ) %>%
  arrange(product, station)

measures <- c("accuracy", "bias", "hit_rate", "far", "ts", "ets", "miss_frac", "hk", "hss")
names(measures) <- c("Accuracy: What fraction of the estimates were correct?
                     (hits + correct negative)/total (1 = perfect)",
                     "Bias: Ratio of number of rain days from estimate over number of rain days from station. 
                     (hits + false alarms)/(hits + misses) (1 = perfect)",
                     "Hit rate (probability of detection) What fraction of the station rain days were correctly estimated? hits/(hits + misses) (1 = perfect)",
                     "False alarm ratio: What fraction of the estimated rain days actually did not rain? false alarms/(hits + false alarms)",
                     "Threat score: How well did the estimate rain days correspond to the observed rain days?",
                     "Equitable threat score: How well did the estimate rain days correspond to the station rain days (accounting for hits due to chance)?",
                     "Proportion of misses: How many rain days in the station were not detected by the estimate?",
                     "Hanssen and Kuipers discriminant: How well did the forecast separate the rainy days from the dry days? Uses all elements in the contingency table.",
                     "Heidke skill score: What was the accuracy of the forecast relative to that of random chance? Measures the fraction of correct estimates after eliminating those estimates which would be correct due purely to random chance."
                     )
for(i in seq(measures)) {
  df <- ver_zm %>% 
    select(product, station, measures[[i]]) %>%
    pivot_wider(id_cols = station, names_from = product, values_from = measures[i])
  names(df)[endsWith(names(df), "rain")] <- substr(names(df)[endsWith(names(df), "rain")],
                                                   1, 
                                                   nchar(names(df)
                                                         [endsWith(names(df), "rain")]
                                                         ) 
                                                   - 5)
  df[nrow(df) + 1, ] <- c(list("(All)"), as.list(as.numeric(colMeans(df[ , -1]))))
  print(cont_tables(df, names(measures)[i]))
}

#' 
