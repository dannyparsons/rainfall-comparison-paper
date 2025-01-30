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

zm <- readRDS(here("data", "station", "cleaned", "westafrica_gridded.RDS"))
# 1 Jan = 1
s_doy_start <- 1
zm <- zm %>% 
  mutate(doy = yday_366(date),
         s_doy = (doy - s_doy_start + 1) %% 366,
         s_doy = ifelse(s_doy == 0, 366, s_doy),
         syear = year,
         syear = ifelse(s_doy > (366 - s_doy_start + 1), syear - 1, syear),
         month = factor(month),
         month_abb = month(date, label = TRUE)
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
  filter(country %in% c("Ghana", "Niger"))
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

## ----prop_complete_days-----------------------------------------------------------------
zm %>% 
  group_by(station) %>% 
  summarise(prop_n = 100* (1 - naflex::na_prop(rain))) %>%
  kable(digits = 1) %>%
  skable()

## ----markov_chain_setup1----------------------------------------------------------------
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
}
predict_stack_all <- bind_rows(predict_stack_lst)

## ----markov_chain_plots_products--------------------------------------------------------
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
    ##! Figure: wa_markov_zero_chirps (graphs for other products not used in paper)
  ggsave(here("results", 
              paste0("westafrica_", "markov_zero_", names(products)[s], ".png")),
         plot = g, width = 12, height = 6)
}
