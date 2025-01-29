library(here)
library(ncdf4)
library(dplyr)
library(lubridate)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggrepel)
library(tidync)
library(gridExtra) # for combining the graphs
library(sf)

maps_data <- "~/Documents/John/PHD 2021/IDEMS work/epicsa-analysis/results/maps"
# Read the data
total_means_agera5_zambia <- readRDS(here(maps_data, "agera5_zambia_total.RDS"))
n_means_agera5_zambia <- readRDS(here(maps_data, "agera5_zambia_n.RDS"))
mean_means_agera5_zambia <- readRDS(here(maps_data, "agera5_zambia_mean.RDS"))

total_means_tamsat_zambia <- readRDS(here(maps_data, "tamsat_zambia_total.RDS"))
n_means_tamsat_zambia <- readRDS(here(maps_data, "tamsat_zambia_n.RDS"))
mean_means_tamsat_zambia <- readRDS(here(maps_data, "tamsat_zambia_mean.RDS"))

total_means_era5_zambia <- readRDS(here(maps_data, "era5_zambia_total.RDS"))
n_means_era5_zambia <- readRDS(here(maps_data, "era5_zambia_n.RDS"))
mean_means_era5_zambia <- readRDS(here(maps_data, "era5_zambia_mean.RDS"))

enacts2019_zambia <- readRDS(here("results", "maps", "enacts_zambia_2019.RDS"))
enacts_old_t_rain <- enacts2019_zambia %>% 
  select(lon, lat, enacts_total) %>%
  filter(!is.na(enacts_total)) %>%
  rename(total_rain = enacts_total)

enacts_old_n_rain <- enacts2019_zambia %>% 
  select(lon, lat, enacts_n) %>%
  filter(!is.na(enacts_n)) %>%
  rename(n_rain = enacts_n)

enacts_old_mean_rain <- enacts2019_zambia %>% 
  select(lon, lat, enacts_mean) %>%
  filter(!is.na(enacts_mean)) %>%
  rename(mean_rain = enacts_mean)

total_means_chirps_zambia <- readRDS(here(maps_data, "chirps_zambia_total.RDS"))
n_means_chirps_zambia <- readRDS(here(maps_data, "chirps_zambia_n.RDS"))
mean_means_chirps_zambia <- readRDS(here(maps_data, "chirps_zambia_mean.RDS"))

zambia_metadata <- readRDS(here("data", "station", "processed", "zambia_stations_metadata_updated.RDS"))

zambia_sf <- ne_countries(country = "Zambia", returnclass = "sf")

filter_data_zambia <- function(data, zambia_sf) {
  data_sf <- st_as_sf(data, coords = c("lon", "lat"), crs = st_crs(zambia_sf), remove = FALSE)
  data_zambia <- st_intersection(data_sf, zambia_sf)
  data_zambia <- as.data.frame(data_zambia)
  return(data_zambia)
}

plot_T_rainfall_map <- function(data, title) {
  ggplot() + 
    geom_raster(data = filter_data_zambia(data, zambia_sf), aes(x = lon, y = lat, fill = total_rain)) +
    geom_sf(data = zambia_sf, fill = NA, colour = "black") +
    geom_text_repel(data = zambia_metadata, aes(x = longitude, y = latitude, label = "")) +
    labs(title = title, fill = "Rainfall (mm)", x = "", y = "") +
    scale_fill_gradientn(
      colours = rev(viridis::viridis(256)),  
      values = scales::rescale(c(500, 2000)),  
      limits = c(500, 2000), 
      oob = scales::squish  
    ) +
    coord_sf() +
    theme(
      panel.background = element_rect(fill = "white", color = "white"),  
      plot.background = element_rect(fill = "white", color = "white"),  
      axis.text = element_blank(), 
      axis.ticks = element_blank(), 
      text = element_text(size = 10)
    )
}

plot_n_rainfall_map <- function(data, title) {
  ggplot() + 
    geom_raster(data = filter_data_zambia(data, zambia_sf), aes(x = lon, y = lat, fill = n_rain)) +
    geom_sf(data = zambia_sf, fill = NA, colour = "black") +
    geom_text_repel(data = zambia_metadata, aes(x = longitude, y = latitude, label = "")) +
    labs(title = title, fill = "Number of Rain Days\n(rainday/year)", x = "", y = "") +
    scale_fill_gradientn(
      colours = rev(viridis::viridis(256)),  
      values = scales::rescale(c(50, 150)), 
      limits = c(50, 150), 
      oob = scales::squish  
    ) +
    coord_sf() + 
    theme(
      panel.background = element_rect(fill = "white", color = "white"), 
      plot.background = element_rect(fill = "white", color = "white"),  
      axis.text = element_blank(),  
      axis.ticks = element_blank(),  
      text = element_text(size = 10)
    )
}

plot_m_rain_rainfall_map <- function(data, title) {
  ggplot() + 
    geom_raster(data = filter_data_zambia(data, zambia_sf), aes(x = lon, y = lat, fill = mean_rain)) +
    geom_sf(data = zambia_sf, fill = NA, colour = "black") +
    geom_text_repel(data = zambia_metadata, aes(x = longitude, y = latitude, label = "")) +
    labs(title = title, fill = "Mean Rain per Rainday", x = "", y = "") +
    scale_fill_gradientn(
      colours = rev(viridis::viridis(256)),  
      values = scales::rescale(c(5, 20)),  
      limits = c(5, 20), 
      oob = scales::squish 
    ) +
    coord_sf() +
    theme(
      panel.background = element_rect(fill = "white", color = "white"),  
      plot.background = element_rect(fill = "white", color = "white"), 
      axis.text = element_blank(),  
      axis.ticks = element_blank(),  
      text = element_text(size = 10)
    )
}

#plot new enacts
files <- list.files("/media/johnbagiliko/Seagate Backup Plus Drive/ENACTS/new_enacts", 
                    pattern = "data_.*\\.nc", full.names = TRUE)
# Extract longitude, latitude
nc <- nc_open(files[1])
lon <- ncvar_get(nc, "X")
lat <- ncvar_get(nc, "Y")
nc_close(nc)

mean_rainfall <- readRDS(here(maps_data, "new_enacts", "array_data", "mean_rainfall_array.rds"))
n_rainfall <- readRDS(here(maps_data, "new_enacts", "array_data", "mean_rainy_days_array.rds"))
mean_r_per_rday <- readRDS(here(maps_data, "new_enacts", "array_data", "mean_rain_per_day_array.rds"))

create_spatial_df <- function(data, lon, lat, name) {
  df <- as.data.frame(as.table(data))
  colnames(df) <- c("lon_idx", "lat_idx", name)
  df$lon <- lon[df$lon_idx]
  df$lat <- lat[df$lat_idx]
  st_as_sf(df, coords = c("lon", "lat"), crs = st_crs(zambia_sf), remove = FALSE)
}

enacts_t_rain_sf <- create_spatial_df(data = mean_rainfall, lon, lat, "total_rain")
enacts_n_rain_sf <- create_spatial_df(data = n_rainfall, lon, lat, "n_rain")
enacts_mean_rain_sf <- create_spatial_df(data = mean_r_per_rday, lon, lat, "mean_rain")

reshape_enacts_for_map <- function(sf_data, val_col_name){
  rainfall_sf_enacts <- sf_data %>% 
    ungroup() %>%
    select(lon, lat, {{ val_col_name }})
  rainfall_sf_enacts$geometry <- NULL
  return(rainfall_sf_enacts)
}

#total rainfall
t_rain_df_enacts <- reshape_enacts_for_map(enacts_t_rain_sf, total_rain)
enacts_old_total_rain <- plot_T_rainfall_map(enacts_old_t_rain, "(a) ENACTS_old")
era5_total_rain <- plot_T_rainfall_map(total_means_era5_zambia, "(b) ERA5")
chirps_total_rain <- plot_T_rainfall_map(total_means_chirps_zambia, "(c) CHIRPS")
tamsat_total_rain <- plot_T_rainfall_map(total_means_tamsat_zambia, "(d) TAMSAT")
agera5_total_rain <- plot_T_rainfall_map(total_means_agera5_zambia, "(e) AGERA5")
enacts_total_rain <- plot_T_rainfall_map(t_rain_df_enacts, "(a) ENACTS")



#number of rainy days
n_rain_df_enacts <- reshape_enacts_for_map(enacts_n_rain_sf, n_rain)
enacts_old_n_rain = plot_n_rainfall_map(enacts_old_n_rain, "(a) ENACTS_old")
era5_n_rain <- plot_n_rainfall_map(n_means_era5_zambia, "(b) ERA5")
chirps_n_rain <- plot_n_rainfall_map(n_means_chirps_zambia, "(c) CHIRPS")
tamsat_n_rain <- plot_n_rainfall_map(n_means_tamsat_zambia, "(d) TAMSAT")
agera5_n_rain <- plot_n_rainfall_map(n_means_agera5_zambia, "(e) AGERA5")
enacts_n_rain <- plot_n_rainfall_map(n_rain_df_enacts, "(f) ENACTS")


#mean rain per rainy day
mean_rain_df_enacts <- reshape_enacts_for_map(enacts_mean_rain_sf, mean_rain)
enacts_old_mean_rain = plot_m_rain_rainfall_map(enacts_old_mean_rain, "(a) ENACTS_old")
era5_mean_rain <- plot_m_rain_rainfall_map(mean_means_era5_zambia, "(b) ERA5")
chirps_mean_rain <- plot_m_rain_rainfall_map(mean_means_chirps_zambia, "(c) CHIRPS")
tamsat_mean_rain <- plot_m_rain_rainfall_map(mean_means_tamsat_zambia, "(d) TAMSAT")
agera5_mean_rain <- plot_m_rain_rainfall_map(mean_means_agera5_zambia, "(e) AGERA5")
enacts_mean_rain <- plot_m_rain_rainfall_map(mean_rain_df_enacts, "(f) ENACTS")

combined_plot_t_rain <- grid.arrange(enacts_old_total_rain, era5_total_rain, chirps_total_rain, tamsat_total_rain, agera5_total_rain, enacts_total_rain, ncol = 2)

combined_plot_n_rain <- grid.arrange(enacts_old_n_rain, era5_n_rain, chirps_n_rain, tamsat_n_rain, agera5_n_rain, enacts_n_rain, ncol = 2)

combined_plot_mean_rain <- grid.arrange(enacts_old_mean_rain, era5_mean_rain, chirps_mean_rain, tamsat_mean_rain, agera5_mean_rain, enacts_mean_rain, ncol = 2)


ggsave(here("results", "maps", "mean_total_rain_combined_plot_new.png"), plot = combined_plot_t_rain, width = 10, height = 8, dpi = 300)

ggsave(here("results", "maps", "mean_n_rain_combined_plot_new.png"), plot = combined_plot_n_rain, width = 10, height = 8, dpi = 300)

ggsave(here("results", "maps", "mean_mean_rain_combined_plot_new.png"), plot = combined_plot_mean_rain, width = 10, height = 8, dpi = 300)

