library(here)
library(dplyr)
library(ggplot2)
library(lubridate)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggrepel)
library(tidync)
library(ncdf4)
library(sp)
library(ggspatial)
library(rnaturalearthhires)

source(here("src", "helper_funs.R"))

husbands <- readRDS(here("data", "station", "cleaned", "husbands_qc.RDS"))
husbands_metadata <- readRDS(here("data", "station", "processed", "barbados_station_metadata.RDS"))

# station location check --------------------------------------------------
husbands_sf <- ne_countries(scale = "large", returnclass = "sf")
ggplot(husbands_sf) +
  geom_sf() +
  geom_point(data = husbands_metadata, aes(x = longitude, y = latitude)) +
  geom_text_repel(data = husbands_metadata, aes(x = longitude, y = latitude, label = station)) +
  coord_sf(xlim = c(-59.7, -59.4), ylim = c(13, 13.4), expand = FALSE)

# CHIRPS Import -----------------------------------------------------------

files <- list.files("D:/data/chirps", pattern = "chirps_husbands", full.names = TRUE)

nc <- nc_open(files[1])
stopifnot(min(nc$dim$X$vals) < min(husbands_metadata$longitude))
stopifnot(max(nc$dim$X$vals) > max(husbands_metadata$longitude))
stopifnot(min(nc$dim$Y$vals) < min(husbands_metadata$latitude))
stopifnot(max(nc$dim$Y$vals) > max(husbands_metadata$latitude))
xs <- nc$dim$X$vals
ys <- nc$dim$Y$vals
resx <- xs[2] - xs[1]
resy <- ys[2] - ys[1]
max_dist <- sqrt((resx/2)^2 + (resy/2)^2)
xy_points <- expand.grid(xs, ys)
xy_extract <- closest_point(points = xy_points, 
                            target = husbands_metadata %>% dplyr::select(longitude, latitude))
xy_extract$station <- husbands_metadata$station
xy_extract$req_longitude <- husbands_metadata$longitude
xy_extract$req_latitude <- husbands_metadata$latitude
xy_extract$dist <- apply(xy_extract, 1, function(r) {
  sp::spDistsN1(matrix(as.numeric(c(r[["Var1"]], r[["Var2"]])), ncol = 2),
                as.numeric(c(r[["req_longitude"]], r[["req_latitude"]])))
}
)
stopifnot(all(xy_extract$dist <= max_dist))
print(nc$dim$`T`$units)
nc_close(nc)

pr_loc <- data.frame(pr = "CHIRPS", latitude = xy_extract$Var2, 
                     longitude = xy_extract$Var1, width = 0.05, stringsAsFactors = FALSE)

chirps_dfs <- list()
pb <- txtProgressBar(min = 0, max = length(files) * nrow(xy_extract), style = 3)
count <- 0
for(f in files) {
  for(i in seq_len(nrow(xy_extract))) {
    df <- tidync(f) %>%
      hyper_filter(X = X == xy_extract[i, 1], 
                   Y = Y == xy_extract[i, 2]) %>%
      hyper_tibble()
    df$station <- xy_extract$station[i]
    df$req_latitude <- xy_extract$req_latitude[i]
    df$req_longitude <- xy_extract$req_longitude[i]
    df$dist <- xy_extract$dist[i]
    chirps_dfs[[length(chirps_dfs) + 1]] <- df
    count <- count + 1
    setTxtProgressBar(pb, count)
  }
}
husbands_chirps <- bind_rows(chirps_dfs) %>%
  mutate(date = as.Date(`T`, origin = structure(-2440588, class = "Date")))

husbands_chirps <- husbands_chirps %>% select(station, date, chirps_rain = prcp)

saveRDS(husbands_chirps, here("data", "station", "cleaned", "husbands_chirps.RDS"))

husbands <- left_join(husbands, husbands_chirps, by = c("station", "date"))

# ERA5 Import -------------------------------------------------------------

files <- list.files("D:/data/era5", pattern = "caribbean", full.names = TRUE)

nc <- nc_open(files[1])
stopifnot(min(nc$dim$longitude$vals) < min(husbands_metadata$longitude))
stopifnot(max(nc$dim$longitude$vals) > max(husbands_metadata$longitude))
stopifnot(min(nc$dim$latitude$vals) < min(husbands_metadata$latitude))
stopifnot(max(nc$dim$latitude$vals) > max(husbands_metadata$latitude))
xs <- nc$dim$longitude$vals
ys <- nc$dim$latitude$vals
resx <- xs[2] - xs[1]
resy <- ys[2] - ys[1]
max_dist <- sqrt((resx/2)^2 + (resy/2)^2)
xy_points <- expand.grid(xs, ys)
xy_extract <- closest_point(points = xy_points, 
                            target = husbands_metadata %>% dplyr::select(longitude, latitude))
xy_extract$station <- husbands_metadata$station
xy_extract$req_longitude <- husbands_metadata$longitude
xy_extract$req_latitude <- husbands_metadata$latitude
xy_extract$dist <- apply(xy_extract, 1, function(r) {
  sp::spDistsN1(matrix(as.numeric(c(r[["Var1"]], r[["Var2"]])), ncol = 2),
                as.numeric(c(r[["req_longitude"]], r[["req_latitude"]])))
}
)
stopifnot(all(xy_extract$dist <= max_dist))
print(nc$dim$time$units)
nc_close(nc)

pr_loc[2, ] <- c(pr = "ERA5", latitude = round(xy_extract$Var2, 3), 
                     longitude = round(xy_extract$Var1, 3), width = 0.25)

era5_dfs <- list()
pb <- txtProgressBar(min = 0, max = length(files) * nrow(xy_extract), style = 3)
count <- 0
for(f in files) {
  for(i in seq_len(nrow(xy_extract))) {
    df <- tidync(f) %>%
      hyper_filter(longitude = longitude == xy_extract[i, 1], 
                   latitude = latitude == xy_extract[i, 2]) %>%
      hyper_tibble()
    df$station <- xy_extract$station[i]
    df$req_longitude <- xy_extract$req_longitude[i]
    df$req_latitude <- xy_extract$req_latitude[i]
    df$dist <- xy_extract$dist[i]
    era5_dfs[[length(era5_dfs) + 1]] <- df
    count <- count + 1
    setTxtProgressBar(pb, count)
  }
}
husbands_era5 <- bind_rows(era5_dfs) %>%
  mutate(date_time = as.POSIXct(time * 60 * 60,
                                origin = as.POSIXct("1900/01/01 00:00:00"), tz = "GMT"),
         date = as.Date(date_time),
         # Assuming 8am rain recording time in Barbados, (Barbados time is UTC - 4)
         # ERA5 record at 9am local time is for period 8am-9am (first recording hour)
         # Formula: lead = start_time - UTC_diff
         # e.g. Barbados start_time = 9am, UTC_diff = -4, lead = 9 - (-4) = 13
         # e.g. Zambia start_time = 9am, UTC_diff = 2, lead = 9 - (2) = 7
         tp_lead_13 = dplyr::lead(tp, 13)) %>%
  group_by(station, date) %>%
  summarise(tp = sum(tp_lead_13) * 1000)

husbands_era5 <- husbands_era5 %>% select(station, date, era5_rain = tp)

saveRDS(husbands_era5, here("data", "station", "cleaned", "husbands_era5.RDS"))

husbands <- left_join(husbands, husbands_era5, by = c("station", "date"))

# IMERG Import ------------------------------------------------------------

files <- list.files("D:/data/imerg", pattern = "caribbean", full.names = TRUE)

nc <- nc_open(files[1])
stopifnot(min(nc$dim$lon$vals) < min(husbands_metadata$longitude))
stopifnot(max(nc$dim$lon$vals) > max(husbands_metadata$longitude))
stopifnot(min(nc$dim$lat$vals) < min(husbands_metadata$latitude))
stopifnot(max(nc$dim$lat$vals) > max(husbands_metadata$latitude))
xs <- nc$dim$lon$vals
ys <- nc$dim$lat$vals
resx <- xs[2] - xs[1]
resy <- ys[2] - ys[1]
max_dist <- sqrt((resx/2)^2 + (resy/2)^2)
xy_points <- expand.grid(xs, ys)
xy_extract <- closest_point(points = xy_points, 
                            target = husbands_metadata %>% dplyr::select(longitude, latitude))
xy_extract$station <- husbands_metadata$station
xy_extract$req_longitude <- husbands_metadata$longitude
xy_extract$req_latitude <- husbands_metadata$latitude
xy_extract$dist <- apply(xy_extract, 1, function(r) {
  sp::spDistsN1(matrix(as.numeric(c(r[["Var1"]], r[["Var2"]])), ncol = 2),
                as.numeric(c(r[["req_longitude"]], r[["req_latitude"]])))
}
)
stopifnot(all(xy_extract$dist <= max_dist))
print(nc$dim$time$units)
nc_close(nc)

pr_loc[3, ] <- c(pr = "IMERG", latitude = round(xy_extract$Var2, 3), 
                 longitude = round(xy_extract$Var1, 3), width = 0.1)

imerg_dfs <- list()
pb <- txtProgressBar(min = 0, max = length(files) * nrow(xy_extract), style = 3)
count <- 0
for(f in files) {
  print(f)
  for(i in seq_len(nrow(xy_extract))) {
    df <- tidync(f) %>%
      hyper_filter(lon = lon == xy_extract[i, 1], 
                   lat = lat == xy_extract[i, 2]) %>%
      hyper_tibble()
    df$station <- xy_extract$station[i]
    df$req_longitude <- xy_extract$req_longitude[i]
    df$req_latitude <- xy_extract$req_latitude[i]
    df$dist <- xy_extract$dist[i]
    imerg_dfs[[length(imerg_dfs) + 1]] <- df
    count <- count + 1
    setTxtProgressBar(pb, count)
  }
}
husbands_imerg <- bind_rows(imerg_dfs) %>%
  mutate(date = as.Date(time, origin = as.Date("1970/01/01")))

husbands_imerg <- husbands_imerg %>% select(station, date, imerg_rain = precipitationCal)

saveRDS(husbands_imerg, here("data", "station", "cleaned", "husbands_imerg.RDS"))

husbands <- left_join(husbands, husbands_imerg, by = c("station", "date"))

saveRDS(husbands, here("data", "station", "cleaned", "husbands_gridded.RDS"))


# Plot areas --------------------------------------------------------------

pr_loc$latitude <- as.numeric(pr_loc$latitude)
pr_loc$longitude <- as.numeric(pr_loc$longitude)
pr_loc$width <- as.numeric(pr_loc$width)
pr_loc <- pr_loc %>% arrange(-width)

ggplot(husbands_sf) +
  geom_sf(fill= "antiquewhite") +
  geom_point(data = husbands_metadata, aes(x = longitude, y = latitude)) +
  geom_rect(data = pr_loc, aes(xmin = longitude - width/2, xmax = longitude + width/2,
                               ymin = latitude - width/2, ymax = latitude + width/2, 
                               fill = pr), alpha = 0.2) +
  annotate(geom = "text", x = -59.62, y = 13.17, 
           label = "Husbands", size = 4) +
  geom_text(data = pr_loc, aes(x = longitude, y = latitude, label = pr),
            size = 3) +
  annotate(geom = "text", x = -59.55, y = 13.2, label = "Barbados", 
           fontface = "italic", color = "grey22", size = 6) +
  coord_sf(xlim = c(-59.9, -59.4), ylim = c(13, 13.4), expand = FALSE) +
  scale_x_continuous(breaks = seq(-59.8, -59.4, 0.1)) +
  scale_y_continuous(breaks = seq(13, 13.3, 0.1)) +
  annotation_scale(location = "bl", width_hint = 0.5) + 
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.5, "in"), pad_y = unit(0.5, "in"), 
                         style = north_arrow_fancy_orienteering) +
  theme_bw() +
  theme(panel.grid.major = element_line(color = gray(0.7), linetype = "dashed", size = 0.5), 
        panel.background = element_rect(fill = "aliceblue"))
ggsave(here("results", "barbados", "barbados_map.png"), width = 12, height = 6)
