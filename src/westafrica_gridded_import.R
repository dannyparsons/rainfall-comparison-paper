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

source(here("src", "helper_funs.R"))

westafrica <- readRDS(here("data", "station", "cleaned", "westafrica_1979_qc.RDS"))
station_metadata <- readRDS(here("data", "station", "processed", "stations_metadata.RDS"))
westafrica_metadata <- station_metadata %>%
  filter(country %in% c("Niger", "Ghana"))


# station location check --------------------------------------------------
sf_wa <- ne_countries(returnclass = "sf")
ggplot(sf_wa) + 
  geom_sf() +
  geom_point(data = westafrica_metadata, aes(x = longitude, y = latitude)) +
  geom_text_repel(data = westafrica_metadata, aes(x = longitude, y = latitude, label = station)) +
  coord_sf(xlim = c(-8, 15), ylim = c(0, 25), expand = FALSE)

rm(sf_wa)

# ARC2 Import -------------------------------------------------------------

files <- list.files("D:/data/arc2", pattern = "arc2-westafrica", full.names = TRUE)

nc <- nc_open(files[1])
stopifnot(min(nc$dim$X$vals) < min(westafrica_metadata$longitude))
stopifnot(max(nc$dim$X$vals) > max(westafrica_metadata$longitude))
stopifnot(min(nc$dim$Y$vals) < min(westafrica_metadata$latitude))
stopifnot(max(nc$dim$Y$vals) > max(westafrica_metadata$latitude))
xs <- nc$dim$X$vals
ys <- nc$dim$Y$vals
resx <- xs[2] - xs[1]
resy <- ys[2] - ys[1]
max_dist <- sqrt((resx/2)^2 + (resy/2)^2)
xy_points <- expand.grid(xs, ys)
xy_extract <- closest_point(points = xy_points, 
                  target = westafrica_metadata %>% select(longitude, latitude))
xy_extract$station <- westafrica_metadata$station
xy_extract$req_longitude <- westafrica_metadata$longitude
xy_extract$req_latitude <- westafrica_metadata$latitude
xy_extract$dist <- apply(xy_extract, 1, function(r) {
  sp::spDistsN1(matrix(as.numeric(c(r[["Var1"]], r[["Var2"]])), ncol = 2),
                as.numeric(c(r[["req_longitude"]], r[["req_latitude"]])))
  }
)
stopifnot(all(xy_extract$dist <= max_dist))
print(nc$dim$`T`$units)
nc_close(nc)

arc_dfs <- list()
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
    arc_dfs[[length(arc_dfs) + 1]] <- df
  }
}
westafrica_arc2 <- bind_rows(arc_dfs) %>%
  mutate(date = as.Date(`T`, origin = as.Date("1960/01/01")))

westafrica_arc2 <- westafrica_arc2 %>% select(station, date, arc2_rain = est_prcp)

saveRDS(westafrica_arc2, here("data", "station", "cleaned", "westafrica_arc2.RDS"))

westafrica <- left_join(westafrica, eswatini_arc2, by = c("station", "date"))

# CHIRPS Import -----------------------------------------------------------

files <- list.files("D:/data/chirps", pattern = "chirps_westafrica", full.names = TRUE)

nc <- nc_open(files[1])
stopifnot(min(nc$dim$X$vals) < min(westafrica_metadata$longitude))
stopifnot(max(nc$dim$X$vals) > max(westafrica_metadata$longitude))
stopifnot(min(nc$dim$Y$vals) < min(westafrica_metadata$latitude))
stopifnot(max(nc$dim$Y$vals) > max(westafrica_metadata$latitude))
xs <- nc$dim$X$vals
ys <- nc$dim$Y$vals
resx <- xs[2] - xs[1]
resy <- ys[2] - ys[1]
max_dist <- sqrt((resx/2)^2 + (resy/2)^2)
xy_points <- expand.grid(xs, ys)
xy_extract <- closest_point(points = xy_points, 
                            target = westafrica_metadata %>% select(longitude, latitude))
xy_extract$station <- westafrica_metadata$station
xy_extract$req_longitude <- westafrica_metadata$longitude
xy_extract$req_latitude <- westafrica_metadata$latitude
xy_extract$dist <- apply(xy_extract, 1, function(r) {
  sp::spDistsN1(matrix(as.numeric(c(r[["Var1"]], r[["Var2"]])), ncol = 2),
                as.numeric(c(r[["req_longitude"]], r[["req_latitude"]])))
}
)
stopifnot(all(xy_extract$dist <= max_dist))
print(nc$dim$`T`$units)
nc_close(nc)

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
westafrica_chirps <- bind_rows(chirps_dfs) %>%
  mutate(date = as.Date(`T`, origin = structure(-2440588, class = "Date")))

westafrica_chirps <- westafrica_chirps %>% select(station, date, chirps_rain = prcp)

saveRDS(westafrica_chirps, here("data", "station", "cleaned", "westafrica_chirps.RDS"))

westafrica <- left_join(westafrica, westafrica_chirps, by = c("station", "date"))

# ERA5 Import -------------------------------------------------------------

files <- list.files("D:/data/era5", pattern = "africa", full.names = TRUE)

nc <- nc_open(files[1])
stopifnot(min(nc$dim$longitude$vals) < min(westafrica_metadata$longitude))
stopifnot(max(nc$dim$longitude$vals) > max(westafrica_metadata$longitude))
stopifnot(min(nc$dim$latitude$vals) < min(westafrica_metadata$latitude))
stopifnot(max(nc$dim$latitude$vals) > max(westafrica_metadata$latitude))
xs <- nc$dim$longitude$vals
ys <- nc$dim$latitude$vals
resx <- xs[2] - xs[1]
resy <- ys[2] - ys[1]
max_dist <- sqrt((resx/2)^2 + (resy/2)^2)
xy_points <- expand.grid(xs, ys)
xy_extract <- closest_point(points = xy_points, 
                            target = westafrica_metadata %>% select(longitude, latitude))
xy_extract$station <- westafrica_metadata$station
xy_extract$req_longitude <- westafrica_metadata$longitude
xy_extract$req_latitude <- westafrica_metadata$latitude
xy_extract$dist <- apply(xy_extract, 1, function(r) {
  sp::spDistsN1(matrix(as.numeric(c(r[["Var1"]], r[["Var2"]])), ncol = 2),
                as.numeric(c(r[["req_longitude"]], r[["req_latitude"]])))
}
)
stopifnot(all(xy_extract$dist <= max_dist))
print(nc$dim$time$units)
nc_close(nc)

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
westafrica_era5 <- bind_rows(era5_dfs) %>%
  mutate(date_time = as.POSIXct(time * 60 * 60,
                                origin = as.POSIXct("1900/01/01 00:00:00"), tz = "GMT"),
         date = as.Date(date_time),
         # Assuming 8am rain recording time in Ghana/Niger
         # Ghana time is UTC + 0, Niger time is UTC + 1
         # ERA5 record at 9am local time is for period 8am-9am (first recording hour)
         # Formula: lead = start_time - UTC_diff
         # e.g. Ghana start_time = 9am, UTC_diff = 0, lead = 9 - 0 = 9
         # e.g. Niger start_time = 9am, UTC_diff = 1, lead = 9 - 1 = 8
         tp_lead = ifelse(station == "Sadore", dplyr::lead(tp, 8), dplyr::lead(tp, 9)),
         ) %>%
  group_by(station, date) %>%
  summarise(tp = sum(tp_lead) * 1000, 
            np0p01 = sum((tp_lead * 1000) > 0.01),
            np0p025 = sum((tp_lead * 1000) > 0.025),
            max_tp = max((tp_lead * 1000)),
            sp0p025 = sum(rle((tp_lead * 1000) > 0.025)$values),
            sp0p01 = sum(rle((tp_lead * 1000) > 0.01)$values))

westafrica_era5_all <- westafrica_era5
westafrica_era5_all <- left_join(westafrica, westafrica_era5_all, by = c("station", "date"))

saveRDS(westafrica_era5_all, here("data", "station", "cleaned", "westafrica_era5_all.RDS"))

westafrica_era5 <- westafrica_era5 %>% select(station, date, era5_rain = tp)

saveRDS(westafrica_era5, here("data", "station", "cleaned", "westafrica_era5.RDS"))

westafrica <- left_join(westafrica, westafrica_era5, by = c("station", "date"))

# IMERG Import ------------------------------------------------------------

files <- list.files("D:/data/imerg", pattern = "imerg", full.names = TRUE)

nc <- nc_open(files[1])
stopifnot(min(nc$dim$lon$vals) < min(westafrica_metadata$longitude))
stopifnot(max(nc$dim$lon$vals) > max(westafrica_metadata$longitude))
stopifnot(min(nc$dim$lat$vals) < min(westafrica_metadata$latitude))
stopifnot(max(nc$dim$lat$vals) > max(westafrica_metadata$latitude))
xs <- nc$dim$lon$vals
ys <- nc$dim$lat$vals
resx <- xs[2] - xs[1]
resy <- ys[2] - ys[1]
max_dist <- sqrt((resx/2)^2 + (resy/2)^2)
xy_points <- expand.grid(xs, ys)
xy_extract <- closest_point(points = xy_points, 
                            target = westafrica_metadata %>% select(longitude, latitude))
xy_extract$station <- westafrica_metadata$station
xy_extract$req_longitude <- westafrica_metadata$longitude
xy_extract$req_latitude <- westafrica_metadata$latitude
xy_extract$dist <- apply(xy_extract, 1, function(r) {
  sp::spDistsN1(matrix(as.numeric(c(r[["Var1"]], r[["Var2"]])), ncol = 2),
                as.numeric(c(r[["req_longitude"]], r[["req_latitude"]])))
}
)
stopifnot(all(xy_extract$dist <= max_dist))
print(nc$dim$time$units)
nc_close(nc)

imerg_dfs <- list()
pb <- txtProgressBar(min = 0, max = length(files) * nrow(xy_extract), style = 3)
count <- 0
for(f in files) {
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
westafrica_imerg <- bind_rows(imerg_dfs) %>%
  mutate(date = as.Date(time, origin = as.Date("1970/01/01")))

westafrica_imerg <- westafrica_imerg %>% select(station, date, imerg_rain = precipitationCal)

saveRDS(westafrica_imerg, here("data", "station", "cleaned", "westafrica_imerg.RDS"))

westafrica <- left_join(westafrica, westafrica_imerg, by = c("station", "date"))

# RFE2 Import -------------------------------------------------------------

files <- list.files("D:/data/rfe2", pattern = "rfe2-westafrica", full.names = TRUE)

nc <- nc_open(files[1])
stopifnot(min(nc$dim$X$vals) < min(westafrica_metadata$longitude))
stopifnot(max(nc$dim$X$vals) > max(westafrica_metadata$longitude))
stopifnot(min(nc$dim$Y$vals) < min(westafrica_metadata$latitude))
stopifnot(max(nc$dim$Y$vals) > max(westafrica_metadata$latitude))
xs <- nc$dim$X$vals
ys <- nc$dim$Y$vals
resx <- xs[2] - xs[1]
resy <- ys[2] - ys[1]
max_dist <- sqrt((resx/2)^2 + (resy/2)^2)
xy_points <- expand.grid(xs, ys)
xy_extract <- closest_point(points = xy_points, 
                            target = westafrica_metadata %>% select(longitude, latitude))
xy_extract$station <- westafrica_metadata$station
xy_extract$req_longitude <- westafrica_metadata$longitude
xy_extract$req_latitude <- westafrica_metadata$latitude
xy_extract$dist <- apply(xy_extract, 1, function(r) {
  sp::spDistsN1(matrix(as.numeric(c(r[["Var1"]], r[["Var2"]])), ncol = 2),
                as.numeric(c(r[["req_longitude"]], r[["req_latitude"]])))
}
)
stopifnot(all(xy_extract$dist <= max_dist))
print(nc$dim$`T`$units)
nc_close(nc)

rfe2_dfs <- list()
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
    rfe2_dfs[[length(rfe2_dfs) + 1]] <- df
    count <- count + 1
    setTxtProgressBar(pb, count)
  }
}
westafrica_rfe2 <- bind_rows(rfe2_dfs) %>%
  mutate(date = as.Date(`T`, origin = as.Date("2000/10/31")))

westafrica_rfe2 <- westafrica_rfe2 %>% select(station, date, rfe2_rain = est_prcp)

eswatini <- left_join(eswatini, eswatini_rfe2, by = c("station", "date"))

saveRDS(westafrica_rfe2, here("data", "station", "cleaned", "westafrica_rfe2.RDS"))

saveRDS(eswatini, here("data", "station", "cleaned", "eswatini_gridded.RDS"))

# TAMSAT3 Import ----------------------------------------------------------

f_tamsat <- "D:/data/tamsat/tamsat_all.nc"

# need to move Saltpond station further inland
westafrica_metadata_tamsat <- westafrica_metadata
westafrica_metadata_tamsat$latitude[1] <- 5.25

nc <- nc_open(f_tamsat)
stopifnot(min(nc$dim$lon$vals) < min(westafrica_metadata_tamsat$longitude))
stopifnot(max(nc$dim$lon$vals) > max(westafrica_metadata_tamsat$longitude))
stopifnot(min(nc$dim$lat$vals) < min(westafrica_metadata_tamsat$latitude))
stopifnot(max(nc$dim$lat$vals) > max(westafrica_metadata_tamsat$latitude))
xs <- nc$dim$lon$vals
ys <- nc$dim$lat$vals
resx <- xs[2] - xs[1]
resy <- ys[2] - ys[1]
max_dist <- sqrt((resx/2)^2 + (resy/2)^2)
xy_points <- expand.grid(xs, ys)
xy_extract <- closest_point(points = xy_points, 
                            target = westafrica_metadata_tamsat %>% select(longitude, latitude))
xy_extract$station <- westafrica_metadata_tamsat$station
xy_extract$req_longitude <- westafrica_metadata_tamsat$longitude
xy_extract$req_latitude <- westafrica_metadata_tamsat$latitude
xy_extract$dist <- apply(xy_extract, 1, function(r) {
  sp::spDistsN1(matrix(as.numeric(c(r[["Var1"]], r[["Var2"]])), ncol = 2),
                as.numeric(c(r[["req_longitude"]], r[["req_latitude"]])))
}
)
stopifnot(all(xy_extract$dist <= max_dist))
print(nc$dim$time$units)
nc_close(nc)

tamsat_dfs <- list()
df <- tidync(f_tamsat)
for(i in seq_len(nrow(xy_extract))) {
  df_tmp <- df %>%
    hyper_filter(lon = lon == xy_extract[i, 1],
                 lat = lat == xy_extract[i, 2]) %>%
    hyper_tibble() %>%
    mutate(date = as.Date(time, origin = as.Date(substr(nc$dim$time$units, 12, 21))))
  df_tmp$station <- xy_extract$station[i]
  df_tmp$req_longitude <- xy_extract$req_longitude[i]
  df_tmp$req_latitude <- xy_extract$req_latitude[i]
  df_tmp$dist <- xy_extract$dist[i]
  tamsat_dfs[[length(tamsat_dfs) + 1]] <- df_tmp
}
westafrica_tamsat <- bind_rows(tamsat_dfs)

westafrica_tamsat <- westafrica_tamsat %>% select(station, date, tamsat_rain = rfe)

saveRDS(westafrica_tamsat, here("data", "station", "cleaned", "westafrica_tamsat.RDS"))

westafrica <- left_join(westafrica, westafrica_tamsat, by = c("station", "date"))

saveRDS(westafrica, here("data", "station", "cleaned", "westafrica_gridded.RDS"))
