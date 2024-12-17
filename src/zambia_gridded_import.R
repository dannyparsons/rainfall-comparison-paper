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

zambia <- readRDS(here("data", "station", "cleaned", "zambia_1979_update_qc.RDS"))
zambia_metadata <- readRDS(here("data", "station", "processed", "zambia_stations_metadata_updated.RDS"))

# station location check --------------------------------------------------
zambia_sf <- ne_countries(country = "Zambia", returnclass = "sf")
ggplot(zambia_sf) + 
  geom_sf() +
  geom_point(data = zambia_metadata, aes(x = longitude, y = latitude)) +
  geom_text_repel(data = zambia_metadata, aes(x = longitude, y = latitude, label = station))


# ARC2 Import -------------------------------------------------------------

files <- list.files("D:/data/arc2", pattern = "arc2_zambia_all", full.names = TRUE)

nc <- nc_open(files[1])
stopifnot(min(nc$dim$X$vals) < min(zambia_metadata$longitude))
stopifnot(max(nc$dim$X$vals) > max(zambia_metadata$longitude))
stopifnot(min(nc$dim$Y$vals) < min(zambia_metadata$latitude))
stopifnot(max(nc$dim$Y$vals) > max(zambia_metadata$latitude))
xs <- nc$dim$X$vals
ys <- nc$dim$Y$vals
resx <- xs[2] - xs[1]
resy <- ys[2] - ys[1]
max_dist <- sqrt((resx/2)^2 + (resy/2)^2)
xy_points <- expand.grid(xs, ys)
xy_extract <- closest_point(points = xy_points, 
                  target = zambia_metadata %>% select(longitude, latitude))
xy_extract$station <- zambia_metadata$station
xy_extract$req_longitude <- zambia_metadata$longitude
xy_extract$req_latitude <- zambia_metadata$latitude
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
zambia_arc2 <- bind_rows(arc_dfs) %>%
  mutate(date = as.Date(`T`, origin = as.Date("1960/01/01")))

zambia_arc2 <- zambia_arc2 %>% select(station, date, arc2_rain = est_prcp)

saveRDS(zambia_arc2, here("data", "station", "cleaned", "zambia_arc2.RDS"))

zambia <- left_join(zambia, zambia_arc2, by = c("station", "date"))

# CHIRPS Import -----------------------------------------------------------

files <- list.files("D:/data/chirps", pattern = "chirps_zambia", full.names = TRUE)

nc <- nc_open(files[1])
stopifnot(min(nc$dim$X$vals) < min(zambia_metadata$longitude))
stopifnot(max(nc$dim$X$vals) > max(zambia_metadata$longitude))
stopifnot(min(nc$dim$Y$vals) < min(zambia_metadata$latitude))
stopifnot(max(nc$dim$Y$vals) > max(zambia_metadata$latitude))
xs <- nc$dim$X$vals
ys <- nc$dim$Y$vals
resx <- xs[2] - xs[1]
resy <- ys[2] - ys[1]
max_dist <- sqrt((resx/2)^2 + (resy/2)^2)
xy_points <- expand.grid(xs, ys)
xy_extract <- closest_point(points = xy_points, 
                            target = zambia_metadata %>% select(longitude, latitude))
xy_extract$station <- zambia_metadata$station
xy_extract$req_longitude <- zambia_metadata$longitude
xy_extract$req_latitude <- zambia_metadata$latitude
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
zambia_chirps <- bind_rows(chirps_dfs) %>%
  mutate(date = as.Date(`T`, origin = structure(-2440588, class = "Date")))

zambia_chirps <- zambia_chirps %>% select(station, date, chirps_rain = prcp)

saveRDS(zambia_chirps, here("data", "station", "cleaned", "zambia_chirps.RDS"))

zambia <- left_join(zambia, zambia_chirps, by = c("station", "date"))

# ERA5 Import -------------------------------------------------------------

files <- list.files("D:/data/era5", pattern = "africa-total_precipitation", full.names = TRUE)

nc <- nc_open(files[1])
stopifnot(min(nc$dim$longitude$vals) < min(zambia_metadata$longitude))
stopifnot(max(nc$dim$longitude$vals) > max(zambia_metadata$longitude))
stopifnot(min(nc$dim$latitude$vals) < min(zambia_metadata$latitude))
stopifnot(max(nc$dim$latitude$vals) > max(zambia_metadata$latitude))
xs <- nc$dim$longitude$vals
ys <- nc$dim$latitude$vals
resx <- xs[2] - xs[1]
resy <- ys[2] - ys[1]
max_dist <- sqrt((resx/2)^2 + (resy/2)^2)
xy_points <- expand.grid(xs, ys)
xy_extract <- closest_point(points = xy_points, 
                            target = zambia_metadata %>% select(longitude, latitude))
xy_extract$station <- zambia_metadata$station
xy_extract$req_longitude <- zambia_metadata$longitude
xy_extract$req_latitude <- zambia_metadata$latitude
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
zambia_era5 <- bind_rows(era5_dfs) %>%
  mutate(date_time = as.POSIXct(time * 60 * 60,
                                origin = as.POSIXct("1900/01/01 00:00:00"), tz = "GMT"),
         date = as.Date(date_time),
         # Using UTC +6 as the standard starting day time used by other products.
         # The value at 7am is for the period 6am-7am therefore 7am is the first recording.
         # This is equivalent to 8am rain recording time in Zambia, (Zambia time is UTC + 2)
         tp_lead_7 = dplyr::lead(tp, 7)) %>%
  group_by(station, date) %>%
  summarise(tp = sum(tp_lead_7) * 1000)

zambia_era5 <- zambia_era5 %>% select(station, date, era5_rain = tp)

saveRDS(zambia_era5, here("data", "station", "cleaned", "zambia_era5.RDS"))

zambia <- left_join(zambia, zambia_era5, by = c("station", "date"))


# ERA5 Land Import --------------------------------------------------------

files <- list.files("D:/data/era5_land", pattern = "africa-total_precipitation", full.names = TRUE)

nc <- nc_open(files[1])
stopifnot(min(nc$dim$longitude$vals) < min(zambia_metadata$longitude))
stopifnot(max(nc$dim$longitude$vals) > max(zambia_metadata$longitude))
stopifnot(min(nc$dim$latitude$vals) < min(zambia_metadata$latitude))
stopifnot(max(nc$dim$latitude$vals) > max(zambia_metadata$latitude))
xs <- nc$dim$longitude$vals
ys <- nc$dim$latitude$vals
resx <- xs[2] - xs[1]
resy <- ys[2] - ys[1]
max_dist <- sqrt((resx/2)^2 + (resy/2)^2)
xy_points <- expand.grid(xs, ys)
xy_extract <- closest_point(points = xy_points, 
                            target = zambia_metadata %>% select(longitude, latitude))
xy_extract$station <- zambia_metadata$station
xy_extract$req_longitude <- zambia_metadata$longitude
xy_extract$req_latitude <- zambia_metadata$latitude
xy_extract$dist <- apply(xy_extract, 1, function(r) {
  sp::spDistsN1(matrix(as.numeric(c(r[["Var1"]], r[["Var2"]])), ncol = 2),
                as.numeric(c(r[["req_longitude"]], r[["req_latitude"]])))
}
)
stopifnot(all(xy_extract$dist <= max_dist))
print(nc$dim$time$units)
nc_close(nc)

era5land_dfs <- list()
pb <- txtProgressBar(min = 0, max = length(files) * nrow(xy_extract), style = 3)
count <- 0
for(f in files) {
  print(f)
  for(i in seq_len(nrow(xy_extract))) {
    df <- tidync(f) %>%
      hyper_filter(longitude = longitude == xy_extract[i, 1], 
                   latitude = latitude == xy_extract[i, 2]) %>%
      hyper_tibble()
    df$station <- xy_extract$station[i]
    df$req_longitude <- xy_extract$req_longitude[i]
    df$req_latitude <- xy_extract$req_latitude[i]
    df$dist <- xy_extract$dist[i]
    era5land_dfs[[length(era5land_dfs) + 1]] <- df
    count <- count + 1
    setTxtProgressBar(pb, count)
  }
}
zambia_era5land_hour <- bind_rows(era5land_dfs) %>%
  mutate(date_time = as.POSIXct(time * 60 * 60,
                                origin = as.POSIXct("1900/01/01 00:00:00"), tz = "GMT"),
         date = as.Date(date_time),
         tp = tp * 1000)

saveRDS(zambia_era5land_hour, here("data", "station", "cleaned", "zambia_era5land_hourly.RDS"))

View(zambia_era5land_hour %>% filter(station == "Livingstone" & year(date) == 2007 & month(date) %in% c(11, 12)))

# Using UTC +6 as the standard starting day time used by other products.
# ERA5 Land recordings are accumulations since 00:00
# Value at 00:00 is the 24 hour accumulation
# Therefore, the 6am-6am total rainfall on day D = (06:00 day D+1 value) + (00:00 day D+1 value) - (06:00 day D+1 value)
# This is equivalent to 8am rain recording time in Zambia, (Zambia time is UTC + 2)
zambia_era5land <- zambia_era5land_hour %>%
  # important to group by station so that lags are correct
  group_by(station) %>%
  mutate(tp_lead_6 = dplyr::lead(tp, 6),
         tp_lag_18 = dplyr::lag(tp, 18)) %>%
  group_by(station, date) %>%
  summarise(tp = dplyr::first(tp) + dplyr::first(tp_lead_6) - dplyr::first(tp_lag_18)) %>%
  mutate(tp = dplyr::lead(tp),
         tp = ifelse(tp < 0, 0, tp))

zambia_era5land <- zambia_era5land %>% select(station, date, era5land_rain = tp)

saveRDS(zambia_era5land, here("data", "station", "cleaned", "zambia_era5land.RDS"))

zambia <- left_join(zambia, zambia_era5land, by = c("station", "date"))


# IMERG Import ------------------------------------------------------------

files <- list.files("D:/data/imerg", pattern = "imerge_africa", full.names = TRUE)

nc <- nc_open(files[1])
stopifnot(min(nc$dim$lon$vals) < min(zambia_metadata$longitude))
stopifnot(max(nc$dim$lon$vals) > max(zambia_metadata$longitude))
stopifnot(min(nc$dim$lat$vals) < min(zambia_metadata$latitude))
stopifnot(max(nc$dim$lat$vals) > max(zambia_metadata$latitude))
xs <- nc$dim$lon$vals
ys <- nc$dim$lat$vals
resx <- xs[2] - xs[1]
resy <- ys[2] - ys[1]
max_dist <- sqrt((resx/2)^2 + (resy/2)^2)
xy_points <- expand.grid(xs, ys)
xy_extract <- closest_point(points = xy_points, 
                            target = zambia_metadata %>% select(longitude, latitude))
xy_extract$station <- zambia_metadata$station
xy_extract$req_longitude <- zambia_metadata$longitude
xy_extract$req_latitude <- zambia_metadata$latitude
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
zambia_imerg <- bind_rows(imerg_dfs) %>%
  mutate(date = as.Date(time, origin = as.Date("1970/01/01")))

zambia_imerg <- zambia_imerg %>% select(station, date, imerg_rain = precipitationCal)

saveRDS(zambia_imerg, here("data", "station", "cleaned", "zambia_imerg.RDS"))

zambia <- left_join(zambia, zambia_imerg, by = c("station", "date"))

# RFE2 Import -------------------------------------------------------------

files <- list.files("D:/data/rfe2", pattern = "zambia", full.names = TRUE)

nc <- nc_open(files[2])
stopifnot(min(nc$dim$X$vals) < min(zambia_metadata$longitude))
stopifnot(max(nc$dim$X$vals) > max(zambia_metadata$longitude))
stopifnot(min(nc$dim$Y$vals) < min(zambia_metadata$latitude))
stopifnot(max(nc$dim$Y$vals) > max(zambia_metadata$latitude))
xs <- nc$dim$X$vals
ys <- nc$dim$Y$vals
resx <- xs[2] - xs[1]
resy <- ys[2] - ys[1]
max_dist <- sqrt((resx/2)^2 + (resy/2)^2)
xy_points <- expand.grid(xs, ys)
xy_extract <- closest_point(points = xy_points, 
                            target = zambia_metadata %>% select(longitude, latitude))
xy_extract$station <- zambia_metadata$station
xy_extract$req_longitude <- zambia_metadata$longitude
xy_extract$req_latitude <- zambia_metadata$latitude
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
zambia_rfe2 <- bind_rows(rfe2_dfs) %>%
  mutate(date = as.Date(`T`, origin = as.Date("2000/10/31")))

zambia_rfe2 <- zambia_rfe2 %>% select(station, date, rfe2_rain = est_prcp)

zambia <- left_join(zambia, zambia_rfe2, by = c("station", "date"))

saveRDS(zambia_rfe2, here("data", "station", "cleaned", "zambia_rfe2.RDS"))

# TAMSAT 3.1 Import -------------------------------------------------------

files <- c()
for(y in 1983:2019) {
  for(m in 1:12) {
    mm <- sprintf("%02d", m)
    files <- c(files, list.files(paste("D:/data/tamsat3.1", y, mm, sep = "/"), 
                                 full.names = TRUE))
  }
}

nc <- nc_open(files[1])
stopifnot(min(nc$dim$lon$vals) < min(zambia_metadata$longitude))
stopifnot(max(nc$dim$lon$vals) > max(zambia_metadata$longitude))
stopifnot(min(nc$dim$lat$vals) < min(zambia_metadata$latitude))
stopifnot(max(nc$dim$lat$vals) > max(zambia_metadata$latitude))
xs <- nc$dim$lon$vals
ys <- nc$dim$lat$vals
resx <- xs[2] - xs[1]
resy <- ys[2] - ys[1]
max_dist <- sqrt((resx/2)^2 + (resy/2)^2)
xy_points <- expand.grid(xs, ys)
xy_extract <- closest_point(points = xy_points, 
                            target = zambia_metadata %>% select(longitude, latitude))
xy_extract$station <- zambia_metadata$station
xy_extract$req_longitude <- zambia_metadata$longitude
xy_extract$req_latitude <- zambia_metadata$latitude
xy_extract$dist <- apply(xy_extract, 1, function(r) {
  sp::spDistsN1(matrix(as.numeric(c(r[["Var1"]], r[["Var2"]])), ncol = 2),
                as.numeric(c(r[["req_longitude"]], r[["req_latitude"]])))
}
)
stopifnot(all(xy_extract$dist <= max_dist))
print(nc$dim$time$units)
nc_close(nc)

tamsat_dfs <- list()
pb <- txtProgressBar(min = 0, max = length(files) * nrow(xy_extract), style = 3)
count <- 0
for(f in files) {
  print(f)
  nc <- nc_open(f)
  df <- tidync(f) 
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
    count <- count + 1
    setTxtProgressBar(pb, count)
  }
  nc_close(nc)
}
zambia_tamsat <- bind_rows(tamsat_dfs)

zambia_tamsat <- zambia_tamsat %>% select(station, date, tamsat_rain = rfe_filled)

saveRDS(zambia_tamsat, here("data", "station", "cleaned", "zambia_tamsat.RDS"))

zambia <- left_join(zambia, zambia_tamsat, by = c("station", "date"))

# TAMSAT 3 Import ---------------------------------------------------------

# files <- list.files("D:/data/tamsat", full.names = TRUE)
# 
# nc <- nc_open(files[1])
# stopifnot(min(nc$dim$lon$vals) < min(zambia_metadata$longitude))
# stopifnot(max(nc$dim$lon$vals) > max(zambia_metadata$longitude))
# stopifnot(min(nc$dim$lat$vals) < min(zambia_metadata$latitude))
# stopifnot(max(nc$dim$lat$vals) > max(zambia_metadata$latitude))
# xs <- nc$dim$lon$vals
# ys <- nc$dim$lat$vals
# resx <- xs[2] - xs[1]
# resy <- ys[2] - ys[1]
# max_dist <- sqrt((resx/2)^2 + (resy/2)^2)
# xy_points <- expand.grid(xs, ys)
# xy_extract <- closest_point(points = xy_points, 
#                             target = zambia_metadata %>% select(longitude, latitude))
# xy_extract$station <- zambia_metadata$station
# xy_extract$req_longitude <- zambia_metadata$longitude
# xy_extract$req_latitude <- zambia_metadata$latitude
# xy_extract$dist <- apply(xy_extract, 1, function(r) {
#   sp::spDistsN1(matrix(as.numeric(c(r[["Var1"]], r[["Var2"]])), ncol = 2),
#                 as.numeric(c(r[["req_longitude"]], r[["req_latitude"]])))
# }
# )
# stopifnot(all(xy_extract$dist <= max_dist))
# print(nc$dim$time$units)
# nc_close(nc)
# 
# tamsat_dfs <- list()
# pb <- txtProgressBar(min = 0, max = length(files) * nrow(xy_extract), style = 3)
# count <- 0
# for(f in files) {
#   print(f)
#   nc <- nc_open(f)
#   df <- tidync(f) 
#   for(i in seq_len(nrow(xy_extract))) {
#     df_tmp <- df %>%
#       hyper_filter(lon = lon == xy_extract[i, 1],
#                    lat = lat == xy_extract[i, 2]) %>%
#       hyper_tibble() %>%
#       mutate(date = as.Date(time, origin = as.Date(substr(nc$dim$time$units, 12, 21))))
#     df_tmp$station <- xy_extract$station[i]
#     df_tmp$req_longitude <- xy_extract$req_longitude[i]
#     df_tmp$req_latitude <- xy_extract$req_latitude[i]
#     df_tmp$dist <- xy_extract$dist[i]
#     tamsat_dfs[[length(tamsat_dfs) + 1]] <- df_tmp
#     count <- count + 1
#     setTxtProgressBar(pb, count)
#   }
#   nc_close(nc)
# }
# zambia_tamsat <- bind_rows(tamsat_dfs)
# 
# zambia_tamsat <- zambia_tamsat %>% select(station, date, tamsat_rain = rfe)
# 
# saveRDS(zambia_tamsat, here("data", "station", "cleaned", "zambia_tamsat.RDS"))
# 
# zambia <- left_join(zambia, zambia_tamsat, by = c("station", "date"))
# 
# saveRDS(zambia, here("data", "station", "cleaned", "zambia_gridded.RDS"))


# ENACTS Import -----------------------------------------------------------

files <- list.files(here("data", "enacts", "zambia", "rain"), full.names = TRUE)
names(files) <- list.files(here("data", "enacts", "zambia", "rain"), full.names = FALSE)
nc <- nc_open(files[1])
stopifnot(min(nc$dim$Lon$vals) < min(zambia_metadata$longitude))
stopifnot(max(nc$dim$Lon$vals) > max(zambia_metadata$longitude))
stopifnot(min(nc$dim$Lat$vals) < min(zambia_metadata$latitude))
stopifnot(max(nc$dim$Lat$vals) > max(zambia_metadata$latitude))
xs <- nc$dim$Lon$vals
ys <- nc$dim$Lat$vals
resx <- xs[2] - xs[1]
resy <- ys[2] - ys[1]
max_dist <- sqrt((resx/2)^2 + (resy/2)^2)
xy_points <- expand.grid(xs, ys)
xy_extract <- closest_point(points = xy_points, 
                            target = zambia_metadata %>% select(longitude, latitude))
xy_extract$station <- zambia_metadata$station
xy_extract$req_longitude <- zambia_metadata$longitude
xy_extract$req_latitude <- zambia_metadata$latitude
xy_extract$dist <- apply(xy_extract, 1, function(r) {
  sp::spDistsN1(matrix(as.numeric(c(r[["Var1"]], r[["Var2"]])), ncol = 2),
                as.numeric(c(r[["req_longitude"]], r[["req_latitude"]])))
}
)
stopifnot(all(xy_extract$dist <= max_dist))
nc_close(nc)

enacts_dfs <- list()
pb <- txtProgressBar(min = 0, max = length(files) * nrow(xy_extract), style = 3)
count <- 0
for(j in seq_along(files)) {
  name <- substr(names(files)[j], 8, 15)
  for(i in seq_len(nrow(xy_extract))) {
    s <- xy_extract$station[i]
    df <- tidync(files[[j]]) %>%
      hyper_filter(Lon = Lon == xy_extract[i, 1], 
                   Lat = Lat == xy_extract[i, 2]) %>%
      hyper_tibble()
    df$station <- s
    df$req_latitude <- xy_extract$req_latitude[i]
    df$req_longitude <- xy_extract$req_longitude[i]
    df$dist <- xy_extract$dist[i]
    enacts_dfs[[paste(name, s, sep = "_")]] <- df
    count <- count + 1
    setTxtProgressBar(pb, count)
  }
}
zambia_enacts <- bind_rows(enacts_dfs, .id = "date") %>%
  mutate(date = as.Date(substr(date, 1, 8), format = "%Y%m%d"))

zambia_enacts <- zambia_enacts %>% select(station, date, enacts_rain = precip)

saveRDS(zambia_enacts, here("data", "station", "cleaned", "zambia_enacts.RDS"))

zambia <- left_join(zambia, zambia_enacts, by = c("station", "date"))


# ENACTS 2019 Import ------------------------------------------------------

files <- list.files(here("data", "enacts", "zambia", "rain_2019"), full.names = TRUE)
names(files) <- list.files(here("data", "enacts", "zambia", "rain_2019"), full.names = FALSE)
nc <- nc_open(files[1])
stopifnot(min(nc$dim$Lon$vals) < min(zambia_metadata$longitude))
stopifnot(max(nc$dim$Lon$vals) > max(zambia_metadata$longitude))
stopifnot(min(nc$dim$Lat$vals) < min(zambia_metadata$latitude))
stopifnot(max(nc$dim$Lat$vals) > max(zambia_metadata$latitude))
xs <- nc$dim$Lon$vals
ys <- nc$dim$Lat$vals
resx <- xs[2] - xs[1]
resy <- ys[2] - ys[1]
max_dist <- sqrt((resx/2)^2 + (resy/2)^2)
xy_points <- expand.grid(xs, ys)
xy_extract <- closest_point(points = xy_points, 
                            target = zambia_metadata %>% select(longitude, latitude))
xy_extract$station <- zambia_metadata$station
xy_extract$req_longitude <- zambia_metadata$longitude
xy_extract$req_latitude <- zambia_metadata$latitude
xy_extract$dist <- apply(xy_extract, 1, function(r) {
  sp::spDistsN1(matrix(as.numeric(c(r[["Var1"]], r[["Var2"]])), ncol = 2),
                as.numeric(c(r[["req_longitude"]], r[["req_latitude"]])))
}
)
stopifnot(all(xy_extract$dist <= max_dist))
nc_close(nc)

enacts_dfs <- list()
pb <- txtProgressBar(min = 0, max = length(files) * nrow(xy_extract), style = 3)
count <- 0
for(j in seq_along(files)) {
  name <- substr(names(files)[j], 8, 15)
  for(i in seq_len(nrow(xy_extract))) {
    s <- xy_extract$station[i]
    df <- tidync(files[[j]]) %>%
      hyper_filter(Lon = Lon == xy_extract[i, 1], 
                   Lat = Lat == xy_extract[i, 2]) %>%
      hyper_tibble()
    df$station <- s
    df$req_latitude <- xy_extract$req_latitude[i]
    df$req_longitude <- xy_extract$req_longitude[i]
    df$dist <- xy_extract$dist[i]
    enacts_dfs[[paste(name, s, sep = "_")]] <- df
    count <- count + 1
    setTxtProgressBar(pb, count)
  }
}
zambia_enacts_2019 <- bind_rows(enacts_dfs, .id = "date") %>%
  mutate(date = as.Date(substr(date, 1, 8), format = "%Y%m%d"))

zambia_enacts_2019 <- zambia_enacts_2019 %>% select(station, date, enacts_rain = precip)

saveRDS(zambia_enacts_2019, here("data", "station", "cleaned", "zambia_enacts_2019.RDS"))

zambia <- left_join(zambia, zambia_enacts_2019, by = c("station", "date"))

# Save all ----------------------------------------------------------------

saveRDS(zambia, here("data", "station", "cleaned", "zambia_gridded.RDS"))
