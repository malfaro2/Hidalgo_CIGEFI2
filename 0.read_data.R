library(R.matlab)
library(tidyverse)
library(fs)
library(here)

data_dir <- "data/datos_1970_1999"
latlon<-read.table("data/datos_1970_1999/latlon_1970_1999.txt")
latlon$station <- 1:199
names(latlon) <- c("lat","lon","station")
latlon$station <- factor(latlon$station, levels=c(1:199))

dat <- data_dir %>% 
        dir_ls(regexp = "\\.mat$") %>% 
        map_dfr(readMat, .id = "source")  %>% 
        mutate(source = as.numeric(str_extract(source, "(\\d+)\\."))) %>%
        rename(station = source) %>% arrange(station) %>% 
        mutate(month = rep(1:360, 199), lat=rep(latlon$lat,each=360), 
               lon=rep(latlon$lon, each=360))
        
glimpse(dat)

## 199 stations, 12, months, lat, lon, 4 monthly variables, 
## 360 months = 71640 obs.

## Add runoff:

load("data/hidro.Rdata")
data.esco <- sapply(1:199, function(i)apply(fs[[i]],1,median))
data.esco <- as_tibble(data.esco)
data.esco.m <- data.esco %>% 
        mutate(year=1:30) %>% 
pivot_longer(-year,names_to = "station", values_to = "runoff") %>% 
        mutate(station=as.double(substr(station, 2, 100)))  
       

save(dat,data.esco.m, file="data/data.Rdata")

