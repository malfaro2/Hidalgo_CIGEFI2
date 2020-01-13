source(file="0.packages.R")
load(file="data/data.Rdata")
estaciones <- read.table("data/datos_1970_1999/latlon_1970_1999.txt", header=F)
estaciones$station <- 1:199
names(estaciones) <- c("lat","lon","station")
estaciones$station <- factor(estaciones$station, levels=c(1:199))
meteo <- dat
map <- rnaturalearth::ne_states(country = c("guatemala", "honduras", "el salvador", "panama", "nicaragua", "costa rica", "belize"), returnclass = "sf")

meteo <- meteo %>% mutate(year=ceiling(meteo$month/12)) %>% 
  group_by(year,station) %>% 
  summarize(pmensual=mean(pmensual),
            Tavgmensual=mean(Tavgmensual),
            PET=mean(PET),aridity=mean(aridity), 
            lat=unique(lat), lon=unique(lon)) %>% 
  left_join(data.esco.m, by=c("station","year")) %>% 
  group_by(station) 
n <- length(unique(meteo$year))
fps <- 20
duration <- 80

map_prep <- meteo %>% 
  mutate(year=year+1969) %>% 
  ggplot(aes(lon, lat)) +
  geom_sf(data = map, inherit.aes = FALSE) +
  coord_sf(ylim = c(7,18), xlim = c(-78, -93)) +
  geom_point(aes(fill = runoff, size = runoff), shape = 21) +
  scale_fill_viridis(name="Runoff",discrete=FALSE) +
  theme_ipsum() + theme(
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 12)) +
  scale_size_area(max_size = 6, guide = "none") +
  labs(title = "Central America: Runoff for {round(frame_time)}", x = "Lon", y = "Lat") +
  transition_time(year)

map_anom <- gganimate::animate(map_prep, 
                               nframes = n, 
                               duration = duration, 
                    width = 480, height = 340)

serie <- meteo %>% 
  mutate(year=year+1969) %>% 
  group_by(year) %>% 
  ggplot(aes(year, runoff)) +
  geom_line(stat = "summary", fun.y = mean) +
  scale_x_continuous("Year") +
  scale_y_continuous("Mean Runoff") +
  labs(title = "Annual Overall Mean") + 
  theme_ipsum() + theme(
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 12)) +
  transition_reveal(year)

serie <- gganimate::animate(serie, 
                            nframes = n, 
                            duration = duration, 
                 width = 480, height = 170)
 
a_mgif <- image_read(map_anom)
b_mgif <- image_read(serie)
new_gif <- image_append(c(a_mgif[1], b_mgif[1]), stack = TRUE)
for(i in 2:n){
  combined <- image_append(c(a_mgif[i], b_mgif[i]), stack = TRUE)
  new_gif <- c(new_gif, combined)
}
save_animation(new_gif, "map-time-runoff.gif")

