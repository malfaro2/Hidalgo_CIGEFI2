## Data and packages

source(file="0.packages.R")
source(file="0.functions.R")
load(file="data/data.Rdata")
estaciones <- read.table("data/datos_1970_1999/latlon_1970_1999.txt", header=F)
estaciones$station <- 1:199
names(estaciones) <- c("lat","lon","station")
estaciones$station <- factor(estaciones$station, levels=c(1:199))
meteo <- dat
map1 <- rnaturalearth::ne_states(
  country = c("guatemala", "honduras", 
              "el salvador", "panama", 
              "nicaragua", "costa rica", 
              "belize", "mexico", 
              "colombia"), returnclass = "sf")


## Center all variables:

meteo <- meteo %>% mutate(year=ceiling(meteo$month/12)) %>% 
  group_by(year,station) %>% 
  summarize(pmensual=mean(pmensual),
            Tavgmensual=mean(Tavgmensual),
            PET=mean(PET),aridity=mean(aridity)) %>% 
  left_join(data.esco.m, by=c("station","year")) %>% 
  group_by(station) %>% 
  mutate(tavg.m=mean(Tavgmensual,na.rm=TRUE),
         ari.m=mean(aridity,na.rm=TRUE),
         pet.m=mean(PET,na.rm=TRUE),
         ro.m=mean(runoff,na.rm=TRUE),
         prec.m=mean(pmensual,na.rm=TRUE)) %>% 
  ungroup() %>% 
  mutate(tavg.c=Tavgmensual-tavg.m,
         ari.c=aridity-ari.m,
         pet.c=PET-pet.m,
         ro.c=runoff-ro.m,
         prec.c=pmensual-prec.m)

## Variable Description - Boxplots

a<-meteo %>%
  mutate(Year = year+69) %>%
  ggplot(aes(x=as.factor(Year), y=aridity)) +
  geom_boxplot() +
  theme_ipsum() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)) +
  xlab("Year") +
  ylab("Mean Aridity") 

b<-meteo %>%
  mutate(Year = year+69) %>%
  ggplot(aes(x=as.factor(Year), y=runoff)) +
  geom_boxplot() +
  theme_ipsum() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)) +
  xlab("Year") +
  ylab("Mean Runoff")

c<-meteo %>%
  mutate(Year = year+69) %>%
  ggplot(aes(x=as.factor(Year), y=PET)) +
  geom_boxplot() +
  theme_ipsum() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)) +
  xlab("Year") +
  ylab("Mean PET")

d<-meteo %>%
  mutate(Year = year+69) %>%
  ggplot(aes(x=as.factor(Year), y=pmensual)) +
  geom_boxplot() +
  theme_ipsum() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)) +
  xlab("Year") +
  ylab("Mean Precipitation")

e<-meteo %>%
  mutate(Year = year+69) %>%
  ggplot(aes(x=as.factor(Year), y=Tavgmensual)) +
  geom_boxplot() +
  theme_ipsum() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)) +
  xlab("Year") +
  ylab("Mean Temperature")

## Trends - Only calculates trends

###### variable 1: pmensual #######
trend.pmensual.c <-meteo %>% arrange(station) %>% 
  mutate(station = factor(station, levels=c(1:199))) %>% 
  nest(-station) %>% 
  mutate(coef = purrr::map(data, ~ 
                lm(prec.c ~ year-1,data = .)$coef)) %>% 
  unnest(coef) %>% left_join(estaciones, by=c("station")) %>% 
  dplyr::select(-data)

###### variable 2: PET #######
trend.PET.c <-meteo %>% arrange(station) %>% 
  mutate(station = factor(station, levels=c(1:199))) %>% 
  nest(-station) %>% 
  mutate(coef = purrr::map(data, ~ 
                lm(pet.c ~ year-1,data = .)$coef)) %>% 
  unnest(coef) %>% left_join(estaciones, by=c("station")) %>% 
  dplyr::select(-data)

###### variable 3: aridity #######
trend.aridity.c <-meteo %>% arrange(station) %>% 
  mutate(station = factor(station, levels=c(1:199))) %>% 
  nest(-"station") %>% 
  mutate(coef = purrr::map(data, ~ 
                lm(ari.c ~ year-1,data = .)$coef)) %>% 
  tidyr::unnest(coef) %>% left_join(estaciones, by=c("station")) %>% 
  dplyr::select(-data)

###### variable 4: tavg #######
trend.Tavgmensual.c <-meteo %>% arrange(station) %>% 
  mutate(station = factor(station, levels=c(1:199))) %>% 
  nest(-station) %>% 
  mutate(coef = purrr::map(data, ~ 
                lm(tavg.c ~ year-1,data = .)$coef)) %>% 
  unnest(coef) %>% left_join(estaciones, by=c("station")) %>% 
  dplyr::select(-data)

###### variable 5: runoff #######
trend.runoff.c <-meteo %>% arrange(station) %>% 
  mutate(station = factor(station, levels=c(1:199))) %>% 
  nest(-station) %>% 
  mutate(coef = purrr::map(data, ~ lm(ro.c ~ year-1,data = .)$coef)) %>% 
  unnest(coef) %>% left_join(estaciones, by=c("station")) %>% 
  dplyr::select(-data)

## Exploratories (not included in the paper)

dotchart((sort(trend.pmensual.c$coef)), 
         labels=trend.pmensual.c$station[order(
           trend.pmensual.c$coef)], 
         cex=0.65, main="Coefficient per station")

dotchart((sort(trend.PET.c$coef)), 
         labels=trend.PET.c$station[order(
           trend.PET.c$coef)], 
         cex=0.65,main="Coefficient per station")

dotchart((sort(trend.aridity.c$coef)), 
         labels=trend.aridity.c$station[order(
           trend.aridity.c$coef)], 
         cex=0.65,main="Coefficient per station")

dotchart((sort(trend.Tavgmensual.c $coef)), 
         labels=trend.Tavgmensual.c$station[order(
           trend.Tavgmensual.c $coef)], 
         cex=0.65,main="Coefficient per station")

dotchart((sort(trend.runoff.c$coef)), 
         labels=trend.runoff.c$station[order(
           trend.runoff.c $coef)], 
         cex=0.65,main="Coefficient per station")

## Arrange all trends in one tibble:

trends <- trend.pmensual.c %>% 
  inner_join(trend.PET.c, by=c("station"="station",
            "lat"="lat","lon"="lon")) %>% 
  inner_join(trend.aridity.c, by=c("station"="station",
            "lat"="lat","lon"="lon")) %>% 
  inner_join(trend.Tavgmensual.c, by=c("station"="station",
            "lat"="lat","lon"="lon")) %>% 
  inner_join(trend.runoff.c, by=c("station"="station",
            "lat"="lat","lon"="lon"))

names(trends) <- c("station","trend.pm.c","lat","lon",
                   "trend.PET.c","trend.aridity.c",
                   "trend.Tavgmensual.c","trend.runoff.c")

trends <- trends[,c("station","lon","lat","trend.pm.c",
                    "trend.PET.c","trend.aridity.c",
                    "trend.Tavgmensual.c","trend.runoff.c" )]

# Figures and tables for the paper:

# Table 1
summary(meteo[,c(3:7)])
# Figure 1
d
# Figure 2
get.map_descriptives(trend.pmensual.c)
# Figure A1
grid.arrange(a,b,c,e,nrow=2)
# Figure A2
get.maps_descriptives(trends)

## Save data for analysis:

save(trends, file="data/trends.Rdata")
save(meteo,estaciones,data.esco.m, file="data/data_for_regression.Rdata")
