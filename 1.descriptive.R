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
  ylab("Mean Aridity")  +
  labs(caption="units: mm/mm")

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
  ylab("Mean Runoff") +
  labs(caption="units: mm/day")

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
  ylab("Mean PET") +
  labs(caption="units: mm/month")

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
  ylab("Mean Precipitation") +
  labs(caption="units: mm/month")

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
  ylab("Mean Temperature") +
  labs(caption="units: Celsius")

## Trends - Only calculates trends
## 

###### variable 1: pmensual #######
trend.pm.c <- meteo %>% arrange(station) %>% 
  mutate(station = factor(station, levels=c(1:199))) %>% 
  nest(data = c(year, pmensual, Tavgmensual, 
                PET, aridity, runoff, tavg.m, 
                ari.m, pet.m, ro.m, prec.m, 
                tavg.c,ari.c, pet.c, ro.c, prec.c)) %>% 
  mutate(coef = purrr::map(data, ~ 
                lm(prec.c ~ year-1,data = .)$coef),
         ci.b =   purrr::map(data, ~ 
                boot_CI(.,prec.c~year-1)),
         p.mk =  purrr::map(data, ~
  as.numeric(MannKendall(ts(.$prec.c))$sl))) %>% 
  tidyr::unnest(coef) %>%  tidyr::unnest_wider(ci.b) %>% 
  tidyr::unnest(p.mk) %>% 
  dplyr::select(station,coef,linf,lsup,p.mk) %>% 
  left_join(estaciones, by=c("station"))

###### variable 2: PET #######
trend.pet.c <-meteo %>% arrange(station) %>% 
  mutate(station = factor(station, levels=c(1:199))) %>% 
  nest(data = c(year, pmensual, Tavgmensual, 
                PET, aridity, runoff, tavg.m, 
                ari.m, pet.m, ro.m, prec.m, 
                tavg.c,ari.c, pet.c, ro.c, prec.c)) %>% 
  mutate(coef = purrr::map(data, ~ 
        lm(pet.c ~ year-1,data = .)$coef),
        ci.b =   purrr::map(data, ~ 
                boot_CI(.,pet.c~year-1)),
        p.mk =  purrr::map(data, ~
        as.numeric(MannKendall(ts(.$pet.c))$sl))) %>% 
  tidyr::unnest(coef) %>%  tidyr::unnest_wider(ci.b) %>% 
  tidyr::unnest(p.mk) %>% 
  dplyr::select(station,coef,linf,lsup,p.mk) %>% 
  left_join(estaciones, by=c("station"))

###### variable 3: aridity #######
trend.ari.c <-meteo %>% arrange(station) %>% 
  mutate(station = factor(station, levels=c(1:199))) %>% 
  nest(data = c(year, pmensual, Tavgmensual, 
                PET, aridity, runoff, tavg.m, 
                ari.m, pet.m, ro.m, prec.m, 
                tavg.c,ari.c, pet.c, ro.c, prec.c)) %>% 
  mutate(coef = purrr::map(data, ~ 
            lm(ari.c ~ year-1,data = .)$coef),
         ci.b =   purrr::map(data, ~ 
         boot_CI(.,ari.c~year-1)),
         p.mk =  purrr::map(data, ~
        as.numeric(MannKendall(ts(.$ari.c))$sl))) %>% 
  tidyr::unnest(coef) %>%  tidyr::unnest_wider(ci.b) %>% 
  tidyr::unnest(p.mk) %>% 
  dplyr::select(station,coef,linf,lsup,p.mk) %>% 
  left_join(estaciones, by=c("station"))

###### variable 4: tavg #######
trend.temp.c <-meteo %>% arrange(station) %>% 
  mutate(station = factor(station, levels=c(1:199))) %>% 
  nest(data = c(year, pmensual, Tavgmensual, 
                PET, aridity, runoff, tavg.m, 
                ari.m, pet.m, ro.m, prec.m, 
                tavg.c,ari.c, pet.c, ro.c, prec.c)) %>% 
  mutate(coef = purrr::map(data, ~ 
                 lm(tavg.c ~ year-1,data = .)$coef),
         ci.b =   purrr::map(data, ~ 
          boot_CI(.,tavg.c~year-1)),
         p.mk =  purrr::map(data, ~
          as.numeric(MannKendall(ts(.$tavg.c))$sl))) %>% 
  tidyr::unnest(coef) %>%  tidyr::unnest_wider(ci.b) %>% 
  tidyr::unnest(p.mk) %>% 
  dplyr::select(station,coef,linf,lsup,p.mk) %>% 
  left_join(estaciones, by=c("station"))

###### variable 5: runoff #######
trend.ro.c <-meteo %>% arrange(station) %>% 
  mutate(station = factor(station, levels=c(1:199))) %>% 
  nest(data = c(year, pmensual, Tavgmensual, 
                PET, aridity, runoff, tavg.m, 
                ari.m, pet.m, ro.m, prec.m, 
                tavg.c,ari.c, pet.c, ro.c, prec.c)) %>% 
  mutate(coef = purrr::map(data, ~ 
                lm(ro.c ~ year-1,data = .)$coef),
         ci.b =   purrr::map(data, ~ 
             boot_CI(.,ro.c~year-1)),
         p.mk =  purrr::map(data, ~
             as.numeric(MannKendall(ts(.$ro.c))$sl))) %>% 
  tidyr::unnest(coef) %>%  tidyr::unnest_wider(ci.b) %>% 
  tidyr::unnest(p.mk) %>% 
  dplyr::select(station,coef,linf,lsup,p.mk) %>% 
  left_join(estaciones, by=c("station"))

## Exploratories (not included in the paper)
## outliers: 

data <- trend.pm.c
## outliers: more than 3 mads away from median.
lsup<-median(data$coef)+5*mad(data$coef)
linf<-median(data$coef)-5*mad(data$coef)
out.pm <- data[which(data$coef<linf|data$coef>lsup),]
title <- "Trend Coefficients Precipitation"
units <- "units: change in mm/month per year"
tt <- "(A)"
aa<-get_dots(data,units,title,tt,out.pm)

data <- trend.pet.c
lsup<-median(data$coef)+5*mad(data$coef)
linf<-median(data$coef)-5*mad(data$coef)
out.pet <- data[which(data$coef<linf|data$coef>lsup),]
title <- "Trend Coefficients PET"
units <- "units: change in mm/month per year"
tt <- "(D)"
bb<-get_dots(data,units,title,tt,out.pet)

data <- trend.ari.c
lsup<-median(data$coef)+5*mad(data$coef)
linf<-median(data$coef)-5*mad(data$coef)
out.ari<- data[which(data$coef<linf|data$coef>lsup),]
title <- "Trend Coefficients Aridity"
units <- "units: change in mm/mm per year"
tt <- "(C)"
cc<-get_dots(data,units,title,tt,out.ari)

data <- trend.temp.c
lsup<-median(data$coef)+5*mad(data$coef)
linf<-median(data$coef)-5*mad(data$coef)
out.temp <- data[which(data$coef<linf|data$coef>lsup),]
title <- "Trend Coefficients Temperature"
units <- "units: change in Celsius per year"
tt <- "(E)"
dd<-get_dots(data,units,title,tt,out.temp)

data <- trend.ro.c
lsup<-median(data$coef)+5*mad(data$coef)
linf<-median(data$coef)-5*mad(data$coef)
out.ro <- data[which(data$coef<linf|data$coef>lsup),]
title <- "Trend Coefficients Runoff"
units <- "units: change in mm/day per year"
tt <- "(B)"
ee<-get_dots(data,units,title,tt,out.ro)

rm(data)
## Arrange all trends in one tibble:

trends <- trend.pm.c %>% 
  inner_join(trend.pet.c, by=c("station"="station",
            "lat"="lat","lon"="lon")) %>% 
  inner_join(trend.ari.c, by=c("station"="station",
            "lat"="lat","lon"="lon")) %>% 
  inner_join(trend.temp.c, by=c("station"="station",
            "lat"="lat","lon"="lon")) %>% 
  inner_join(trend.ro.c, by=c("station"="station",
            "lat"="lat","lon"="lon"))

names(trends) <- c("station",
      "trend.pm.c","li.trend.pm.c",
      "ls.trend.pm.c","mk.trend.pm.c",
                   "lat","lon",
      "trend.pet.c","li.trend.pet.c",
      "ls.trend.pet.c","mk.trend.pet.c",
      "trend.ari.c","li.trend.ari.c",
      "ls.trend.ari.c","mk.trend.ari.c",
      "trend.temp.c","li.trend.temp.c",
      "ls.trend.temp.c","mk.trend.temp.c",
      "trend.ro.c","li.trend.ro.c",
      "ls.trend.ro.c","mk.trend.ro.c")

# Figures and tables for the paper:

# Table 1
summary(meteo[,c(3:7)])
# Figure 1
d 
# Figure A1
grid.arrange(a,b,c,e,nrow=2) # agregar unidad de medida
# Figure 3
grid.arrange(aa,ee,cc,bb,dd,nrow=2) 

## Save data for analysis:
## 
outliers <- list(out.ari,out.pet,out.pm,
                 out.ro, out.temp)

save(trends,outliers, file="data/trends.Rdata")
save(meteo,estaciones,data.esco.m, file="data/data_for_regression.Rdata")
