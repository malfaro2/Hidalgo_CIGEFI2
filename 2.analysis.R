## Data and packages

load(file="data/trends.Rdata")

map1 <- rnaturalearth::ne_states(
  country = c("guatemala", "honduras", 
              "el salvador", "panama", 
              "nicaragua", "costa rica", 
              "belize", "mexico", 
              "colombia"), returnclass = "sf")

all_out <- unique(c(outliers[[1]]$station,
outliers[[3]]$station,outliers[[4]]$station))

## Set coordinates and projections:
trend_WOO <- trends[!trends$station%in%all_out,]
trend_WO  <- trends
coordinates(trend_WOO) <- ~lon+lat
trends <- trend_WOO
n <- dim(trends)[1]
n0 <- dim(trend_WO)[1]
source(file="0.packages.R")
source(file="0.functions.R")
print(c(n,n0))

## Estimate variograms for each variable

###### variable 1: aridity #######

var=variogram(trends$trend.ari.c~1, data=trends,cutoff=10)
plot(var)
varmodel <- vgm(0.0000025,"Gau",10,.0000005)
plot(var, model=varmodel) 
fitmodel.3 <- fit.variogram(var, model=varmodel, fit.sills=FALSE)
mean(trends$trend.ari.c)
fitmodel.3$psill[1]/fitmodel.3$psill[2]
fitmodel.3
# model   psill    range
# 1   Nug 5.0e-07 0.000000
# 2   Gau 2.5e-06 6.259119

weig<-var$np/var$dist^2  
SSTot<- sum((weig*(var$gamma-weighted.mean(var$gamma, weig)))^2)
SSErr <- attr(fitmodel.3,"SSErr")
R2 <- 1-(SSErr)/SSTot;R2

print(xyplot(gamma ~ dist, var, pch = 3, 
             type = 'b', lwd = 2, col = 'darkblue',
             panel = function(x, y, ...) {
               for (i in 1:100) {
                 trends$random = sample(trends$trend.ari.c)
                 v = variogram((random) ~ 1, trends,cutoff=10)
                 llines(v$dist, v$gamma, col = 'grey')
               }
               panel.xyplot(x, y, ...)
             },
             ylim = c(0, 0.0000025), 
             xlab = 'distance', ylab = 'semivariance'
))

###### variable 2: PET #######

var=variogram(trends$trend.pet.c~1, data=trends)
plot(var)
varmodel <- vgm(psill=0.005,model="Gau",5,nugget=0)
varmodel2 <- vgm(psill=0.005,model="Exp",5,nugget=0)
par(mfrow=c(1,2))
plot(var, model=varmodel, main="PET trends - Gaussian")
plot(var, model=varmodel2, main="PET trends - Exponential")
fitmodel.2 <- fit.variogram(var, model=varmodel, fit.ranges=FALSE)
mean(trends$trend.pet.c)
fitmodel.2$psill[1]/fitmodel.2$psill[2]
fitmodel.2
# model       psill range
# 1   Nug 0.000000000     0
# 2   Gau 0.004734831     5

weig<-var$np/var$dist^2  
SSTot<- sum((weig*(var$gamma-weighted.mean(var$gamma, weig)))^2)
SSErr <- attr(fitmodel.2,"SSErr");SSErr
R2 <- 1-(SSErr)/SSTot;R2

print(xyplot(gamma ~ dist, var, pch = 3, 
             type = 'b', lwd = 2, col = 'darkblue',
             panel = function(x, y, ...) {
               for (i in 1:100) {
                 trends$random = sample(trends$trend.pet.c)
                 v = variogram((random) ~ 1, trends,cutoff=10)
                 llines(v$dist, v$gamma, col = 'grey')
               }
               panel.xyplot(x, y, ...)
             },
             
             ylim = c(0, 0.004), xlab = 'distance', 
             ylab = 'semivariance'
))
plot(var, model=varmodel, add=TRUE) 


###### variable 3: precipitation #######

var=variogram(trends$trend.pm.c~1, data=trends,cutoff=10)
plot(var)
varmodel <- vgm(.13,"Gau",10,.06)
plot(var, model=varmodel) 
fitmodel.1 <- fit.variogram(var, model=varmodel, fit.sills=FALSE)
fitmodel.1$psill[1]/fitmodel.1$psill[2]
weig<-var$np/var$dist^2  
SSTot<- sum((weig*(var$gamma-weighted.mean(var$gamma, weig)))^2)
SSErr <- attr(fitmodel.1,"SSErr");SSErr
R2 <- 1-(SSErr)/SSTot;R2
fitmodel.1
# model psill    range
# 1   Nug  0.06 0.000000
# 2   Gau  0.13 7.859011

print(xyplot(gamma ~ dist, var, pch = 3, type = 'b', 
             lwd = 2, col = 'darkblue',
             panel = function(x, y, ...) {
               for (i in 1:100) {
                 trends$random = sample(trends$trend.pm.c)
                 v = variogram((random) ~ 1, trends,cutoff=10)
                 llines(v$dist, v$gamma, col = 'grey')
               }
               panel.xyplot(x, y, ...)
             },
             ylim = c(0, 0.17), xlab = 'distance', 
             ylab = 'semivariance'
))

###### variable 4: runoff #######

var=variogram(trends$trend.ro.c~1, data=trends,cutoff=10)
plot(var)
varmodel <- vgm(0.00010,"Gau",10,.00004)
plot(var, model=varmodel) 
fitmodel.5 <- fit.variogram(var, model=varmodel, fit.ranges=FALSE)
mean(trends$trend.ro.c)
fitmodel.5$psill[1]/fitmodel.5$psill[2]
fitmodel.5
# model        psill range
# 1   Nug 3.784385e-05     0
# 2   Gau 1.104516e-04    10

weig<-var$np/var$dist^2  
SSTot<- sum((weig*(var$gamma-weighted.mean(var$gamma, weig)))^2)
SSErr <- attr(fitmodel.5,"SSErr")
R2 <- 1-(SSErr)/SSTot;R2

print(xyplot(gamma ~ dist, var, pch = 3, type = 'b', 
             lwd = 2, col = 'darkblue',
             panel = function(x, y, ...) {
               for (i in 1:100) {
                 trends$random = sample(trends$trend.ro.c)
                 v = variogram((random) ~ 1, trends,cutoff=10)
                 llines(v$dist, v$gamma, col = 'grey')
               }
               panel.xyplot(x, y, ...)
             },
             ylim = c(0, 0.00012), xlab = 'distance', 
             ylab = 'semivariance'
))



###### variable 5: temperature #######

var=variogram(trends$trend.temp.c~1, data=trends)
plot(var)
varmodel <- vgm(psill=0.00009,model="Gau",7,nugget=0)
plot(var, model=varmodel) 
fitmodel.4 <- fit.variogram(var, model=varmodel, fit.ranges=FALSE)
mean(trends$trend.temp.c)
fitmodel.4$psill[1]/fitmodel.4$psill[2]
fitmodel.4
# model        psill range
# 1   Nug 5.042764e-08     0
# 2   Gau 8.051668e-05     7

weig<-var$np/var$dist^2  
SSTot<- sum((weig*(var$gamma-weighted.mean(var$gamma, weig)))^2)
SSErr <- attr(fitmodel.4,"SSErr")
R2 <- 1-(SSErr)/SSTot;R2

print(xyplot(gamma ~ dist, var, pch = 3, type = 'b', 
             lwd = 2, col = 'darkblue',
             panel = function(x, y, ...) {
               for (i in 1:100) {
                 trends$random = sample(trends$trend.temp.c)
                 v = variogram((random) ~ 1, trends,cutoff=10)
                 llines(v$dist, v$gamma, col = 'grey')
               }
               panel.xyplot(x, y, ...)
             },
             ylim = c(0, 0.000045), xlab = 'distance', 
             ylab = 'semivariance'
))

## t-test wrongly assuming independence (temporal and spatial)

set.seed(1818)
trend_WO <- as_tibble(trend_WO)

trends.sig <- trend_WO %>% 
  pivot_longer(-c(station,lat,lon),
               names_to = c("type","variable"),
               names_pattern = "(..).*\\.(.*?)\\.",
               values_to="value") %>% 
  pivot_wider(c(station,lat,lon,variable),
              names_from = type, 
              values_from = value) %>% 
  group_by(variable,station) %>% 
  mutate(trend.sig = ls < 0 | 0 < li)  %>% 
  mutate(sig_coef = trend.sig*tr) 

## visualisation

maps0 <- trends.sig %>% 
  split(.$variable) 

map.mk.ar<-mapping(map1,maps0$ari)
map.mk.pm<-mapping(map1,maps0$pm)
map.mk.pet<-mapping(map1,maps0$pet)
map.mk.ro<-mapping(map1,maps0$ro)
map.mk.tem<-mapping(map1,maps0$temp)

## Second, correct by the spatial dependency:

# Step 1: reconstruct Sigma 

Sigma.pm<-get.var.cov(fitmodel.1$psill[1],fitmodel.1$psill[2],
            fitmodel.1$range[2])
Sigma.pet<-get.var.cov(fitmodel.2$psill[1],fitmodel.2$psill[2],
                     fitmodel.2$range[2])
Sigma.ari<-get.var.cov(fitmodel.3$psill[1],fitmodel.3$psill[2],
                       fitmodel.3$range[2])
Sigma.tem<-get.var.cov(fitmodel.4$psill[1],fitmodel.4$psill[2],
                       fitmodel.4$range[2])
Sigma.ro<-get.var.cov(fitmodel.5$psill[1],fitmodel.5$psill[2],
                      fitmodel.5$range[2])

# Setp 2: construct the empirical distributions for each variable
# Step 3: map those trends that have conf intervals that do 
# not overlap with 0

set.seed(18)
trend.sig.pm<-get.sig(Sigma.pm,trend_WO$trend.pm.c)
trend.sig.pet<-get.sig(Sigma.pet,trend_WO$trend.pet.c)
trend.sig.ari<-get.sig(Sigma.ari,trend_WO$trend.ari.c)
trend.sig.tem<-get.sig(Sigma.tem,trend_WO$trend.temp.c)
trend.sig.ro<-get.sig(Sigma.tem,trend_WO$trend.ro.c)

trend.sig.pm.ind<-get.sig.ind(trend_WO$trend.pm.c)
trend.sig.pet.ind<-get.sig.ind(trend_WO$trend.pet.c)
trend.sig.ari.ind<-get.sig.ind(trend_WO$trend.ari.c)
trend.sig.tem.ind<-get.sig.ind(trend_WO$trend.temp.c)
trend.sig.ro.ind<-get.sig.ind(trend_WO$trend.ro.c)

signif<- tibble(
  sig=c(trend.sig.ari[[1]],trend.sig.pet[[1]],trend.sig.pm[[1]],
                 trend.sig.ro[[1]],trend.sig.tem[[1]]),
  sig.ind=c(trend.sig.ari.ind[[1]],trend.sig.pet.ind[[1]],
            trend.sig.pm.ind[[1]],
           trend.sig.ro.ind[[1]],trend.sig.tem.ind[[1]]),
  variable=c(rep(c("trend.aridity.c","trend.PET.c","trend.pm.c",
                   "trend.runoff.c","trend.Tavgmensual.c"),each=n0)))
dim(signif)

trends.sig2 <- trend_WO %>% 
  pivot_longer(-c(station,lat,lon),
               names_to = c("type","variable"),
               names_pattern = "(..).*\\.(.*?)\\.",
               values_to="value") %>% 
  pivot_wider(c(station,lat,lon,variable),
              names_from = type, 
              values_from = value) %>% 
  arrange(variable) %>% 
  bind_cols(signif) %>% 
  mutate(trend.sig=sig*tr,trend.sig.ind=sig.ind*tr)

## visualisation

maps <- trends.sig2 %>% 
  split(.$variable) 

map.cor.ar<-mapping.dep(map1,maps$ari,-1)
map.cor.ro<-mapping.dep(map1,maps$ro,1)
map.cor.pet<-mapping.dep(map1,maps$pet,1)
map.cor.tem<-mapping.dep(map1,maps$temp,1)
map.cor.pm<-mapping.dep(map1,maps$pm,-1)

map.ind.ar<-mapping.ind(map1,maps$ari)
map.ind.pm<-mapping.ind(map1,maps$pm)
map.ind.pet<-mapping.ind(map1,maps$pet)
map.ind.ro<-mapping.ind(map1,maps$ro)
map.ind.tem<-mapping.ind(map1,maps$temp)

data_CI <- maps0$ari
interval_ind <- trend.sig.ari.ind[[2]]
interval_dep <- trend.sig.ari[[2]]
plotCI(data_CI,interval_ind,interval_dep)

data_CI <- maps0$pm
interval_ind <- trend.sig.pm.ind[[2]]
interval_dep <- trend.sig.pm[[2]]
plotCI(data_CI,interval_ind,interval_dep)

data_CI <- maps0$pet
interval_ind <- trend.sig.pet.ind[[2]]
interval_dep <- trend.sig.pet[[2]]
plotCI(data_CI,interval_ind,interval_dep)

data_CI <- maps0$ro
interval_ind <- trend.sig.ro.ind[[2]]
interval_dep <- trend.sig.ro[[2]]
plotCI(data_CI,interval_ind,interval_dep)

data_CI <- maps0$temp
interval_ind <- trend.sig.tem.ind[[2]]
interval_dep <- trend.sig.tem[[2]]
plotCI(data_CI,interval_ind,interval_dep)


