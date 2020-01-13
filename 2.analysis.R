## Data and packages

source(file="0.packages.R")
source(file="0.functions.R")
load(file="data/trends.Rdata")
map1 <- rnaturalearth::ne_states(
  country = c("guatemala", "honduras", 
              "el salvador", "panama", 
              "nicaragua", "costa rica", 
              "belize", "mexico", 
              "colombia"), returnclass = "sf")

## Set coordinates and projections:

coordinates(trends) <- ~lon+lat

## Estimate variograms for each variable

###### variable 1: pm #######

var=variogram(trends$trend.pm.c~1, data=trends,cutoff=20)
plot(var)
varmodel <- vgm(.8,"Gau",20,.1)
plot(var, model=varmodel) 
fitmodel.1 <- fit.variogram(var, model=varmodel, fit.sills=FALSE)
# model psill    range
# 1   Nug  0.10   0.0000
# 2   Gau  0.80   15.25

###### variable 2: PET #######

var=variogram(trends$trend.PET.c*10~1, data=trends)
plot(var)
varmodel <- vgm(psill=0.5,model="Gau",5,nugget=0)
plot(var, model=varmodel) 
fitmodel.2 <- fit.variogram(var, model=varmodel, fit.ranges=FALSE)
# model        psill range
# 1   Nug 0.0005155825     0
# 2   Gau 0.4763303455     5

###### variable 3: aridity #######

var=variogram(trends$trend.aridity.c*100~1, data=trends,cutoff=20)
plot(var)
varmodel <- vgm(.15,"Gau",20,.010)
plot(var, model=varmodel) 
fitmodel.3 <- fit.variogram(var, model=varmodel, fit.sills=FALSE)
# model psill    range
# 1   Nug 0.010   0.0000
  # 2   Gau 0.15  17.09

###### variable 4: tavg #######

var=variogram(trends$trend.Tavgmensual.c*100~1, data=trends)
plot(var)
varmodel <- vgm(psill=0.4,model="Gau",4.5,nugget=0)
plot(var, model=varmodel) 
fitmodel.4 <- fit.variogram(var, model=varmodel, fit.ranges=FALSE)
# model     psill range
# 1   Nug 0.0000000   0.0
# 2   Gau 0.3957163   4.5

###### variable 5: runoff #######

var=variogram(trends$trend.runoff.c*100~1, data=trends,cutoff=30)
plot(var)
varmodel <- vgm(2.5,"Gau",10,.8)
plot(var, model=varmodel) 
fitmodel.5 <- fit.variogram(var, model=varmodel, fit.ranges=FALSE)
# model     psill range
# 1   Nug 0.8231868     0
# 2   Gau 2.7433864    10
  
## t-test vs assuming dependence

## First, assume independence and have t-test:

x <- matrix(rep(1:360,195),ncol=195)
dis<-as.matrix(dist(coordinates(trends)))
trends <- as.tibble(trends)

trends.sig <- trends %>% 
  pivot_longer(-c(station,lat,lon),
               names_to = "variable",values_to = "trend") %>% 
  group_by(variable) %>% 
  mutate(var.trend = sqrt(var(trend))/length(trend)) %>% 
  mutate(trend.sig = trend>mean(trend)+1.96*var.trend | 
                     trend<mean(trend)-1.96*var.trend) %>% 
  mutate(sig_coef = trend.sig*trend) 

## visualisation

maps <- trends.sig %>% 
  split(.$variable) 

map.t.ar<-mapping(map1,maps$trend.aridity.c)
map.t.pm<-mapping(map1,maps$trend.pm.c)
map.t.pet<-mapping(map1,maps$trend.PET.c)
map.t.ro<-mapping(map1,maps$trend.runoff.c)
map.t.tem<-mapping(map1,maps$trend.Tavgmensual.c)

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

trend.sig.pm<-get.sig(Sigma.pm,trends$trend.pm.c)
trend.sig.pm.ind<-get.sig.ind(trends$trend.pm.c)
trend.sig.pet<-get.sig(Sigma.pet,trends$trend.PET.c*10)
trend.sig.pet.ind<-get.sig.ind(trends$trend.PET.c*10)
trend.sig.ari<-get.sig(Sigma.ari,trends$trend.aridity.c*100)
trend.sig.ari.ind<-get.sig.ind(trends$trend.aridity.c*100)
trend.sig.tem<-get.sig(Sigma.tem,trends$trend.Tavgmensual.c*100)
trend.sig.tem.ind<-get.sig.ind(trends$trend.Tavgmensual.c*100)
trend.sig.ro<-get.sig(Sigma.tem,trends$trend.runoff.c*100)
trend.sig.ro.ind<-get.sig.ind(trends$trend.runoff.c*100)

signif<- tibble(
  sig=c(trend.sig.ari,trend.sig.pet,trend.sig.pm,
                 trend.sig.ro,trend.sig.tem),
  sig.ind=c(trend.sig.ari.ind,trend.sig.pet.ind,trend.sig.pm.ind,
           trend.sig.ro.ind,trend.sig.tem.ind),
  variable=c(rep(c("trend.aridity.c","trend.PET.c","trend.pm.c",
                   "trend.runoff.c","trend.Tavgmensual.c"),each=199)))

trends.sig <- trends %>% 
  pivot_longer(-c(station,lat,lon),
               names_to = "variable",values_to = "trend") %>% 
  arrange(variable) %>% 
  bind_cols(signif) %>% 
  mutate(trend.sig=sig*trend,trend.sig.ind=sig.ind*trend)

## visualisation

maps <- trends.sig %>% 
  split(.$variable) 

map.cor.ar<-mapping.dep(map1,maps$trend.aridity.c)
map.cor.pm<-mapping.dep(map1,maps$trend.pm.c)
map.cor.pet<-mapping.dep(map1,maps$trend.PET.c)
map.cor.ro<-mapping.dep(map1,maps$trend.runoff.c)
map.cor.tem<-mapping.dep(map1,maps$trend.Tavgmensual.c)

map.ind.ar<-mapping.ind(map1,maps$trend.aridity.c)
map.ind.pm<-mapping.ind(map1,maps$trend.pm.c)
map.ind.pet<-mapping.ind(map1,maps$trend.PET.c)
map.ind.ro<-mapping.ind(map1,maps$trend.runoff.c)
map.ind.tem<-mapping.ind(map1,maps$trend.Tavgmensual.c)



