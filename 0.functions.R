# Your custom code is a bunch of functions.

get.maps<-function(trend){
  
  data2<- as.data.frame(cbind(trend$station,
                              estaciones$lon,
                              estaciones$lat, 
                              trend$coef))
  names(data2)<- c("station","lon","lat","coef")
  d.max <- max(dist(data2[,c('lon',"lat")]))
  d.max #km
  
  ## Interpolate the coefficients:
  coordinates(data2) <- ~lon+lat
  proj4string(data2) <- CRS('+proj=longlat +datum=WGS84')
  TA <- CRS('+proj=longlat +datum=WGS84')
  data2<-remove.duplicates(data2)
  aq <- spTransform(data2, TA)
  ca <- spTransform(data2, TA)
  r <- raster::raster(ca)
  res(r) <- 0.1  # 10 km if your CRS's units are in km
  g <- as(r, 'SpatialGrid')
  gs <- gstat(formula=data2$coef~1, locations=aq)
  v <- variogram(gs, width=20)
  fve <- fit.variogram(v, vgm("Exp"))
  plot(variogramLine(fve, 400), type='l')
  points(v[,2:3], pch=20, col='red')
  k <- gstat(formula=data2$coef~1, locations=remove.duplicates(aq), model=fve)
  # predicted values
  kp <- predict(k,g)
  ## [using ordinary kriging]
  smooth.map<-spplot(kp)
  
  ok <- brick(kp)
  ok <- mask(ok, map1)
  names(ok) <- c('prediction','variance')
  smooth.map2<-spplot(ok)
  
  ## Recalculo los coeficientes signif:
  x <- matrix(rep(1:360,195),ncol=195)
  dis<-as.matrix(dist(coordinates(data2)))
  Sigma<-(fve$psill[2]-fve$psill[1])*(1-exp(-3*(dis/fve$range[2])^2))+fve$psill[1]
  var.beta <- Sigma 
  uni<-sqrt(unique(diag(Sigma)))
  
  trend2<-trend %>% 
    mutate(sig = coef>mean(coef)+1.96*uni|coef<mean(coef)-1.96*uni) %>% 
    mutate(sig_coef = sig*coef)
  
  map_trends <- trend2 %>% 
    ggplot(aes(lon, lat)) +
    geom_sf(data = map1, inherit.aes = FALSE) +
    coord_sf(ylim = c(7,18), xlim = c(-78, -93)) +
    geom_point(aes(fill = sig_coef, size = sig_coef), shape = 21) +
    scale_fill_divergent("sig_coef") +
    theme_bw()+
    scale_size_area(max_size = 6, guide = "none") +
    labs(title = "Centroamérica 1970 a 1999: Coeficientes temporales", 
         x = "Lon", y = "Lat") 
  
  return(list(smooth.map,smooth.map2, map_trends))
}

get.maps2<-function(trend){
  
  data2<- as.data.frame(cbind(trend$station,
                              estaciones$lon,
                              estaciones$lat, 
                              trend$coef))
  names(data2)<- c("station","lon","lat","coef")
  d.max <- max(dist(data2[,c('lon',"lat")]))
  d.max #km
  
  ## Interpolate the coefficients:
  coordinates(data2) <- ~lon+lat
  proj4string(data2) <- CRS('+proj=longlat +datum=WGS84')
  TA <- CRS('+proj=longlat +datum=WGS84')
  data2<-remove.duplicates(data2)
  aq <- spTransform(data2, TA)
  ca <- spTransform(data2, TA)
  r <- raster::raster(ca)
  res(r) <- 0.1  # 10 km if your CRS's units are in km
  g <- as(r, 'SpatialGrid')
  gs <- gstat(formula=data2$coef~1, locations=aq)
  v <- variogram(gs, width=20)
  fve <- fit.variogram(v, vgm("Sph"))
  plot(variogramLine(fve, 400), type='l')
  points(v[,2:3], pch=20, col='red')
  k <- gstat(formula=data2$coef~1, locations=remove.duplicates(aq), model=fve)
  # predicted values
  kp <- predict(k,g)
  ## [using ordinary kriging]
  smooth.map<-spplot(kp)
  
  ok <- brick(kp)
  ok <- mask(ok, map1)
  names(ok) <- c('prediction','variance')
  smooth.map2<-spplot(ok)
  
  ## Recalculo los coeficientes signif:
  x <- matrix(rep(1:360,195),ncol=195)
  dis<-as.matrix(dist(coordinates(data2)))
  #Sigma<-(fve$psill[2]-fve$psill[1])*(1-exp(-3*(dis/fve$range[2])^2))+
  #fve$psill[1]
  uni <- var.beta <- var(trend$coef)
  #uni<-sqrt(unique(diag(Sigma)))
  
  trend2<-trend %>% 
    mutate(sig = coef>mean(coef)+3.96*uni|coef<mean(coef)-3.96*uni) %>% 
    mutate(sig_coef = sig*coef)
  
  map_trends <- trend2 %>% 
    ggplot(aes(lon, lat)) +
    geom_sf(data = map1, inherit.aes = FALSE) +
    coord_sf(ylim = c(7,18), xlim = c(-78, -93)) +
    geom_point(aes(fill = sig_coef, size = sig_coef), shape = 21) +
    scale_fill_divergent("sig_coef") +
    theme_bw()+
    scale_size_area(max_size = 6, guide = "none") +
    labs(title = "Centroamérica 1970 a 1999: Coeficientes temporales", 
         x = "Lon", y = "Lat") 
  
  return(list(smooth.map,smooth.map2, map_trends))
}

get.maps_descriptives<-function(trends){
  
  data2<- as.data.frame(cbind(trends$station,
                              estaciones$lon,
                              estaciones$lat, 
                             # trends$trend.pm.c,
                              trends$trend.PET.c,
                              trends$trend.aridity.c,
                              trends$trend.Tavgmensual.c,
                              trends$trend.runoff.c))
  names(data2)<- c("station","lon","lat","coef.1","coef.2",
                   "coef.3","coef.4")#,"coef.5")
  d.max <- max(dist(data2[,c('lon',"lat")]))
  d.max #km
  
  map_trends <- data2 %>% 
    pivot_longer(names_to ="Variable",values_to ="Coef", 
                 -c(station,lon,lat)) %>% 
    ggplot(aes(lon, lat)) +
    theme_ipsum() +
    facet_wrap(.~Variable, nrow=2,
               labeller = labeller(Variable = 
                   c("coef.1" = "(A)",
                     "coef.2" = "(B)",
                     "coef.3" = "(C)",
                     "coef.4" = "(D)")
    )) +
    geom_sf(data = map1, inherit.aes = FALSE) +
    coord_sf(ylim = c(7,18), xlim = c(-78, -93)) +
    geom_point(aes(fill = Coef), size = 2, shape = 21) +
    scale_fill_viridis(name="Trend",discrete=FALSE) +
    theme(
      panel.spacing = unit(0.1, "lines"),
      strip.text.x = element_text(size = 12)
    ) +
    scale_size_area(max_size = 6, guide = "none") +
    labs(x = "Lon", y = "Lat") 
  
  return(map_trends)
}

get.map_descriptives<-function(trend){
  
  data2<- as.data.frame(cbind(trend$station,
                              estaciones$lon,
                              estaciones$lat, 
                              trend$coef))
  names(data2)<- c("station","lon","lat","coef")
  d.max <- max(dist(data2[,c('lon',"lat")]))
  d.max #km
  
  map_trends <- data2 %>% 
    ggplot(aes(lon, lat)) +
    theme_ipsum() +
    geom_sf(data = map1, inherit.aes = FALSE) +
    coord_sf(ylim = c(7,18), xlim = c(-78, -93)) +
    geom_point(aes(fill = coef), size = 3, shape = 21) +
    scale_fill_viridis(name="Trend",discrete=FALSE) +
    theme(
      panel.spacing = unit(0.1, "lines"),
      strip.text.x = element_text(size = 12)
    ) +
    scale_size_area(max_size = 6, guide = "none") +
    labs(x = "Lon", y = "Lat") 
  
  return(map_trends)
}


## Construct the var-cov matrix
get.var.cov<-function(c0,ce,ae){
  # Define distance matrix of zeros
  dist.mat <- matrix(0,n,n)
  # Fill in distances between each sampling point
  #dist.mat<-distances(data.frame(trends[,c("lat","lon")]))
  for (i in 1:n) {
    dist.mat[i,] <- sqrt( (trends$lat-trends$lat[i])^2 +
                            (trends$lon-trends$lon[i])^2 )
  }
  # Define elements of variance-covariance matrix for
  # exponential covariance function
  # Initialize to zeros
  exp.sigma <- matrix(0,n,n)
  # set parameters
  exp.sigma <- c0+ce*exp(-dist.mat^2/ae^2)
  return(exp.sigma)
}

## Null hypothesis assuming dependence

get.sig<-function(sigma.hat,y){
  set.seed(18)
null.hyp<-sapply(1:5000,function(i)
mvrnorm(n = 1, rep(0,n), sigma.hat))
qq<-apply(null.hyp,1,function(x)quantile(x,c(0.025,0.975)))
interval<-apply(qq,1,mean)
aa<-mapply(function(i)(findInterval(y[i],interval)!=1),1:n0)
return(list(aa,interval))
}

## Null hypothesis assuming dependence

get.sig.ind<-function(y){
  set.seed(18)
  null.hyp<-sapply(1:5000,function(i)
    rnorm(n = n, 0, sqrt(var(y))))
  qq<-apply(null.hyp,1,function(x)quantile(x,c(0.025,0.975)))
  interval<-apply(qq,1,mean)
  aa<-mapply(function(i)(findInterval(y[i],interval)!=1),1:n0)
  return(list(aa,interval))
}


mapping <- function(map1,data){
  map_coef<-data %>% 
    filter(mk<0.05) %>% 
    ggplot(aes(lon, lat)) +
    geom_sf(data = map1, inherit.aes = FALSE) +
    coord_sf(ylim = c(7,18), xlim = c(-78, -93)) +
    geom_point(aes(fill = sig_coef, 
                   size = sig_coef), 
               shape = 21) +
    scale_fill_gradient2(name="Trend",mid="white",
                       limits=range(data$tr)) +
    theme_ipsum() +
    theme(
      panel.spacing = unit(0.1, "lines"),
      strip.text.x = element_text(size = 8)) +
    scale_size_area(max_size = 6, guide = "none") +
    labs(x = "Lon", y = "Lat")
  return(map_coef)
}

mapping.dep <- function(map1,data,neg){
  map_coef<-data %>% ggplot(aes(lon, lat)) +
    geom_sf(data = map1, inherit.aes = FALSE) +
    coord_sf(ylim = c(7,18), xlim = c(-78, -93)) +
    geom_point(aes(fill = data$trend.sig, 
                   size = neg*data$trend.sig), 
                   shape = 21, na.rm=TRUE ) +
    scale_fill_gradient2(name="Trend",mid="white",
                         limits=range(data$tr)) +
    theme_ipsum() +
    theme(#plot.background=element_rect(fill="#f7f7f7"),
      #panel.background=element_rect(fill="#f7f7f7"),
      panel.spacing = unit(0.1, "lines"),
      strip.text.x = element_text(size = 8)) +
      scale_size_area(guide = "none") +
    labs(x = "Lon", y = "Lat") 
  return(map_coef)
}

mapping.ind <- function(map1,data){
  map_coef<-data %>% ggplot(aes(lon, lat)) +
    geom_sf(data = map1, inherit.aes = FALSE) +
    coord_sf(ylim = c(7,18), xlim = c(-78, -93)) +
    geom_point(aes(fill = data$trend.sig.ind, 
                   size = data$trend.sig.ind), 
                   shape = 21) +
    scale_fill_gradient2(name="Trend",mid="white",
                         limits=range(data$tr)) +
    theme_ipsum() +
    theme(
      panel.spacing = unit(0.1, "lines"),
      strip.text.x = element_text(size = 8)) +
    scale_size_area(max_size = 6, guide = "none") +
    labs(x = "Lon", y = "Lat") 
  return(map_coef)
}

get_dots <- function(data, units, title,tt,outliers){
  theme_set(theme_classic())
plotA <- ggplot(data, aes(x=station, y=coef)) + 
  geom_point(col="cyan4", size=1) +   # Draw points
  geom_segment(aes(x=station, 
                   xend=station, 
                   y=min(coef), 
                   yend=max(coef)), 
               linetype="dashed", 
               size=0.01) +   # Draw dashed lines
  geom_hline(yintercept = 0, linetype="dashed",size=0.5) +
  labs(title=tt, 
       #subtitle="Central American Stations", 
       caption=units,
       y=title) +  
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(axis.title.x=element_text()) +
  coord_flip()+ geom_point(data=outliers,
          aes(x=station, y=coef), color="lightpink3")
return(plotA)
}

plotCI <- function(data_CI,interval_ind,interval_dep){
  gg <- ggplot(data_CI, 
               aes(x=li,xend=ls, y=station, 
                   color=trend.sig,group=station)) + 
    geom_point(aes(x=tr, y=station), color="black")+
    geom_dumbbell( ) + 
    labs(x=NULL, 
         y=NULL) +
    geom_vline(xintercept = interval_dep, 
               color = "steelblue4", linetype="dotted")+
    geom_vline(xintercept = interval_ind, 
               color = "white", linetype="dotted")+
    geom_vline(xintercept = 0, 
               color = "black")+
    theme(plot.title = element_text(hjust=0.5, face="bold"),
          plot.background=element_rect(fill="#f7f7f7"),
          panel.background=element_rect(fill="#f7f7f7"),
          #panel.grid.minor=element_blank(),
          #panel.grid.major.y=element_blank(),
          #panel.grid.major.x=element_line(),
          axis.ticks=element_blank(),
          legend.position="none",
          panel.border=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  return(gg)
}

boot_CI <- function(data,formula){
  est <- function(formula, data, indices) {
    d <- data[indices,] # allows boot to select sample
    fit <- lm(formula, data=d)
    return(summary(fit)$coef[,"Estimate"])
  }
  # bootstrapping with 1000 replications
  results <- boot(data=data, statistic=est,
                  R=1000, formula=formula)
  aa<-boot.ci(results, type="bca")
  return(list(linf=aa$bca[c(4)],lsup=aa$bca[c(5)]))
}

mK_CI <- function(data,formula){
datos <- ts(data)
MannKendall(datos)
  
  return(list(linf=aa$bca[c(4)],lsup=aa$bca[c(5)]))
}
