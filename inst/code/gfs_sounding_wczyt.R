library(readr)
library(dplyr)
library(tidyr)
library(sounding)
library(parallel)


df <- read_delim("~/gfs/sondaz.gz", delim=",", col_names=FALSE, col_types='TTccddd') 
colnames(df) <- (c("date1", "date2", "var", "lev", "lon", "lat", "value")) 
df <- df[-which(df$lev=="surface" & df$var=="TMP"), ] # pozbywamy sie niepotrzebnej temperatury surface

# zamieniamy 10m i surface na 2m.
df[which(df$lev=="10 m above ground"),"lev"] <- "2 m above ground"
df[which(df$lev=="surface"),"lev"] <- "2 m above ground"
df <- filter(df, var!="APTMP", #var!="PRMSL", 
             lev!="5 mb", lev!="7 mb", lev!="10 mb", lev!="20 mb", lev!="30 mb",  lev!="50 mb",
             lev!="PV=2e-06 (Km^2/kg/s) surface", lev!="mean sea level")

df <- df %>% #filter(date2>=Sys.time()) %>%  # przefiltrowanie po dacie tylko dla potrzeb testow lub operacyjnego dzialania
        arrange(var, date2) %>%
        spread(data=., key=var, value=value)
df$TMP <- df$TMP-273.15
df$DPT <- df$DPT-273.15

# konwersja wiatru:
ws <- sqrt((df$UGRD^2) + (df$VGRD^2))
df$WD <- round((180/pi * atan2(df$UGRD, df$VGRD) + 180))
df$WS <- ws
rm(ws) # usuniecie smieci

# konwersja RH na dewpoint
b <- 17.67   # Bolton
c <- 243.5   # Bolton
# now calculate the dew-point temperature
fact <- log(df$RH/100)+b*(df$TMP)/(c+df$TMP)
df$DPT <- round(c*fact/(b-fact),2)
df <- dplyr::select(df, -UGRD, -VGRD, -RH, -date1)
rm(list=c("b","c","fact"))

# oznaczmy schematycznie najnizszy poziom jako 2000 mb
df[which(df$lev=="2 m above ground"),"lev"] <- "2000 mb" #
# i przekonwertujmy na tryb numeryczny
df$lev <- as.numeric(gsub("[^0-9\\.]", "", df$lev) )

#dla kolumny 'PRES' od razu doprowadzmy do finalnej postaci
df$PRES <- round(df$PRES/100,1)

# konwersja do listy:
library(data.table)
df <- as.data.table(df)
#test <- split(df, list(df$lon, df$lat, df$date2))

library(microbenchmark)
test <- split(df, by=c("lon", "lat", "date2"))

# tutaj rozpoczac czyszczenie:
#library(microbenchmark)
# b <- microbenchmark(czyszczenie(a), times=50)
#lista <- microbenchmark(mclapply(test2, FUN=czyszczenie, mc.cores = 2), times = 5)
#lista


#czyszcenie 
#a <- df %>% filter(lon==14, lat==49, date2==min(date2))

        czyszczenie <- function(a){
                ind_2m <- which(a$lev==2000)
                a$lev[ind_2m] <- a$PRES[ind_2m]
                alt <- a$HGT[ind_2m]
                a %>% dplyr::filter(HGT>=alt) %>% dplyr::arrange(HGT) %>% 
                        dplyr::select(.,-PRES) 
        }

#czyszczenie(a)

# 
 p3 <- Sys.time()
 lista <- mclapply(test, FUN=czyszczenie, mc.cores = 10)
 p3 - Sys.time()
 
 as.numeric(difftime(p3, Sys.time(), units="secs"))
# 
# p2 <- Sys.time()
# lista <- mclapply(test, FUN=czyszczenie, mc.cores = 2)
# p2 - Sys.time()
# 
# p1 <- Sys.time()
# lista <- lapply(test, FUN=czyszczenie)
# p1 - Sys.time()

        
# library(lineprof)
# tmp <- tempfile()
# Rprof(tmp, interval = 0.0001)
# lista <- microbenchmark(lapply(test2, FUN=czyszczenie), times = 5)
# Rprof(NULL)
# summaryRprof(filename = tmp)



#sounding(pressure=a$lev, a$HGT, a$TMP, a$DPT, a$WD, a$WS)
#rysuj_skewt(a$lev, a$HGT, a$TMP, a$DPT, a$WD, a$WS)
# dodanie copylefta:
rect(-27.8,-1.9,-21,-0.2, lwd = 1, col="white")
text(-28,-1.2,"enwo.pl",cex=0.7, pos=4)

# 
# library('tidyverse')
# library('parallel')
# library('methods')
# coord_df <- expand.grid(lat=seq(from=49, to=55, by=0.25),
#                         lon=seq(from=14, to=25, by=0.25))
# mc_plot <- function(id, gfs, coord_df){
#         print_meteorogram(coord_df[[1]][[id]], coord_df[[2]][[id]], gfs=gfs)
# }


source("R/sounding_plot.R")
dev.off()
system.time(print(p1))
source("R/station_symbol.R") # choragiewki


plotuj <- function(x){
        print(p1)
        skewt.lines(temp = x$TMP, pressure = x$lev, "brown", grubosc = 2)
        skewt.lines(temp = x$DPT, pressure = x$lev, "blue", grubosc = 2)
        v <- skewty(x$lev)
        v[c(1,2)] <- c(0.2,1.1)
        u <- rep(28, times=length(v))

        # dodanie choragiewek:
        station.symbol(cx = u-1, cy=v, direction = x$WD, speed = x$WS, cex = 0.7,circle = F)
        text(u+3, v, sprintf("%.1f",round(x$WS,1)), cex = 0.6)
        text(u[1]+3,45.4, "[m/s]", cex=0.6) # + jednostki

        # dodanie copylefta:
        rect(-27.8,-1.9,-21,-0.2, lwd = 1, col="white")
        text(-27.9,-1.2,"enwo.pl",cex=0.7, pos=4)

        # dodanie informacji o lon/lat + czasie
        rect(-27.8, 44.2,-18, 45.8, lwd = 1, col="white")
        text(-28.6,44.9,paste0(sprintf("%.2f", x$lon[1]), "E ", sprintf("%.2f",x$lat[1]),"N"),cex=0.7, pos=4)

        rect(-27.8, 42.2,-13, 43.8, lwd = 1, col="white")
        text(-28.5,42.9,paste0(format(x$date2[1], "%Y-%m-%d %H:%M"), " UTC"),cex=0.7, pos=4)
}


mcplot <- function(listewka){
        
        plik <- paste0(format(listewka$date2[1], "%Y%m%d_%Hh-"),  
                       sprintf("%.2f", listewka$lon[1]), "E", sprintf("%.2f",listewka$lat[1]),"N",".svg")
        print(plik)
        svg(plik)
        plotuj(listewka)
        dev.off()
}

system.time(mclapply(lista, mcplot, mc.cores = 4))
