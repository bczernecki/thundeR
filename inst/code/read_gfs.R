#' Wczyt danych z sondowaniami z GFSa
#' 
#' Funkcja jest jedynie zrównolegleniem funkcji 'plotuj()' za pomocą mclapply'a.
#' Przy wywoływaniu należy pamiętać o podaniu odpowiedniej liczy rdzeni (jak podano w przykładzie).
#' 
#' @param plik_gfs  Ścieżka do pliku .gz z sondowaniami z GFSa
#' @param cores Liczba rdzeni do przygotowania podlisty z uzyciem mclapply
#'
#' @import dplyr
#' @import tidyr
#' @import parallel
#' @importFrom readr read_delim
#'
#' @export
#' @examples
#' # sposób wywołanie funkcji:
#' # dla przykładowych danych wrzuconych do paczki:
#' # dane <- read_gfs(plik_gfs="inst/data/sondaz2.gz")
#' 


read_gfs <- function(plik_gfs, cores=10){
        # library(readr)
        # library(dplyr)
        # library(tidyr)
        # library(parallel)

        df <- read_delim(file = plik_gfs, delim=",", col_names=FALSE, col_types='TTccddd') 
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
        test <- split(df, by=c("lon", "lat", "date2"))
        
        czyszczenie <- function(a){
                ind_2m <- which(a$lev==2000)
                a$lev[ind_2m] <- a$PRES[ind_2m]
                alt <- a$HGT[ind_2m]
                a %>% dplyr::filter(HGT>=alt) %>% dplyr::arrange(HGT) %>% 
                        dplyr::select(.,-PRES) 
        }
        
        lista <- mclapply(test, FUN=czyszczenie, mc.cores = cores)
        return(lista)
} 
