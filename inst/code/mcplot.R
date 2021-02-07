#' Plotuje diagramy Skew-T z listy w sposób zrównoleglony
#'
#' Funkcja automatyzująca plotowanie wczytanej listy z GFS do diagramów Skew-T do katalogu 
#' z wykorzystaniem mclapply.
#' Jest to funkcja wewnętrza jedynie do zrównoleglenia 'plotuj()' za pomocą mclapply'a.
#' Przy wywoływaniu należy pamiętać o podaniu odpowiedniej liczy rdzeni (jak w przykładzie) oraz poprawnej ścieżki.
#' 
#' @param path katalog w ktorym beda zapisywane pliki; domyślnie 'skewt'
#' @param format format graficzny pliku wyjsciowego; svg lub png
#' @param wczytane_holdxy TRUE/FALSE: czy load("data/holdxy.Rdata") zostalo wczesniej wywolane w celu przyspieszenia plotowania
#' @param wczyt_szablon Czy wczytywac pp <- readRDS("data/pp.rds"); Jesli zostalo wczesniej wywolane w celu przyspieszenia plotowania ustawic na 'FALSE' 
#' 
#' @import parallel
#' @export
#' @examples
#' # przykładowe wywołanie funkcji z uruchomienie mcplota w pliku mcplot.R na dole
#' #library(parallel)
#' 


mcplot <- function(x, path="skewt", format="svg", wczytane_holdxy=TRUE){
        dir.create(path, showWarnings = FALSE)
        plik <- paste0(path,"/", format(x$date2[1], "%Y%m%d_%Hh-"),  
                       sprintf("%.2f", x$lon[1]), "E", sprintf("%.2f",x$lat[1]),"N.",format)
        if(format=="svg") svg(plik, width=9)
        if(format=="png") png(plik)
        plotuj(x, wczyt_szablon=FALSE, wczytane_holdxy = wczytane_holdxy)
        dev.off()
}



#dane <- read_gfs(plik_gfs="inst/data/sondaz2.gz")
#dane <- readRDS("inst/data/gfs_lista.rds")
#source("R/sounding_plot.R")
#source("R/station_symbol.R") # choragiewki
#data("pp") # pozwala na unikniecie wczytywania kodu zrodlowego; mocno przyspiesza:
#x <- dane[[3000]]
#system.time(plotuj(x, wczyt_szablon=T))
#system.time(plotuj(x, wczyt_szablon=F))

#system.time(lapply(dane[1:10], mcplot, wczytane_holdxy=FALSE))

# przygotwanie do nieco wiekszej proby - najlepiej jesli szablony zostana wczytane wczesniej w celu optymalizacji:
#data("holdxx")
#data("holdyy")
#data("pp")
#data("p_hold")
#system.time(mclapply(dane[1:100], mcplot, format="svg", wczytane_holdxy=TRUE, mc.cores = 4))

