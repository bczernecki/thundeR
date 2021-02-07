#' Plotowanie diagramu Skew-T z danych z sondowaniami w formacie z GFSa
#' 
#' Przykłady na dole:
#' 
#' @param x  Ramka danych w formacie z sondowaniami z GFSa
#' @param wczyt_szablon Czy wyliczac szablon Skew-T(TRUE/FALSE); wymagane jeśli brak obiektu 'pp' w global env. ;
#' @param etykiety_hgt Czy dopisywac wysokosci glownych poziomow izobarycznych w metrach (TRUE/FALSE)
#' @param zaznacz_cape Czy rysowac poligon/punktu z miejscem CAPE'a/CINa (TRUE/FALSE)
#' @param wczytane_holdxy Czy wcześniej wczytano wyliczone macierzy holdx/holdy oraz poziomow 'p' dla pseudoadiabat przyspiesza wyznaczanie obszaru CAPE/CIN?
#' @param ramka Czy wyrysowac ramke danych z indeksami?
#'
#' @export
#' @examples
#' # sposób wywołanie funkcji:
#' # dla przykładowych danych wrzuconych do paczki:
#' x <- data.frame(lev=c(1000, 855, 700, 500, 300, 100, 10),
#' HGT=c(0, 1500, 2500, 6000, 8500, 12000, 25000),
#' TMP=c(25, 10, 0, -15, -30, -50, -92),
#' DPT = c(20, 5, -5, -30, -55, -80, -99),
#' WD = c(0, 90, 135, 180, 270, 350, 0),
#' WS = c(5, 10, 20, 30, 40, 5, 0),
#' lon=52,
#' lat=19,
#' date2=Sys.time())
#' 
#' #data("holdxy")
#' #data("holdxx")
#' #data("pp")
#' 
#' #plotuj(x, wczytane_holdxy=TRUE)
#' 
#' # przygotowanie szablonu:
#' #pp <- readRDS("data/pp.rds")
#' #system.time(plotuj(x, wczyt_szablon=T, wczytane_holdxy=T))
#' #system.time(plotuj(x, wczyt_szablon=F, wczytane_holdxy=T))


           
plotuj <- function(x, wczyt_szablon=TRUE, etykiety_hgt=TRUE, zaznacz_cape=TRUE, wczytane_holdxy=TRUE, ramka=TRUE){
        
        indeksy <- sounding(x$lev, x$HGT, x$TMP, x$DPT, x$WD, x$WS)
        
        # if(wczyt_szablon==TRUE){
        # data("pp")
        # data("holdxx")
        # data("holdyy")
        # data("p_hold")
        # }
        
        sounding_plot() # make a template
        #print(pp)
        sounding::skewt.lines(temp = x$TMP, pressure = x$lev, "red", grubosc = 2.5)
        sounding::skewt.lines(temp = x$DPT, pressure = x$lev, "blue", grubosc = 2.5)
        v <- skewty(x$lev)
        v[c(1,2)] <- c(0.2,1.1)
        u <- rep(28, times=length(v))
        
        # dodanie choragiewek:
        choragiewki(cx = u-2, cy=v, direction = x$WD, speed = x$WS, cex = 0.7)
        text(u+2, v, sprintf("%.1f",round(x$WS,1)), cex = 0.6)
        text(u[1]+2,45.4, "[m/s]", cex=0.6) # + jednostki
        
        # dodanie wysokosci w metrach
        if(etykiety_hgt==TRUE){
           pressure_linie<- c(1050, 1000, 850, 700, 500, 400, 300, 200, 150, 100)
           y <- skewty(pressure_linie)
           etykiety_hgt <- round(spline(x$lev,x$HGT, xout = pressure_linie)$y)
           text(-28.5,y[-c(1, length(y))]+0.4,paste(etykiety_hgt[-c(1, length(y))],"m"), cex=0.8, pos=4)
        }
        
        # dodanie poligonu kejpowego:
        if(zaznacz_cape==TRUE){
         
                          EL <- as.numeric(indeksy[122])
                          LFC <- as.numeric(indeksy[121])
                          LFC_temp <- (approx(x$lev,  x$TMP,xout=LFC)$y)
                          EL_temp <- (approx(x$lev,  x$TMP,xout=EL)$y)
                          y_lfc <- skewty(LFC)
                          x_lfc <- skewtx(temp = LFC_temp, ycoord = y_lfc)
                          points(x_lfc, y_lfc, pch=21, col="blue")
                          # i to samo dla EL:
                          y_el <- skewty(EL)
                          x_el <- skewtx(temp = EL_temp, ycoord = y_el)
                          points(x_el, y_el, pch=19, col="red")
                          
                          koordynaty <- c(x_lfc, y_lfc)  # koordynaty diagramu dla LCL'a
                          
                          if(indeksy[20]>30){
                          # wklejam ponownie kod do tworzenia macierzy na pseudoadiabaty:
                          #p <- round(seq(from = indeksy[120], to = indeksy[122], by = -5))
                          if( (LFC-EL)>0 & (EL<x$lev[1]+50)) {
                                  
                                  if(wczytane_holdxy==F){
                                          p_hold <- round(seq(from = 1050, to = 100, by=-5))
                                          npts <- length(p_hold)
                                          sy <- skewty(p_hold)
                                          sx <- double(length = npts)
                                          
                                          pseudo <- seq(from=32,to=0,by=-0.05)
                                          NPSEUDO <- length(pseudo)
                                          holdxx <- matrix(0, nrow = npts, ncol = NPSEUDO)
                                          holdyy <- matrix(0, nrow = npts, ncol = NPSEUDO)
                                          for (ipseudo in 1:NPSEUDO) {
                                                  for (ilen in 1:npts) {
                                                          moist <- satlft(pseudo[ipseudo], p_hold[ilen])
                                                          sx[ilen] <- skewtx(moist, sy[ilen])
                                                  }
                                                  inds <- (sx < -27.8)
                                                  sx[inds] <- NA
                                                  sy[inds] <- NA
                                                  holdxx[, ipseudo] <- sx
                                                  holdyy[, ipseudo] <- sy
                                          }
                                          
                                          for (ipseudo in 1:NPSEUDO) {
                                                  sx <- holdxx[, ipseudo]
                                                  sy <- holdyy[, ipseudo]
                                                  #                          lines(sx, sy, lty = 1, col = "red")
                                          }
                                        #save(list = c("p_hold", "holdxx", "holdyy"), file = "data/holdxy.Rdata")
                                  }
                                  
                                  # przycinanie macierzy do odpowiednich poziomow                
                                  holdx2 <- holdxx[(p_hold>=EL & p_hold<=LFC),]
                                  holdy2 <- holdyy[(p_hold>=EL & p_hold<=LFC),]
                                  
                                  if(!is.null(dim(holdx2))){ # na wypadek, gdyby tylko jeden wymiary pozostal po przycinaniu powyzej
                                          ind2 <- which.min(abs(holdx2[1,]-koordynaty[1])) # szukam rzedu z najblizszym x
                                          #opcjonalne wyrysowanie samej linii:
                                          #lines(  holdx2[,ind2],  holdy2[,1] ,lwd=4, col="gray") # w
                                          
                                          # czas na wyrysowanie polygonu zamiast linii
                                          dobre <- which(x$lev<=LFC & x$lev>=EL)
                                          
                                          iksy1 <- skewtx(c(LFC_temp, x$TMP[dobre], EL_temp),
                                                          skewty(c(LFC,x$lev[dobre], EL))) # iksy dla sondazu
                                          # brakuje jeszcze punktu poczatkowego dla LCLa i EL,
                                          # a przeciez one sa juz policzone!
                                          #iksy2 <- holdx[-nrow(holdx),ind2] # iksy dla pseudoadiabaty CAPEowej, ale najlepiej bez ostatniego elementu
                                          iksy2 <- holdx2[,ind2] # iksy dla pseudoadiabaty CAPEowej, ale najlepiej bez ostatniego elementu
                                          iksy <- c(iksy2,rev(iksy1))
                                          
                                          igreki1 <- skewty(c(LFC, x$lev[dobre], EL)) # igreki dla sondazu
                                          #igreki2=holdy[-nrow(holdy),1] # igreki dla pseudoadiabty capeowej, ale najlepiej bez ostatniej
                                          igreki2 <- holdy2[,1] # igreki dla pseudoadiabty capeowej, ale najlepiej bez ostatniej
                                          igreki <- c(igreki2,rev(igreki1))
                                          
                                          polygon(iksy,igreki,border = T, col = "#FF005050")
                                  } # koniec ifa sprawdzajacego holdx2
                                  
                                  points(x_el, y_el, pch=19, col="red") # jeszcze raz u góry kropa
                          }
                          
                          # 
                          # # dodanie statystyk
                          # rect(-27.5,30,-13.5,44, col = "#FFFFFF95",lwd = 2)
                          # text(-28,41,paste("ML CAPE", round(indeksy[36],1),sep="\t"), pos=4, col="red", cex=0.8)
                          # text(-28,39,paste("ML CIN", round(indeksy[38],1),sep="\t"), pos=4, col="red", cex=0.8)
                          # text(-28,37,paste("ML LFC", round(LFC,1),sep="\t"), pos=4, col="red", cex=0.8)
                          # text(-28,35,paste("ML EL", round(EL,1),sep="\t"), pos=4, col="red", cex=0.8)
        }
        
        } # koniec dla dodawania poligonu kejpowego
        
        
        # dodanie copylefta:
        # rect(-27.8,-1.9,-21,-0.2, lwd = 1, col="white")
        # text(-27.9,-1.2,"enwo.pl",cex=0.7, pos=4)
        
        # dodanie informacji o lon/lat + czasie
        rect(-27.8, 44.2,-18, 45.8, lwd = 1, col="white")
        text(-28.6,44.9,paste0(sprintf("%.2f", x$lon[1]), "E ", sprintf("%.2f",x$lat[1]),"N"),cex=0.7, pos=4)
        
        rect(-27.8, 42.2,-13, 43.8, lwd = 1, col="white")
        text(-28.5,42.9,paste0(format(x$date2[1], "%Y-%m-%d %H:%M"), " UTC"),cex=0.7, pos=4)
        
        
        
        # dodaje ramke z indeksami:
        if(ramka==T){
                par(pty = "m")
                par(fig=c(0.80, 0.99, 0.08, 0.94), new=T, mar=c(0, 0, 0, 0), oma=c(0,0,0,0))
                plot(0:10,0:10, xaxt='n', yaxt='n', col='white', frame.plot = F)
                rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "#FFFF5005")
                
                # wyrzucmy z automatu wszystkie wartosci wirtualne:
                ind2 <- (c(3, 11:19, 27:35, 43:53))
                indeksy <- indeksy[ind2]
                indeksy <- indeksy[order(names(indeksy))]
                
                text(rep(0, times=length(ind2)), y=seq(from=par("usr")[1]+0.2, to=par("usr")[2]-0.2, length.out = length(ind2)),
                     names(indeksy), cex=0.6, pos = 4)
                text(rep(7.5, times=length(ind2)), y=seq(from=par("usr")[1]+0.2, to=par("usr")[2]-0.2, length.out = length(ind2)),
                     round(indeksy,1), cex=0.6, pos = 4)
        }
        
        }