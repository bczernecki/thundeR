library(parallel)
dane <- read_gfs(plik_gfs="~/gfs/sondaz.gz", cores=12)
data("pp")
data("holdxx")
data("holdyy")
data("p_hold")

x <- dane[[3000]]
system.time(plotuj(x, wczyt_szablon=F))

system.time(lapply(dane[1:2000], mcplot, wczytane_holdxy=T))

# przygotwanie do nieco wiekszej proby - najlepiej jesli szablony zostana wczytane wczesniej w celu optymalizacji:
data("holdxx")
data("holdyy")
data("p_hold")
system.time(mclapply(dane, mcplot, format="svg", wczytane_holdxy=TRUE, mc.cores = 12))