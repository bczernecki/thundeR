#!/bin/bash

# Run R CMD check with valgrind to collect and detect reported errors.
#
# It is used by the valgrind GH actions workflow and can also be run locally,
# provided Docker is available.

mkdir -p valgrind-check

# Note that we use -l and --install-args="--clean" in order to not have in the
# valgrind-check directory (uploaded as artifact in GH actions) some large
# directories and files produced as part of R CMD check
docker run --rm -v $(pwd):/thunder bczernecki/thunder3 bash -c ' \
  R -e "install.packages(\"remotes\")" \
  R -e "remotes::install_github(\"bczernecki/climate\", dependencies = TRUE)" \
  R -e "remotes::install_deps(\"thunder\", dependencies = TRUE)" \
  R CMD build /thunder \
  R CMD INSTALL thunder*.tar.gz \
  R -d "valgrind --dsymutil=yes --leak-check=full --track-origins=yes" -e "pressure = c(1000, 855, 700, 500, 300, 100, 10);altitude = c(0, 1500, 2500, 6000, 8500, 12000, 25000);    temp = c(25, 10, 0, -15, -30, -50, -92) ;    dpt = c(20, 5, -5, -30, -55, -80, -99) ;    wd = c(0, 90, 135, 180, 270, 350, 0) ;    ws = c(5, 10, 20, 30, 40, 5, 0) ;    accuracy = 2 ;    thunder:::sounding_default(pressure, altitude, temp, dpt, wd, ws, 0, accuracy, interpolate_step = 1000, meanlayer_bottom_top = c(0,500), storm_motion = c(999, 999))" \
  '
