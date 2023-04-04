#' @keywords internal
#' @noRd

sounding_wyoming = function(wmo_id, 
                            yy, mm, dd, hh, min = 00, 
                            bufr = FALSE,
                            allow_failure = TRUE) {
  
  if (allow_failure) {
    tryCatch(sounding_wyoming_bp(wmo_id, 
                                 yy, mm, dd, hh, min, 
                                 bufr = bufr), 
             error = function(e) {
               message(paste("Problems with downloading data.",
                             "Run function with argument allow_failure = FALSE",
                             "to see more details"))})
  } else {
    sounding_wyoming_bp(wmo_id, 
                        yy, mm, dd, hh, min, 
                        bufr = bufr)
  }
}

#' @keywords internal
#' @noRd
sounding_wyoming_bp = function(wmo_id,
                               yy, mm, dd, hh, min, 
                               bufr = bufr) {
  
  if (length(yy) != 1 || length(mm) != 1 || length(dd) != 1 || length(hh) != 1) {
    stop("The function supports downloading data for a given day. Please change arguments yy, mm, dd, hh to single values")
  }
  
  if (length(wmo_id) != 1) {
    stop("The function supports downloading data for one station at the time. Please change the `wmo_id` argument to a single value")
  }
  
  mm = formatC(mm, width = 2, format = "d", flag = "0")
  dd = formatC(dd, width = 2, format = "d", flag = "0")
  hh = formatC(hh, width = 2, format = "d", flag = "0")
  min = formatC(min, width = 2, format = "d", flag = "0")
  
  if (bufr) {
    url = paste0("http://weather.uwyo.edu/cgi-bin/bufrraob.py?src=bufr&datetime=", 
                 yy, "-", mm, "-", dd, "+", hh, ":", min, ":00&id=", 
                 sprintf("%05d", wmo_id), "&type=TEXT:LIST")
  } else {
    url = paste0("http://weather.uwyo.edu/cgi-bin/sounding?TYPE=TEXT%3ALIST&YEAR=",
                 yy, "&MONTH=", mm, "&FROM=", dd, hh, "&TO=", dd, hh, "&STNM=",
                 sprintf("%05d", wmo_id))
  }
  
  temp = tempfile()
  test_url(url, temp)
  
  # run only if downloaded file is valid
  df = NULL
  if (!is.na(file.size(temp)) & (file.size(temp) > 800)) { 
    
    txt = read.fwf(file = temp, widths = 1000)
    sects = grep(pattern = "PRE>", x = txt$V1)
    if (length(sects) == 0) {
      stop("HTTP status was '503 Service Unavailable'. Have you provided a correct station id?
      Please check wmo_id numbers at:
      http://weather.uwyo.edu/upperair/sounding.html")
    }
    
    df = read.fwf(file = temp, skip = sects[1] + 4, widths = rep(7, 11),
                  n = (sects[2] - (sects[1] + 5)))
    
    colnames(df) = c("PRES", "HGHT", "TEMP", "DWPT", "RELH", "MIXR",
                     "DRCT", "SKNT", "THTA", "THTE", "THTV")
    
    if (bufr == FALSE) {
      # the section below is not valid for BUFR decoded data:
      txt = read.fwf(file = temp, skip = sects[2] + 1, widths = 1000,
                     n = (sects[3] - (sects[2] + 2)), stringsAsFactors = FALSE)$V1
      df2 = as.data.frame(matrix(data = unlist(strsplit(txt, split = ": ")), ncol = 2, byrow = TRUE))
      colnames(df2) = c("parameter"," value")
    } else {
      # for bufr data try to read only the most essential metadata
      ind = grep(pattern = "Observations", txt$V1)
      df2 = data.frame(bufr_metadata = gsub("<.*?>", "", txt$V1[ind:(ind + 1)]), 
                       stringsAsFactors = FALSE)
      # and convert m/s to knots to stay in alignment with the default format used:
      df$SKNT = round(df$SKNT * 1.9438, 1)
    }
    df = list(df, df2)
    
  } else { # end of checking file size / problems with internet connection
    message(paste0("Service not working or wmo_id or date not correct. Check url:\n", url)) 
  }
  
  unlink(temp)
  return(df)
}
