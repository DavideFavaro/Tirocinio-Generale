
library(shiny)
library(shinyWidgets)
library(leaflet)
library(mapview)
library(curl)
#library(getSpatialData)
#library(sf)
library(sp)
#library(rgdal)
library(DT)
library(htmlwidgets)
library(rpostgis)

source("getSpatialData_dev.R")

Sys.setlocale("LC_MESSAGES", 'en_GB.UTF-8')
Sys.setenv(LANG = "en_US.UTF-8")
months.list<-1:12
names(months.list)<-format(ISOdate(2000, 1:12, 1), "%B")

addResourcePath('global', "../staticfiles/")

platforms.tables<-list("Sentinel-1"="s1", 
                       "Sentinel-2"="s2", 
                       "Sentinel-3"="s3" )
orbit.direction<-list("ASCENDING"="ASC", "DESCENDING"="DESC")
polarization.mode<-list("HH" ,  "HH HV"	,  "HV" , "VH"	, "VH VV", "VV"	, "VV VH")	 
product.type<-list("GRD",  "OCN",	 "RAW", "SLC"	 )
operationalmode<-list("EW", "IW",	 "SM", "WV")	 

platforms.vars<-list("Sentinel-1"=list("Mode"="operationalmode",
                                       "Type"="producttype",
                                       "Polar."="polarisationmode",
                                       "Orbit"="orbitdirection",
                                       "Date"="beginposition"), 
                       "Sentinel-2"="s2", 
                       "Sentinel-3"="s3" )

products<-NULL
dt.data.sp.df<-NULL
search.bounds<-NULL
download.script<-"wget --no-check-certificate  --trust-server-names  --content-disposition --continue --user=%s --password=%s   \"https://scihub.copernicus.eu/dhus/odata/v1/Products('%s')/\\$value\""
 
columns.lut<-list("Sentinel-1"=list(
  "lut"=rbind( 
      "identifier"=c("title", NA, F),
      "uuid"=c("uuid", NA, F),
      #  "url"=c("url", "tbIconCol", T),
      "size"= c("Size",NA,T),  
      "beginposition" =  c("Date", "as.Date",T) ,
      "sensoroperationalmode"= c("Mode", "as.factor",T),
      "producttype" = c("Product.Type", "as.factor",T),
      "polarisationmode" = c("P.Mode", "as.factor",T),
      "orbitdirection" =  c("Orbit", "as.factor",T)
  ), 
  "rownames"="uuid"
),
"Sentinel-2"=list(
  "lut"=rbind( 
    "identifier"=c("title", NA, F),
    "uuid"=c("uuid", NA, F),
    #  "url"=c("url", "tbIconCol", T),
    "size"= c("Size",NA,T),  
    "beginposition" =  c("Date", "as.Date",T) ,
    "sensoroperationalmode"= c("Mode", "as.factor",T),
    "producttype" = c("Product.Type", "as.factor",T),
    "orbitdirection" =  c("Orbit", "as.factor",T)
  ), 
  "rownames"="uuid"
)
)



tbIconCol <- function(p){
  sprintf("<a style=\"padding-right:5px;\" title=\"Download dataset!\" href=\"%s\" target=\"_blank\">%s</a>
                             <a style=\"padding-right:5px;\" title=\"Zoom to tile %s\" onclick=\"zoomLayer('%s');\">%s</a>", 
          p$url,
          as.character(icon("download", "fa-lg")),
          rownames(p),
          rownames(p),
          as.character(icon("search", "fa-lg")) 
  )
}


