module Global

#= libraries
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
=#

using DataFrames
using Dates


source("getSpatialData_dev.R")

Sys.setlocale("LC_MESSAGES", "en_GB.UTF-8")
Sys.setenv(LANG = "en_US.UTF-8")


dates = [ join([ Date(2000, x, 1), " ", Time(12) ]) for x in  1:12 ]
months_list = Dict(  dates .=> collect(1:12) )
addResourcePath("global", "../staticfiles/")

platforms_tables = Dict( "Sentinel-1" => "s1", "Sentinel-2" => "s2", "Sentinel-3" => "s3" )
orbit_direction = Dict( "ASCENDING" => "ASC", "DESCENDING" => "DESC" )
polarization_mode = [ "HH" ,  "HH HV"	,  "HV" , "VH"	, "VH VV", "VV"	, "VV VH" ]	 
product_type = [ "GRD",  "OCN",	 "RAW", "SLC"	]
operationalmode = [ "EW", "IW",	 "SM", "WV" ]	 

platforms_vars = Dict( "Sentinel-1" => Dict( "Mode" => "operationalmode",
                                             "Type" => "producttype",
                                             "Polar." => "polarisationmode",
                                             "Orbit" => "orbitdirection",
                                             "Date" => "beginposition" ), 
                       "Sentinel-2" => "s2", 
                       "Sentinel-3" => "s3" )

products = nothing
dt_data_sp_df = nothing
search_bounds = nothing
download_script = "wget --no-check-certificate  --trust-server-names  --content-disposition --continue --user=%s --password=%s   \"https://scihub.copernicus.eu/dhus/odata/v1/Products('%s')/\\$value\""

#+ AS ROWS
 #columns_lut = Dict(
 #  "Sentinel-1" => Dict(
 #    :lut => DataFrame([
 #      ( Row = "identifier", c1 = "title", c2 = missing, c3 = false ),
 #      ( Row = "uuid", c1 = "uuid", c2 = missing, c3 = false ),
 #      #  ( Row = "url", c1 = "url", c2 = "tbIconCol", c3 = true ),
 #      ( Row = "size", c1 = "Size", c2 = missing, c3 = true ),  
 #      ( Row = "beginposition", c1 = "Date", c2 = "as_Date", c3 = true ),
 #      ( Row = "sensoroperationalmode", c1 = "Mode", c2 = "as_factor", c3 = true ),
 #      ( Row = "producttype", c1 = "Product_Type", c2 = "as_factor", c3 = true ),
 #      ( Row = "polarisationmode", c1 = "P_Mode", c2 = "as_factor", c3 = true ),
 #      ( Row = "orbitdirection", c1 = "Orbit", c2 = "as_factor", c3 = true )
 #    ]),
 #    :rownames => "uuid"
 #  ),
 #  "Sentinel-2" => Dict(
 #    :lut => DataFrames([
 #      ( Row = "identifier", c1 = "title", c2 = missing, c3 = false ),
 #      ( Row = "uuid", c1 = "uuid", c2 = missing, c3 = false ),
 #      #  ( Row = "url", c1 = "url", "tbIconCol", c3 = true ),
 #      ( Row = "size", c1 = "Size", c2 = missing, c3 = true ),
 #      ( Row = "beginposition", c1 = "Date", c2 = "as_Date", c3 = true ) ,
 #      ( Row = "sensoroperationalmode", c1 = "Mode", c2 = "as_factor", c3 = true ),
 #      ( Row = "producttype", c1 = "Product.Type", c2 = "as_factor", c3 = true ),
 #      ( Row = "orbitdirection", c1 = "Orbit", c2 = "as_factor", c3 = true )
 #    ]), 
 #    :rownames => "uuid"
 #  )
 #)
#

#+ AS COLUMNS 
 columns_lut = Dict(
  "Sentinel-1" => Dict(
    :lut => DataFrame( 
        identifier = [ "title",  missing, false ],
        uuid = [ "uuid",  missing, false ],
        # url = [ "url", "tbIconCol", true ],
        size = [ "Size", missing,true ],  
        beginposition =  [ "Date", "as.Date",true ] ,
        sensoroperationalmode = [ "Mode", "as.factor",true ],
        producttype = [ "Product.Type", "as.factor",true ],
        polarisationmode = [ "P.Mode", "as.factor",true ],
        orbitdirection =  [ "Orbit", "as.factor",true ]
    ), 
    :rownames => "uuid"
  ),
  "Sentinel-2" => Dict(
    :lut => DataFrame( 
      identifier = [ "title", missing, false ],
      uuid = [ "uuid", missing, false ],
      #  url = [ "url", "tbIconCol", true ],
      size = [ "Size",missing,true ],  
      beginposition =  [ "Date", "as.Date",true ] ,
      sensoroperationalmode = [ "Mode", "as.factor",true ],
      producttype = [ "Product.Type", "as.factor",true ],
      orbitdirection =  [ "Orbit", "as.factor",true ]
    ), 
  :rownames => "uuid"
  )
)

function tbIconCol(p)
  println("<a style=\"padding-right:5px;\" title=\"Download dataset!\" href=\"$(p[:url])\" target=\"_blank\">$(string(icon("download", "fa-lg")))</a>
           <a style=\"padding-right:5px;\" title=\"Zoom to tile $(names(p))\" onclick=\"zoomLayer('$(names(p))');\">$(string(icon("search", "fa-lg")))</a>")
end

end # module