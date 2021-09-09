module Esacron


"""imports
  source("getSpatialData_dev.R")
  library(sf)
  library(sp)
  library(rgeos)
  library(rpostgis) 
  library(rgdal)
"""

using ArchGDAL, CategoricalArrays, CSV, DataFrames, Dates, LibPQ

include("getSpatialDataDev.jl")
include("ValuableAdditions.jl")

#+ postgresql://"user":"Password"@"host":"port"/"database"
conn = LibPQ.Connection("postgresql://****:*********@localhost:27017/esa")
conn2 = LibPQ.Connection("postgresql://****:*********@geolab.vs-ix.net:5432/esa")

startdate = Dict( "Sentinel-1" => Date("2024-04-01"),
                  "Sentinel-2" => Date("2015-06-01")#=,
                  "Sentinel-3" => Date("2016-02-01") =#)
tables = Dict( "Sentinel-1" => "s1",
               "Sentinel-2" => "s2",
               "Sentinel-3" => "s3" )

platforms_lastdate = DataFrame(CVS.File("platforms.lastdate.csv"))

#platforms.lastdate$`Sentinel-1`<-platforms.lastdate$`Sentinel-1`-2
#platforms.lastdate$`Sentinel-2`<-platforms.lastdate$`Sentinel-2`-2

## 1 day  timespan
timespan = 1
platform = "Sentinel-2"  
login_CopHub( username = "fpirotti", password = "libero" ) #asks for password or define 'password'
set_archive("/tmp")

reloaday = false

for platform in names(platforms_startdate)
  if reloaday
    lns = readlines( platform*".log" )
    rename!(lns, replace.( lns, "problem" => "" ) ) 
    
    ww = findall(occursin( "problem", lns ))
    
    for date in names(lns[ww])
      println(date)
      date = Date(date)
      time_range =  [ string(date-timespan), string(date) ]
      products = getSentinel_query2( time_range = time_range,  platform = platform )
      products_tmp = products 
      
      if isnothing(products) || nrow(products) < 1
        lns[string(date)] = string(date)*" -- 0 righe "  
        
        lns_out = write( platform*".log", lns )
      end
      
      products = products[ findall( !ismissing, products[:footprint] ) , : ]

      insert!( products, 1, :index => products[:uuid] ) 

      products[:format] = CategoricalArray(products[:format])
      products[:format] = CategoricalArray(products[:format])
      replace!.( products[:orbitdirection], Ref("ASCENDING" => "ASC"), Ref("DESCENDING" => "DESC") )
      
      products[:orbitnumber] = convert.( Int64, products[:orbitnumber] )
      if :lastorbitnumber in names(products)
        products[:lastorbitnumber] = convert( Int64, products[:lastorbitnumber] )
      end
      if :lastrelativeorbitnumber in names(products)
        products[:lastrelativeorbitnumber] = convert( Int64, products[:lastrelativeorbitnumber] )
      end
      products[:polygon] = nothing
      products[:gmlfootprint] = nothing
      products[:title] = nothing
      products[:filename] = nothing 
      products[:summary] = nothing
      
      products[:platformname] = nothing
      #products[:platformname] = convert( Int64, products[:platformname] )
      products[:instrumentname] = nothing
      #products[:instrumentname] = convert( Int64, products[:instrumentname] )
      products[:instrumentshortname] = platform[10:12]
      products[:instrumentshortname] = convert( Int64, products[:instrumentshortname] )
      
      
      if "status" in names(products) 
        products[:status] = CategoricalArray(products[:status])
        products[:status] = as.numeric(products[:status]) #ARCHIVED
      end
      ## icon just add to url.alt + Products('Quicklook')/$value
      ## url just add to url.alt + $value
      products[:url] = nothing
      products[:url_icon] = nothing
      products[:url_alt] = nothing

      products[:month] = format( Date(products[:beginposition]), "%m") 
      tmp_size = split( products[:size], " " )
      products[:size] = map( x -> { 
        ret = missing
        if x[2] == "MB"
          ret = as.numeric(x[1])
        end
        if x[2] == "GB"
          ret = as.numeric(x[1]) * 1000
        end
        if x[2] == "KB"
          ret = as.numeric(x[1]) / 1000
        end
        return ret
      }, tmp_size )
      dt_data = map( x -> {
        tg = ArchGDAL.fromWKT(products[x,"footprint"])
        pl = tg.polygons[1]
        pl.ID = string.(products[x,"uuid"])
        return pl
      }, products[:,1] )
      
      #ftp = products[:footprint]
      products[:footprint] = nothing
      
      dt_data_sp = SpatialPolygons( dt_data, proj4string = CRS("+proj=longlat +datum=WGS84") )
      dt_data_sp_df = SpatialPolygonsDataFrame( dt_data_sp, products )
      nPolys = map( x -> length(x.Polygons), dt_data_sp.polygons )
      dd = dt_data_sp_df[ findall( nPolys > 1 ), : ] 
      outins = pgInsert( conn, partial_match = true, name = ( "eo_sathub", platforms_tables[platform] ), 
                         data_obj = dt_data_sp_df, upsert_using = "uuid,identifier" )
      outins = pgInsert( conn2, partial_match = true, name = ( "eo_sathub", platforms_tables[platform] ), 
                         data_obj = dt_data_sp_df, upsert_using = "uuid,identifier" )
      
      if outins
        lns[string(date)] = string(date)
      end
      println(date)
    
      
      
      log_con  =  open( platform*".log", "w" )
      lns_out = write( log_con, join(lns) ) 
      close(log_con) 
    end
  end









         
  ##remove multi polygon!!
 # sql = "SELECT   beginposition  FROM eo_sathub.s2 where ST_NumGeometries(geom)>1 group by beginposition order by beginposition"
  
  #redo = dbGetQuery(conn, sql)
  #for(date in redo[:beginposition])
  #Dates.today() - 1 per evitare ingestion / beginposition differenza
  while platforms_lastdate[platform] < Dates.today() 
    
    println(platform)
    
    println(platforms_lastdate[platform])
    
    date = platforms_lastdate[platform]
    
    time_range = [ string(date), string( date + timespan ) ]
    
    # products = getSentinel_query2( time_range = time_range,  platform = platform )

    #, extra = "processinglevel:Level-2Ap"
    
    products = getSentinel_query2( time_range = time_range,  platform = platform,
                                   aoi = nothing, username = nothing, 
                                   password = nothing, hub = "auto", verbose = true, ingestiondate = true )
    
    
    #CSV.write( "products.jld2", products ) 
    products_tmp = products
    platforms_lastdate[platform] += timespan 
    
    if isnothing(products) || nrow(products) < 1

      CSV.write( "platforms_lastdate.csv", platforms_lastdate )
      
     # log_con = open( platform*".log", "r" )
      open( platform*".log", "w" ) do io
        write( io, string( platforms_lastdate[platform] )[1]*" problem0" )
      end
      println( platforms_lastdate[platform] )
    end
    
    products = products[ findall( !ismissing, products[:footprint] ), : ]
    unique!(products)
    insert!( products, 1, :uuid => products[:uuid] )
    
    #uuids = dbGetQuery( conn, "select uuid from eo_sathub.s2 where processinglevel='Level-2Ap' " )
    #tokeep = which( !(rownames(products) %in% uuids$uuid)  )
    
    products[:format] = CategoricalArray(products[:format])
    replace!.( products[:orbitdirection], Ref("ASCENDING" => "ASC"), Ref("DESCENDING" => "DESC") )
    
    products[:orbitnumber] = convert( Int64, products[:orbitnumber] )
    if "lastorbitnumber" in names(products)
      products[:lastorbitnumber] = convert( Int64, products[:lastorbitnumber] )
    end
    if "lastrelativeorbitnumber" in names(products)
      products[:lastrelativeorbitnumber] = convert( Int64, products[:lastrelativeorbitnumber] )
    end
    products[:polygon] = nothing
    products[:gmlfootprint] = nothing
    products[:title] = nothing
    products[:filename] = nothing 
    products[:summary] = nothing
  
    products[:platformname] = nothing
    #products[:platformname] = convert( Int64, products[:platformname] )
    products[:instrumentname] = nothing
    #products[:instrumentname] = convert( Int64, products[:instrumentname] )
    products[:instrumentshortname] = platform[10:12]
    products[:instrumentshortname] = convert( Int64, products[:instrumentshortname] ) 

    
    if "status" in names(products) 
      products[:status] = CategoricalArray(products[:status])
      products[:status] = as.numeric(products[:status]) #ARCHIVED
    end
    ## icon just add to url.alt + Products('Quicklook')/[:value]
    ## url just add to url.alt + [:value
    products[:url] = nothing
    products[:(url_icon)] = nothing
    products[:(url_alt)] = nothing
     
    products[:month] = convert( Int64, Dates.format( Date(products[:beginposition]), "%m") )
    tmp_size = split( products[:size], " ", keepempty = false )
    products[:size] = map( x -> { 
      ret = missing
      if x[2] == "MB"
        ret = as.numeric(x[1])
      end
      if x[2] == "GB"
        ret = as.numeric(x[1]) * 1000 
      end
      if x[2] == "KB"
        ret = as.numeric(x[1]) / 1000
      end
      return ret
      }, tmp_size )
    
    error = nothing
    nempty = findall( occursin.( "EMPTY", products[:footprint] ) )
    if length(nempty) > 0
      error = "$(length(nempty)) righe senza footprint"
      println(error)
      products = products[ -nempty, : ]
    end
    
    dt_data = map( x -> {
 
      fp = products[x,:footprint] 

      tgt = ArchGDAL.fromWKT(fp)
      tg = tgt 
      if length( tg.polygons[1].Polygons ) != 1
        dat = map( xx -> unique( xx.coords ), tg.polygons[1].Polygons )
        dat2 = vcat( dat )
        map( xx -> {  
          if dat2[xx,1] < 0
            dat2[xx,1] += 360
          end
        }, 1:nrow(dat2) )
        
        ch = chull( dat2 ) 
        coords = dat2[ [ch, ch[1]], : ]  # closed polygon  
        pl = Polygons(list(Polygon(coords)), ID=1) 
      else 
        pl = tg.polygons[1]
      end
      
      pl.ID = string( products[x,"uuid"] )
      return pl
    }, products[:,1] )
    
    #products
    #ftp = products[:footprint
    products[:footprint] = nothing
   # dt.data = plyr::compact(dt.data)
    dt_data_sp = SpatialPolygons(dt_data, proj4string=CRS("+proj=longlat +datum=WGS84 +pm=-360") )
    
    insert!( products, 1, :uuid => products[:uuid] )
    
    dt_data_sp_df =  SpatialPolygonsDataFrame( dt_data_sp, products )  
    outins = pgInsert( conn, partial.match = true, name = ["eo_sathub",platforms_tables[platform]], data.obj = dt_data_sp_df, upsert.using = ["uuid","identifier"] )
    if !outins  
      error = "ERRORE in upsert conn"  
      println(error)
    end
    outins = pgInsert(conn2, partial.match = true, name = ["eo_sathub",platforms_tables[platform]], data.obj = dt_data_sp_df, upsert.using = ["uuid","beginposition"] )
    if !outins   
      error = "ERRORE in upsert conn2"  
      println(error)
    end
    
    log_con = open( platform*".log", open="a" )
    if !isnothing(error)
      write( log_con, sprintf("last problem = %s" , error) )
    else 
      CSV.write( platforms_lastdate, "platforms_lastdate.csv" )
    end
    
    #print(platforms_lastdate[platform])
    write( log_con, string.(platforms_lastdate[platform])[1] )
    write( log_con, "\n" )
    close(log_con) 
    
    println("Finito data ===")
    println(platforms_lastdate[platform])
  end
  
  println("Finito =====")
  println(platform)
  
end

#dbGetQuery(conn, "VACUUM FULL ANALYZE eo_sathub.s1")
#dbGetQuery(conn, "VACUUM FULL ANALYZE eo_sathub.s2")
#dbGetQuery(conn, "REFRESH MATERIALIZED VIEW eo_sathub.veneto_subset WITH DATA")


#dbGetQuery(conn2, "VACUUM FULL ANALYZE eo_sathub.s1")
#dbGetQuery(conn2, "VACUUM FULL ANALYZE eo_sathub.s2")
#dbGetQuery(conn2, "REFRESH MATERIALIZED VIEW eo_sathub.veneto_subset WITH DATA")

for i in dbListConnections(pg)
  dbDisconnect(i)
end

end # module