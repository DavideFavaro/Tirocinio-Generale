source("getSpatialData_dev.R")
library(sf)
library(sp)
library(rgeos)
library(rpostgis) 
library(rgdal)

pg = dbDriver("PostgreSQL")
for(i in dbListConnections(pg)) dbDisconnect(i)
conn <- dbConnect(pg, dbname = "esa", port= 27017, host = "localhost", user = "****", password = "*********")
conn2 <- dbConnect(pg, dbname = "esa", port= 5432, host = "geolab.vs-ix.net", user = "****", password = "*********")

platforms.startdate<-list("Sentinel-1"=as.Date("2024-04-01"), 
     "Sentinel-2"=as.Date("2015-06-01") )#, 
   #  "Sentinel-3"=as.Date("2016-02-01") )

platforms.tables<-list("Sentinel-1"="s1", 
                          "Sentinel-2"="s2", 
                          "Sentinel-3"="s3" )


platforms.lastdate<-readRDS("platforms.lastdate.rds") 
#platforms.lastdate$`Sentinel-1`<-platforms.lastdate$`Sentinel-1`-2
#platforms.lastdate$`Sentinel-2`<-platforms.lastdate$`Sentinel-2`-2

 ## 1 day  timespan
timespan<- 1
platform<-"Sentinel-2"  
login_CopHub(username = "fpirotti", password="libero") #asks for password or define 'password'
set_archive("/tmp")

reloaday<-F

for(platform in names(platforms.startdate))
{  
   
  if(reloaday)
  {
    log_con <- file(paste0(sep="", platform,".log"), open="r")
    lns<-readLines(log_con)
    names(lns)<- gsub("problem", "", lns)
    close(log_con) 
     
    ww<-grep("problem", lns)
    
    
    for(date in names(lns[ww]) )
    {
      print(date)
      date<-as.Date(date)
      time_range =  c(as.character(date-timespan) ,
                      as.character(date) )
      products <- getSentinel_query2(time_range=time_range,  platform = platform)
      products.tmp<-products 
      
      if(is.null(products) || nrow(products)<1 ){
        
        lns[[as.character(date)]]<-paste0( collapse=" -- 0 righe ",as.character(date))  
        
        log_con <- file(paste0(sep="", platform,".log"), open="w")
        lns.out<-writeLines(paste0(unlist(lns)), log_con) 
        close(log_con) 
        
        next
      }
      
      
      
      products<-products[ !is.na(products$footprint), ]
      rownames(products)<-products$uuid 
      products$format<-as.factor(products$format)
      products[ products$orbitdirection == 'ASCENDING', "orbitdirection"]<-'ASC'
      products[ products$orbitdirection == 'DESCENDING', "orbitdirection"]<-'DESC'
      
      products$orbitnumber<-as.integer(products$orbitnumber)
      if(  "lastorbitnumber" %in% names(products) )  products$lastorbitnumber<-as.integer(products$lastorbitnumber)
      
      if(  "lastrelativeorbitnumber" %in% names(products) ) products$lastrelativeorbitnumber<-as.integer(products$lastrelativeorbitnumber)
      products$polygon<-NULL
      products$gmlfootprint<-NULL
      products$title<-NULL
      products$filename<-NULL 
      products$summary<-NULL
      
      products$platformname<-NULL
      #products$platformname<-as.integer(products$platformname)
      products$instrumentname<-NULL
      #products$instrumentname<-as.integer(products$instrumentname)
      products$instrumentshortname<-substr(platform,10,12)
      products$instrumentshortname<-as.integer(products$instrumentshortname) 
      
      
      if(  "status" %in% names(products) )  { 
        products$status<-factor(products$status)
        products$status<-as.numeric(products$status) #ARCHIVED
      }
      ## icon just add to url.alt + Products('Quicklook')/$value
      ## url just add to url.alt + $value
      products$url<-NULL
      products$url.icon<-NULL
      products$url.alt<-NULL
      
      products$month<-format(as.Date(products$beginposition), "%m") 
      tmp.size<-stringr::str_split(products$size," ")
      products$size<- sapply(tmp.size, function(x){ 
        ret<-NA
        if(x[[2]]=='MB') ret<-as.numeric(x[[1]]) 
        if(x[[2]]=='GB') ret<-as.numeric(x[[1]])*1000 
        if(x[[2]]=='KB') ret<-as.numeric(x[[1]])/1000
        ret
      })
      dt.data<-lapply(rownames(products), function(x){
        tg<-readWKT(products[x,"footprint"])
        pl<-tg@polygons[[1]]
        pl@ID<-as.character(products[x,"uuid"])
        pl
      })
      
      #ftp<-products$footprint
      products$footprint<-NULL
      
      dt.data.sp<-SpatialPolygons(dt.data, proj4string=CRS("+proj=longlat +datum=WGS84") )
      dt.data.sp.df<- SpatialPolygonsDataFrame(dt.data.sp, products)
      nPolys <- sapply(dt.data.sp@polygons, function(x)length(x@Polygons))
      dd<-dt.data.sp.df[ which(nPolys>1),] 
      outins<-pgInsert(conn, partial.match = T, name=c("eo_sathub",platforms.tables[[platform]]), 
                       data.obj=dt.data.sp.df, upsert.using="uuid,identifier" )
      outins<-pgInsert(conn2, partial.match = T, name=c("eo_sathub",platforms.tables[[platform]]), 
                       data.obj=dt.data.sp.df, upsert.using="uuid,identifier" )
      
      if(outins) {  lns[[as.character(date)]]<-as.character(date)      }
      print(date)
    
      
      
      log_con <- file(paste0(sep="", platform,".log"), open="w")
      lns.out<-writeLines(paste0(unlist(lns)), log_con) 
      close(log_con) 
      
      
    }
  }
  
   
   
  
  
  
   
  
  
  
  ##remove multi polygon!!
 # sql<-"SELECT   beginposition  FROM eo_sathub.s2 where ST_NumGeometries(geom)>1 group by beginposition order by beginposition"
  
  #redo<-dbGetQuery(conn, sql)
  #for(date in redo$beginposition)
  #Sys.date - 1 per evitare ingestion / beginposition differenza
  while( platforms.lastdate[[platform]] < (Sys.Date()) ) 
   {
    
    print(platform)
    
    print(platforms.lastdate[[platform]])
    
    date<-platforms.lastdate[[platform]]
    
    time_range =  c(as.character(date) , as.character(date+timespan) )
    
    # products <- getSentinel_query2(time_range=time_range,  platform = platform)

    #, extra="processinglevel:Level-2Ap"
    
    products <- getSentinel_query2(time_range=time_range,  platform = platform,
                                   aoi = NULL, username = NULL, 
                                   password = NULL, hub = "auto", verbose = TRUE, ingestiondate = T)
    
 
    #saveRDS(products,  "products.rds") 
    products.tmp<-products
    platforms.lastdate[[platform]]<-(platforms.lastdate[[platform]]+timespan) 
    
    if(is.null(products) || nrow(products)<1 ){
      saveRDS(platforms.lastdate,"platforms.lastdate.rds")
      
     # log_con <- file(paste0(sep="", platform,".log"), open="r")
      cat( paste0(as.character(platforms.lastdate[[platform]])[[1]], " problem0") , file = paste0(sep="", platform,".log")    )
      print(platforms.lastdate[[platform]])
      next
    }
    
    products<-products[ !is.na(products$footprint), ]
    products<-products[ !duplicated(products), ]
    rownames(products)<-products$uuid 
    
    #uuids<-dbGetQuery(conn, "select uuid from eo_sathub.s2 where processinglevel='Level-2Ap' ")
    #tokeep<-  which( !(rownames(products) %in% uuids$uuid)  )
    
    products$format<-as.factor(products$format)
    products[ products$orbitdirection == 'ASCENDING', "orbitdirection"]<-'ASC'
    products[ products$orbitdirection == 'DESCENDING', "orbitdirection"]<-'DESC'
    
    products$orbitnumber<-as.integer(products$orbitnumber)
    if(  "lastorbitnumber" %in% names(products) )  products$lastorbitnumber<-as.integer(products$lastorbitnumber)
    
    if(  "lastrelativeorbitnumber" %in% names(products) ) products$lastrelativeorbitnumber<-as.integer(products$lastrelativeorbitnumber)
    products$polygon<-NULL
    products$gmlfootprint<-NULL
    products$title<-NULL
    products$filename<-NULL 
    products$summary<-NULL
  
    products$platformname<-NULL
    #products$platformname<-as.integer(products$platformname)
    products$instrumentname<-NULL
    #products$instrumentname<-as.integer(products$instrumentname)
    products$instrumentshortname<-substr(platform,10,12)
    products$instrumentshortname<-as.integer(products$instrumentshortname) 

    
    if(  "status" %in% names(products) )  { 
      products$status<-factor(products$status)
      products$status<-as.numeric(products$status) #ARCHIVED
    }
    ## icon just add to url.alt + Products('Quicklook')/$value
    ## url just add to url.alt + $value
    products$url<-NULL
    products$url.icon<-NULL
    products$url.alt<-NULL
     
    products$month<-as.integer(format(as.Date(products$beginposition), "%m") )
    tmp.size<-stringr::str_split(products$size," ")
    products$size<- sapply(tmp.size, function(x){ 
      ret<-NA
      if(x[[2]]=='MB') ret<-as.numeric(x[[1]]) 
      if(x[[2]]=='GB') ret<-as.numeric(x[[1]])*1000 
      if(x[[2]]=='KB') ret<-as.numeric(x[[1]])/1000
      ret
      })
    
    error<-NULL
    nempty<-which(grepl("EMPTY", products$footprint, fixed = T))
    if(length(nempty)>0){
      error<-sprintf("%d righe senza footprint", length(nempty))  
      print(error)
      products<-products[-nempty,]
    }
    
    dt.data<-lapply(rownames(products), function(x){
 
      fp<-products[x,"footprint"] 

      tgt<-readWKT(fp)
      tg<-tgt 
      if(length(tg@polygons[[1]]@Polygons)!=1){
        dat<-lapply(tg@polygons[[1]]@Polygons, function(xx){ unlist(unique(xx@coords)) }) 
        dat2<-do.call(rbind,dat) 
        sapply(1:nrow(dat2), function(xx) {  
          if(dat2[xx,1]<0) dat2[xx,1]<<- dat2[xx,1]+360
        } )
        
        ch <- chull( (dat2) ) 
        coords <- dat2[c(ch, ch[1]), ]  # closed polygon  
        pl <-  Polygons(list(Polygon(coords)), ID=1) 
      }
      else { 
        pl<-tg@polygons[[1]]
      }
      
      pl@ID<-as.character(products[x,"uuid"])
      pl
    })
    
    #products
    #ftp<-products$footprint
    products$footprint<-NULL
   # dt.data<-plyr::compact(dt.data)
    dt.data.sp<-SpatialPolygons(dt.data, proj4string=CRS("+proj=longlat +datum=WGS84 +pm=-360") )
    
    rownames(products)<-products$uuid
    
    dt.data.sp.df<- SpatialPolygonsDataFrame(dt.data.sp, products )  
    outins<-pgInsert(conn, partial.match = T, name=c("eo_sathub",platforms.tables[[platform]]), data.obj=dt.data.sp.df,  upsert.using=c("uuid","identifier")  )
    if(!outins) {   
      error<-"ERRORE in upsert conn"  
      print(error)
      }
    outins<-pgInsert(conn2, partial.match = T, name=c("eo_sathub",platforms.tables[[platform]]), data.obj=dt.data.sp.df, upsert.using=c("uuid","beginposition") )
    if(!outins) {   
      error<-"ERRORE in upsert conn2"  
      print(error)
      }  
    
    log_con <- file(paste0(sep="", platform,".log"), open="a")
    if(!is.null(error)) {    cat( sprintf("last problem = %s" , error), file = log_con)    }
    else saveRDS(platforms.lastdate,"platforms.lastdate.rds")
    
    #print(platforms.lastdate[[platform]])
    cat( as.character(platforms.lastdate[[platform]])[[1]] , file = log_con)
    cat( "\n" , file = log_con)
    close(log_con) 
    
    print("Finito data ===")
    print(platforms.lastdate[[platform]])
  }
  
  print("Finito =====")
  print(platform)
  
}

#dbGetQuery(conn, "VACUUM FULL ANALYZE eo_sathub.s1")
#dbGetQuery(conn, "VACUUM FULL ANALYZE eo_sathub.s2")
#dbGetQuery(conn, "REFRESH MATERIALIZED VIEW eo_sathub.veneto_subset WITH DATA")


#dbGetQuery(conn2, "VACUUM FULL ANALYZE eo_sathub.s1")
#dbGetQuery(conn2, "VACUUM FULL ANALYZE eo_sathub.s2")
#dbGetQuery(conn2, "REFRESH MATERIALIZED VIEW eo_sathub.veneto_subset WITH DATA")

for(i in dbListConnections(pg)) dbDisconnect(i)