
# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  ##our very own getSentinel modified!
  download.script.selected<-NULL
 
  
  output$mymap <-  renderLeaflet({
    leaflet() %>%  
      onRender("function(el, x) {
               mapElement=this;
               this.on('layeradd', onLayerAddFunction);
               initAutocomplete(this); 
  }") %>%
      addProviderTiles(providers$Stamen.TonerLite,
                       options = providerTileOptions(noWrap = TRUE), group = "Toner Lite" ) %>%
      addProviderTiles( providers$OpenStreetMap, group = "Streets" ) %>%
      addLayersControl(
        baseGroups = c("Streets",   "Toner Lite"), 
        options = layersControlOptions(collapsed = FALSE) )  %>% 
      addMouseCoordinates(style = "basic")
})
  
  observeEvent( input$platform, {
    if(input$platform=="") return(NULL);
    
  #  if(input$platform=="s1") return(NULL);
    
    
    
    
    
    
    
    
    con <- dbConnect("PostgreSQL", 
                     dbname = "esa", 
                     host = "localhost", 
                     user = "****", 
                     password = "*********")
    
    queryst<-sprintf("SELECT count(*) as n, min(beginposition) as min , max(beginposition) 
                     as max from eo_sathub.%s",platforms.tables[[ input$platform ]] )
    
    
    queryst2<-sprintf("SELECT count(*) as n from eo_sathub.%s",platforms.tables[[ input$platform ]] )
    
    df_postgres <- dbGetQuery(con,  queryst)
    updateDateRangeInput(session, 'daterange', start = df_postgres$min, end = df_postgres$max )
    dbDisconnect(con)
    
    
    # print(input$mymap_zoom)
    # if(input$mymap_zoom > 5) {
    #   leafletProxy('mymap') %>%  
    #     clearGroup("Sentinel Tiles") %>%  
    #     addWMSTiles(baseUrl = "http://cartografia01.cogeme.com/cgi-bin/mapserv?map=/srv/shiny-server/esa.explorer/mapserver_wms.map" ,
    #                 options = WMSTileOptions(format = "image/png", transparent = T),
    #                 layers=platforms.tables[[ input$platform ]], group = "Sentinel Tiles" ) %>%
    #     addLayersControl(
    #       baseGroups = c("Streets",   "Toner Lite"),
    #       overlayGroups = c("Sentinel Tiles"),
    #       options = layersControlOptions(collapsed = FALSE) )  
    # }
    # else {
    #   leafletProxy('mymap') %>%   clearGroup("Sentinel Tiles")  %>%
    #     addLayersControl(
    #       baseGroups = c("Streets",   "Toner Lite"),
    #       options = layersControlOptions(collapsed = FALSE) )    
    # }
    
    
    
  })
  
  # observe({ 
  #   #req(input$mytable_rows_selected)
  #   dt<- unique(c(input$mytable_rows_selected, input$mytable_row_last_clicked))
  #   dt<-input$mytable_row_last_clicked
  #   # highlightLayer('%s');
  #   for(i in dt ) shinyjs::runjs( sprintf("highlightLayer('%s')",  dt.data.sp.df@data[ i, columns.lut[[input$platform]][['rownames']] ] ))
  # })
  
  
  
  
  observeEvent(input$search, {
    search.bounds<<-input$mymap_bounds 
 
      if(is.null(input$platform) || input$platform=="") {
           sendSweetAlert(session=session, text="Please choose platform", html = T) 
            return(NULL)
            }
   
      if(is.null(input$password) || input$password==""){ 
           sendSweetAlert(session=session, "",
                          text=htmltools::HTML('Please provide password to ESA SciHub portal 
                          (register <a href="https://scihub.copernicus.eu/dhus/#/self-registration" 
                          target="_blank">here</a> if you are not registered)!'), html = TRUE)  
        return(NULL)
           
      }
    if(is.null(input$user) || input$user==""){            
           sendSweetAlert(session=session, "",
                          text=  HTML('Please provide username to ESA SciHub portal 
                          (register <a href="https://scihub.copernicus.eu/dhus/#/self-registration" 
                          target="_blank">here</a> if you are not registered)!'), html = TRUE)   
      return(NULL)
           
           
    } 
    withProgress(message = 'Querying Database',
                 detail = 'This may take a while...', value = 0, {
                   
                   extras<-c()
                   print(input$orbit.dir )
                   print(input$day )
                   
                   if(!is.null(input$month) && input$month!="") {
                     extras<-c(extras, sprintf("month in (%s)", 
                                               paste(input$month, collapse=",", sep="")) ) 
                   }
                   
                   if(!is.null(input$day) && input$relativeorbitnumber!="") {
                     extras<-c(extras, sprintf("relativeorbitnumber in (%s)", 
                                               paste(input$relativeorbitnumber, collapse=",", sep="")) ) 
                   } 
                   if(!is.null(input$day) && input$day!="") {
                     extras<-c(extras, sprintf("day in ('%s')", 
                                               paste(input$day, collapse="','", sep="")) ) 
                   } 
                   
                   if(!is.null(input$orbit.dir) && input$orbit.dir!="") {
                     extras<-c(extras, sprintf("orbitdirection in ('%s')", 
                                               paste(input$orbit.dir, collapse="','", sep="")) ) 
                   } 
                   
                   if(platforms.tables[[ input$platform ]]=='s1'){
                     if(!is.null( input$sensoroperationalmode ) && input$sensoroperationalmode!="") {
                       extras<-c(extras, sprintf("sensoroperationalmode in ('%s')", 
                                                 paste(input$sensoroperationalmode, collapse="','", sep="")) ) 
                     } 
                     if(!is.null( input$polarisationmode ) && input$polarisationmode!="") {
                       extras<-c(extras, sprintf("polarisationmode in ('%s')", 
                                                 paste(input$polarisationmode, collapse="','", sep="")) ) 
                     } 
                     if(!is.null( input$producttype ) && input$producttype!="") {
                       extras<-c(extras, sprintf("producttype in ('%s')", 
                                                 paste(input$producttype, collapse="','", sep="")) ) 
                     }
                   }
                   
                   if(platforms.tables[[ input$platform ]]=='s2'){
                     if(!is.null( input$cloudcoverpercentage ) && input$cloudcoverpercentage!="") {
                       extras<-c(extras, sprintf("cloudcoverpercentage < %f ", 
                                                 input$cloudcoverpercentage ) ) 
                     } 
 
                   }
                   con <- dbConnect("PostgreSQL", 
                                    dbname = "esa", 
                                    host = "localhost", 
                                    user = "****", 
                                    password = "*********")
                   
                   if(length(extras) > 0 ) extras<-paste(sep="", " AND ",  paste(extras, collapse=" AND ", sep="")  )
                   else extras <-  ""
                   
                   queryst<-sprintf("WHERE  beginposition > '%s'::timestamp AND beginposition < '%s'::timestamp  AND st_intersects(geom,  ST_MakeEnvelope (   %f, %f,  %f, %f,   4326) ) %s limit 1000",
                                    input$daterange[[1]], input$daterange[[2]],
                                    input$mymap_bounds$west,  input$mymap_bounds$south,
                                    input$mymap_bounds$east,  input$mymap_bounds$north,
                                    extras 
                                    )

                                                       
                   print(queryst)
                   tryCatch( 
                      dt.data.sp.df <<- pgGetGeom(con, c('eo_sathub', platforms.tables[[ input$platform ]] ),
                                         other.cols = rownames(columns.lut[[input$platform]][['lut']]), 
                                         clauses=queryst),
                      error=function(e) { dt.data.sp.df<<-NULL   } )
                   
                   
                   dbDisconnect(con)
 
                   if(is.null(dt.data.sp.df) || nrow(dt.data.sp.df@data)==0)
                   { 
                     sendSweetAlert(session=session, "No product found for query", sprintf("No %s product found!", input$platform) )
                     return(NULL) 
                   }
                    
                   names(dt.data.sp.df@data)<-columns.lut[[input$platform]][['lut']][,1]
                   
                   rownames(dt.data.sp.df@data) <- dt.data.sp.df@data[, columns.lut[[input$platform]][['rownames']] ]
                   
                   apply(columns.lut[[input$platform]][['lut']], 1, function(x){
                     if(!is.na( x[[2]] )) {
                       if( substr(x[[2]],1,3)=="as." ) 
                         dt.data.sp.df@data[[ x[[1]] ]]<<-get( x[[2]] )( dt.data.sp.df@data[[ x[[1]] ]] )
                       else
                         dt.data.sp.df@data[[ x[[1]] ]]<<-get( x[[2]] )( dt.data.sp.df@data )
                       
                     }
                   })
                   updateTabsetPanel(session, "mytabset", selected = "tabletab")
                   output$mytable = DT::renderDataTable({
                     DT::datatable(
                       dt.data.sp.df@data, filter = "top",
                       rownames= FALSE,   
                       class = 'compact hover stripe', escape = FALSE,
                       options = list(columnDefs = 
                                        list(
                                          list(visible=FALSE, 
                                                  targets= which( !as.logical(columns.lut[[ input$platform ]][['lut']][,3]) )-1 )),
                                      initComplete = JS(
                                        "function(settings, json) {",
                                        "$(this.api().table().header()).css({'background-color': '#949DA6', 'color': '#000'});",
                                        "}")
                       )
                     )
                   })
                 })
    
  })
 
  observeEvent(input$mytable_rows_selected, { 
    s = input$mytable_rows_selected
    download.script.selected<<-c()
    if (length(s)) {
      
      lapply(s, function(x){
        download.script.selected<<- c(download.script.selected, 
                        sprintf(download.script, input$user, input$password,dt.data.sp.df@data[x,"uuid"], dt.data.sp.df@data[x,"uuid"] )              
                                      )
      })
      cat('These rows were selected:\n\n') 
      bounds<-paste0(sep="", "[  [", search.bounds$west,",", search.bounds$east," ], [", search.bounds$south,",", search.bounds$north," ]  ]")
      
      nm<-paste0("&nbsp;<i onclick=\"mapElement.bounds(",bounds,")\" class=\"fa fa-search\"></i>&nbsp;Search area", sep="")
      
 
      leafletProxy("mymap",data=dt.data.sp.df[s,]) %>%
        clearShapes() %>%  
        fitBounds(dt.data.sp.df[s,]@bbox["x","min"], dt.data.sp.df[s,]@bbox["y","min"],
                  dt.data.sp.df[s,]@bbox["x","max"], dt.data.sp.df[s,]@bbox["y","max"]) %>%
        addPolygons(  layerId = rownames(dt.data.sp.df[s,]@data), 
                      weight = 1, color="red",
                      fillOpacity = 0.0, opacity=0.5,
                      highlightOptions = highlightOptions(color = "yellow", 
                                                          weight = 5, opacity=1.,
                                                          bringToFront = TRUE),
                      popup = ~sprintf("<img width=\"200\" src=\"https://scihub.copernicus.eu/dhus/odata/v1/Products('%s')/Products('Quicklook')/\"><br><a href=\"https://scihub.copernicus.eu/dhus/odata/v1/Products('%s')/$value\" target=\"_blank\">Download</a><br>Size:%s", uuid, uuid, size),
                      group = "Footprints") %>%
        addRectangles(
          lng1=search.bounds$west, lat1=search.bounds$south,
          lng2=search.bounds$east, lat2=search.bounds$north,
          fillColor = "transparent", group = nm
        ) %>%
            addLayersControl(
              baseGroups = c("Streets",   "Toner Lite"),
              overlayGroups = c("Footprints", nm),
              options = layersControlOptions(collapsed = FALSE) )
    }
  })
  
  filterProducts <-function(){
    products_filtered1 <<- products[which(products$processinglevel == "Level-1C"),] #filter by Level
    products_filtered2 <<- products_filtered1[products_filtered1$cloudcoverpercentage <= 30, ] #filter by clouds
  }
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("data-", Sys.Date(), ".sh", sep="")
    },
    content = function(file) {
      write.csv(download.script.selected, file, sep = " ",  row.names = F, col.names = F, quote = F )
    }
  )
  
  
  
  
  }
