module  Server

using CSV, DataFrames, Dates, LibPQ, Printf

export server #function

# Define server logic required to draw a histogram
function server(input, output, session)
  ##our very own getSentinel modified!
  download_script_selected = nothing
 
  
  output[:mymap] = renderLeaflet({
    leaflet() |>  
      onRender( "function(el, x) {
                 mapElement=this;
                 this.on('layeradd', onLayerAddFunction);
                 initAutocomplete(this); }" )
      |>
      addProviderTiles( providers[:Stamen_TonerLite], options = providerTileOptions(noWrap = true), group = "Toner Lite" )
      |>
      addProviderTiles( providers[:OpenStreetMap], group = "Streets" )
      |>
      addLayersControl( baseGroups = ["Streets", "Toner Lite"], options = layersControlOptions(collapsed = false) )
      |> 
      addMouseCoordinates(style = "basic")
  })
  
  observeEvent( input[:platform], {
    if input[:platform] == ""
      return nothing
    end
      
    #  if input[:platform == "s1"]
    #    return nothing
    #  end

#+ postresql://"user":"password"@"host"/dbname
#+ user = "****", password = "*********", host = "localhost", dbname = "esa"
    con = LibPQ.Connection("postgresql://****:*********@localhost/esa")
    queryst = "SELECT count(*) as n, min(beginposition) as min , max(beginposition) as max from eo_sathub.$(platforms[:tables][ input[:platform] ])"
    

    queryst2 = "SELECT count(*) as n from eo_sathub.$(platforms_tables[ input[:platform] ])"

    df_postgres = execute(con,  queryst)
    updateDateRangeInput(session, "daterange", start = df[:postgres][:min], nd = df[:postgres][:max] )
    close(con)
      
      # println(input[:mymap_zoom])
      # if input[:mymap_zoom] > 5
      #   leafletProxy("mymap") |>  
      #     clearGroup("Sentinel Tiles")
      #     |>  
      #     addWMSTiles( baseUrl = "http://cartografia01.cogeme.com/cgi-bin/mapserv?map=/srv/shiny-server/esa.explorer/mapserver_wms.map",
      #                  options = WMSTileOptions(format = "image/png", transparent = true ),
      #                  layers = platforms[:tables][ input[:platform] ], group = "Sentinel Tiles" )
      #     |>
      #     addLayersControl( baseGroups = ["Streets", "Toner Lite"],
      #                       overlayGroups = ["Sentinel Tiles"],
      #                       options = layersControlOptions(collapsed = false) )  
      # else
      #   leafletProxy("mymap") |>
      #     clearGroup("Sentinel Tiles")
      #     |>
      #     addLayersControl( baseGroups = ["Streets", "Toner Lite"],
      #                       options = layersControlOptions(collapsed = false) )    
      # end

  })
  
  observe({ 
    #req(input[:mytable_rows_selected])
    dt = unique( [ input[:mytable_rows_selected], input[:mytable_row_last_clicked] ] )
    dt = input[:mytable_row_last_clicked]
    # highlightLayer("%s");
    for i in dt
      shinyjs.runjs( "highlightLayer('$(dt_data_sp_df[:data][ i, columns_lut[[ input[:platform] ]][ ["rownames"] ] ])')" )
    end
  })
  
  
  
  
  observeEvent(input[:search], {
    search_bounds = input[:mymap_bounds] 
 
    if isnothing(input[:platform]) || input[:platform] == ""
      sendSweetAlert( session = session, text = "Please choose platform", html = true ) 
      return nothing
    end
    if isnothing(input[:password]) || input[:password] == ""
      sendSweetAlert( session = session, "",
                      text = htmltools.HTML("Please provide password to ESA SciHub portal 
                      (register <a href=\"https://scihub.copernicus.eu/dhus/#/self-registration\" 
                      target=\"_blank\">here</a> if you are not registered)!"), html = true )  
      return nothing     
    end
    if isnothing(input[:user]) || input[:user] == ""          
      sendSweetAlert( session = session, "",
                      text = HTML("Please provide username to ESA SciHub portal 
                      (register <a href=\"https://scihub.copernicus.eu/dhus/#/self-registration\" 
                      target=\"_blank\">here</a> if you are not registered)!"), html = true)    
      return nothing          
    end 
    withProgress( message = "Querying Database",
                  detail = "This may take a while...", value = 0, {

                   extras = []
                   println(input[:orbit_dir])
                   println(input[:day])
                   
                   if !isnothing(input[:month]) && input[:month] != ""
                     push!( extras, "month in ($( join(input[:month], ", ") ))" )
                   end
                   if !isnothing(input[:day]) && input[:relativeorbitnumber] != ""
                    push!( extras, "relativeorbitnumber in ($( join(input[:relativeorbitnumber], ", ") ))" )
                   end 
                   if !isnothing(input[:day]) && input[:day] != ""
                    push!( extras, "day in ('$( join(input[:day], "', '") )')" )
                   end
                   if !isnothing(input[:orbit_dir]) && input[:orbit_dir] != ""
                    push!( extras, "orbitdirection in ('$( join(input[:orbit_dir], "', '") )')" )
                   end
                   if platforms_tables[ input[:platform] ] == "s1"
                     if !isnothing( input[:sensoroperationalmode] ) && input[:sensoroperationalmode] != ""
                      push!( extras, "sensoroperationalmode in ('$( join(input[:sensoroperationalmode], "', '") )')" )
                     end
                     if !isnothing( input[:polarisationmode] ) && input[:polarisationmode] != ""
                      push!( extras, "polarisationmode in ('$( join(input[:polarisationmode], "', '") )')" ) 
                     end
                     if !isnothing( input[:producttype] ) && input[:producttype] != ""
                      push!( extras, "producttype in ('$( join(input[:producttype], "', '") )')" )
                     end
                   end
                   
                   if platforms_tables[ input[:platform] ] == "s2"
                     if !isnothing( input[:cloudcoverpercentage] ) && input[:cloudcoverpercentage] != ""
                       extras = [ extras, "cloudcoverpercentage < $(input[:cloudcoverpercentage])" ]
                     end
                   end

                  
                   con = LibPQ.Connection("postgresql://****:*********@localhost/esa")
                   
                   extras = length(extras) > 0 ? " AND $( join( extras, " AND " ) )" : ""
                   
                   queryst = "WHERE  beginposition > '$(input[:daterange][1])'::timestamp AND
                                     beginposition < '$(input[:daterange][2])'::timestamp AND
                                     st_intersects(geom,  ST_MakeEnvelope ( $(input[:mymap_bounds][:west]),
                                                                            $(input[:mymap_bounds][:south]),
                                                                            $(input[:mymap_bounds][east]),
                                                                            $(input[:mymap_bounds][:north]), 4326) ) $extras limit 1000"
                      
                   println(queryst)
                   try 
                      dt_data_sp_df = pgGetGeom( con, ["eo_sathub", platforms_tables[ input[:platform] ] ],
                                                       other_cols = names( columns_lut[ input[:platform] ]["lut"]), 
                                                       clauses = queryst )
                    catch e
                      dt_data_sp_df
                    end
                   
                   close(con)
 
                   if isnothing(dt_data_sp_df) || nrow(dt_data_sp_df[:data]) == 0
                     sendSweetAlert( session = session, "No product found for query", "No $(input[:platform]) product found!" )
                     return nothing
                   end
                    
                   rename!( dt_data_sp_df[:data], columns_lut[ input[:platform] ]["lut"][:,1] )
                   
                   insert!( dt_data_sp_df[:data], 1, :Rownames => dt_data_sp_df[:data][:, columns_lut[ input[:platform] ]["rownames"] ] )
                   
                   map( function (x)
                          if !ismissing(x[2])
                            if x[2][1:3] == "as." 
                              dt_data_sp_df[:data][ x[1] ] = get( x[2] )( dt_data_sp_df[:data][ x[1] ] )
                            else
                              dt_data_sp_df[:data][ x[1] ] = get( x[2] )( dt_data_sp_df[:data] )
                            end
                          end
                        end, columns_lut[ input[:platform] ]["lut"] )
                   updateTabsetPanel( session, "mytabset", selected = "tabletab" )
                   output[:mytable] = DT.renderDataTable({
                     DT.datatable(
                       dt_data_sp_df[:data],
                       filter = "top",
                       rownames = false,
                       class = "compact hover stripe",
                       escape = false,
                       options = Dict( :columnDefs => [ Dict( :visible => false,
                                                              :targets => findall( columns_lut[ input[:platform] ]["lut"][:,3] .!= true ) - 1 ) ],
                                       :initComplete => JS(
                                         "function(settings, json) {",
                                         "$(this.api().table().header()).css({'background-color': '#949DA6', 'color': '#000'});",
                                         "}" ) )
                     )
                   })
                 })
  })
  
  observeEvent( input[:mytable_rows_selected], { 
    s = input[:mytable_rows_selected]
    download_script_selected = []
    if length(s)
      
      map( x -> push!( download_script_selected,  sprintf( download_script,
                                                           input[:user],
                                                           input[:password],
                                                           dt_data_sp_df[:data][x,"uuid"],
                                                           dt_data_sp_df[:data][x,"uuid"] ) ), s )
      println("These rows were selected:\n\n")
      bounds = "[  [ $(search_bounds[:west]), $(search_bounds[:east]) ], [ $(search_bounds[:south]), $(search_bounds[:north]) ]  ]"
      
      nm = "&nbsp;<i onclick=\"mapElement.bounds($bounds)\" class=\"fa fa-search\"></i>&nbsp;Search area"
      
 
      leafletProxy("mymap", data = dt_data_sp_df[s,:]) |>
        clearShapes()
        |>  
        fitBounds( dt_data_sp_df[s,:][:bbox]["x","min"], dt_data_sp_df[s,:][:bbox]["y","min"],
                   dt_data_sp_df[s,:][:bbox]["x","max"], dt_data_sp_df[s,:][:bbox]["y","max"] )
        |>
        addPolygons( layerId = dt_data_sp_df[s,:][:data][:,1], 
                     weight = 1,
                     color = "red",
                     fillOpacity = 0.0,
                     opacity = 0.5,
                     highlightOptions = highlightOptions( color = "yellow", 
                                                          weight = 5,
                                                          opacity=1.0,
                                                          bringToFront = true ),
                     popup = "<img width=\"200\" src=\"https://scihub.copernicus.eu/dhus/odata/v1/Products('$uuid')/Products('Quicklook')/\"><br><a href=\"https://scihub.copernicus.eu/dhus/odata/v1/Products('$uuid')/\$value\" target=\"_blank\">Download</a><br>Size:$size",
                     group = "Footprints" )
        |>
        addRectangles( lng1 = search_bounds[:west], lat1 = search_bounds[:south],
                       lng2 = search_bounds[:east], lat2 = search_bounds[:north],
                       fillColor = "transparent", group = nm ) |>
            addLayersControl(
              baseGroups = ["Streets", "Toner Lite"],
              overlayGroups = ["Footprints", nm],
              options = layersControlOptions(collapsed = false) )
    end
  })

  function filterProducts()
    products_filtered1 = products[ findall( products[:processinglevel] .== "Level-1C" ), : ] #filter by Level
    products_filtered2 = products_filtered1[ products_filtered1[:cloudcoverpercentage] .<= 30, ] #filter by clouds
  end

  output[:downloadData] = downloadHandler(
    filename = () -> return "data-$(today()).sh",
    content = file -> CSV.write( file, download_script_selected, writeheader = false, quotestrings = false )
  )

end

end # module