module GetSpatialDataDevInternal

using DataFrames

export options, #Options dictionary
       checkCmd, convMODIS_names, copHub_select, EE_ds, EE_preview, EE_query, ERS_login, ERS_logout, ESPA_order, ESPA_download,
       gSD_get, gSD_download, gSD_post, isFalse, isTrue, make_aoi, out #Functions



#' Outputs errors, warnings and messages
#'
#' @param input character
#' @param type numeric, 1 = message/cat, 2 = warning, 3 = error and stop
#' @param msg logical. If \code{TRUE}, \code{message} is used instead of \code{cat}. Default is \code{FALSE}.
#' @param sign character. Defines the prefix string.
#'
#' @keywords internal
#' @noRd 
#'
function out( input; type = 1, ll = nothing, msg = false, sign = "", verbose = options[:gSD_verbose] )
  if isnothing(ll)
    ll = verbose ? 1 : 2
  end
  if type == 2 && ll <= 2
    @warn sign.*input
  else
    if type == 3
      throw( DomainError( input ) )
    else
      if ll == 1
        if !msg
          println(join( sign.*input, "\n" ))
        else
          @info sign.*input
        end
      end
    end
  end
end



#' Simplifies check of variables being FALSE
#'
#' @param evaluate variable or expression to be evaluated
#'
#' @keywords internal
#' @noRd

function isFalse(evaluate)
  #if evaluate == false
  #  return true
  #else
  #  return false
  #end
  return evaluate != true
end



#' Simplifies check of variables being TRUE
#'
#' @param evaluate variable or expression to be evaluated
#'
#' @keywords internal
#' @noRd
function isTrue(evaluate)
  #if evaluate == true
  #  return true
  #else
  #  return false
  #end
  return evaluate == true
end


#' Checks, if specific command is available
#'
#' @param cmd command
#' @importFrom devtools system_check
#' @keywords internal
#' @noRd
function checkCmd(cmd)
  #sc = try 
  #       system_check(cmd, quiet = true) 
  #     catch exc 
  #       if exc isa SystemError
  #         return false
  #       else
  #         return true
  #       end
  #     end
  try
    system_check(cmd, quiet = true) 
  catch exception
      return !isa(exception, SystemError)
  end
end


#' gSD.get
#' @param url url
#' @param username user
#' @param password pass
#' @param dir.file output file path
#' @param prog show or not show progress console
#' @importFrom httr GET stop_for_status warn_for_status message_for_status progress
#' @keywords internal
#' @noRd
function gSD_get(url; username = nothing, password = nothing, dir_file = nothing, prog = false)

  #x <- NULL # needed due to checks
  #get.str <-"x <- GET(url"
  #if(!is.null(username)) get.str <- paste0(get.str, ", authenticate(username, password)")
  #if(!is.null(dir.file)) get.str <- paste0(get.str, ", write_disk(dir.file)")
  #if(is.TRUE(prog)) get.str <- paste0(get.str, ", progress()")
  #get.str <- paste0(get.str, ")")
  #eval(parse(text = get.str))

  #getfield( """QUALUNQUE SIA IL MODULO DELLA FUNZIONE""", Symbol(get_str) )()
  
  x = nothing # needed due to checks
  if !isnothing(username)
    a = authenticate(username, password)
  end
  if !isnothing(dir_file)
    w = write_disk(dir_file)
  end
  if prog
    p = progress()
  end
  x = GET( url, a, w, p )
  stop_for_status(x, "connect to server.")
  warn_for_status(x)
  #message_for_status(x); cat("\n")
  return x
end



#' gSD.post
#' @param url url
#' @param username user
#' @param password pass
#' @param body body
#' @importFrom httr POST stop_for_status warn_for_status message_for_status progress
#' @keywords internal
#' @noRd
function gSD_post( url; username = nothing, password = nothing, body = false )
  
  #post.str <-"x <- POST(url"
  #if(!is.null(username)) post.str <- paste0(post.str, ", authenticate(username, password)")
  #post.str <- paste0(post.str, ", body = body)")
  #eval(parse(text = post.str))
  if !isnothing(username)
    a = authenticate(username, password)
  end
  x = POST( url, a, body = body )
  
  if !isnothing(username)
    x = POST( url, authenticate(username, password), body = body )
  else
    x = POST( url, body = body )
  end
  stop_for_status(x, "connect to server.")
  warn_for_status(x)
  #message_for_status(x); cat("\n")}
  return x
end



#' gSD.download
#' @param name name
#' @param url url
#' @param file file
#' @importFrom tools md5sum
#' @keywords internal
#' @noRd
function gSD_download( name, url_file, file; url_checksum = nothing )
  out( join([ "Attempting to download '", name, "' to '", file, "'..." ]), msg = true )
  file_tmp = tempfile( tmpdir = join( split( file, "/" )[1:end-1], "/" ) ) #, fileext = ".tar.gz") 
  gSD_get( url_file, dir_file = file_tmp, prog = true )
  if isnothing( url_checksum )
    md5 = split( content( gSD_get(url_checksum), as = "text", encoding = "UTF-8" ), " " )[1]
    if string( md5sum(file_tmp) ) == lowercase(md5)
      out( "Successfull download, MD5 check sums match.", msg = true )
    else
      out( "Download failed, MD5 check sums do not match. Will retry.", type = 2 )
      rm(file_tmp)
      return false
    end
  end #else out("Download finished. MD5 check sums not available (file integrity could not be checked).", msg = T)

  mv( file_tmp, file )
  return false
end



#' get Copernicus Hub API url and credentials from user input
#'
#' @param x API keyword or URL
#' @param p platform
#' @param user user name
#' @param pw password
#' @keywords internal
#' @noRd
function copHub_select(x, p, user, pw) #cophub_api
  if x == "auto"
    if p == "Sentinel-1" || p == "Sentinel-2"
      x = "operational"
    else
      x = "pre-ops"
    end
  end
  if x == "operational"
    x = (options[:gSD_api])[:dhus]
  end
  if x == "pre-ops"
    x = (options[:gSD_api])[:s3]
    user = "s3guest"
    pw = "s3guest"
  end
  return [user, pw, x]
end



#' get ERS API key from user input
#'
#' @param username username
#' @param password password
#' @keywords internal
#' @noRd
function ERS_login(username, password)
  x = POST( join([ (options[:gSD_api])[:ee], "login?jsonRequest={\"username\":\"", username, "\",\"password\":\"", password, "\",\"authType\":\"EROS\",\"catalogId\":\"EE\"}" ]))
  stop_for_status(x, "connect to server.")
  warn_for_status(x)
  (content(x))[:data]
end



#' logout from ERS with API key
#'
#' @param api.key api.key
#' @keywords internal
#' @noRd
function ERS_logout(api_key)
  x = gSD_get(join([ (options[:gSD_api])[:ee], "logout?jsonRequest={\"apiKey\":\"", api_key, "\"}" ]))
  stop_for_status(x, "connect to server.")
  warn_for_status(x)
  (content(x))[:data]
end



#' get EE products
#'
#' @param api.key api.key
#' @param wildcard wildcard
#' @keywords internal
#' @noRd
function EE_ds(api_key; wildcard = nothing)
  q = join([ (options[:gSD_api])[:ee], "datasets?jsonRequest={\"apiKey\":\"", api_key, "\"}" ]) #, if(isnothing(wildcard)) "}" else  ",\"datasetName\":\"", wildcard, "\"}")
  if !isnothing(wildcard) 
    q = gsub("}", join([ ",\"datasetName\":\"", wildcard, "\"}" ]), q)
  end
  x = gSD_get(q)
  getindex.( content(x)[:data], :datasetName)
end



#' query EE
#'
#' @param aoi aoi
#' @param time_range time_range
#' @param name name
#' @param api.key api.key
#' @param meta.fields meta.fields
#'
#' @importFrom sf st_bbox st_as_text
#' @importFrom xml2 as_list
#'
#' @keywords internal
#' @noRd
function EE_query(aoi, time_range, name, api_key; meta_fields = nothing )

  spatialFilter = join([ "\"spatialFilter\":{\"filterType\":\"mbr\",\"lowerLeft\":{\"latitude\":", st_bbox(aoi)[:ymin], ",\"longitude\":", st_bbox(aoi)[:xmin], "},\"upperRight\":{\"latitude\":", st_bbox(aoi)[:ymax], ",\"longitude\":", st_bbox(aoi)[:xmin], "}}" ])
  temporalFilter = join([ "\"temporalFilter\":{\"startDate\":\"", time_range[1], "\",\"endDate\":\"", time_range[2], "\"}" ])
  
  out("Searching USGS EarthExplorer for available products...")
  query = map( (x, ak = api_key, sf = spatialFilter, tf = temporalFilter) -> gSD_get(join([ options[:gSD_api][:ee], "search?jsonRequest={\"apiKey\":\"",ak,"\",\"datasetName\":\"",x,"\",",sf,',',tf,",\"startingNumber\":1,\"sortOrder\":\"ASC\",\"maxResults\":50000}" ])), name)
  query_cont = content.(query)
  if length(name) == 1 
    if query_cont[1][:error] != ""
      out("Invalid query. This dataset seems to be not available for the specified time range.", type = 3)
    end
  end
  query_use = map( x -> (x[:error] == "") && (length(x[:data][:results]) != 0), query_cont )
  query_cont = query_cont[query_use]
  query_names = name[query_use]

  query_results = map( x -> x[:data][:results], query_cont )
  if length(query.results) != 0
    #query.df <- unlist(mapply(y = query.results, n = query.names, function(y, n) lapply(y, function(x, ds_name = n){
    #  x.names <- names(x)
    #  x.char <- as.character(x)
    #  
    #  # Make sf polygon filed from spatialFootprint
    #  spf.sub <- grep("spatialFoot", x.names)
    #  spf <- unlist(x[spf.sub])
    #  spf <- as.numeric(spf[grep("coordinates", names(spf))])
    #  spf.sf <- .make_aoi(cbind(spf[seq(1, length(spf), by = 2)], spf[seq(2, length(spf), by = 2)]), type = "sf", quiet = T)
    #  
    #  df <- rbind.data.frame(x.char, stringsAsFactors = F)
    #  colnames(df) <- x.names
    #  df[,spf.sub] <- st_as_text(spf.sf)
    #  df <- cbind.data.frame(df, ds_name, stringsAsFactors = F)
    #  colnames(df)[ncol(df)] <- "product"
    #  return(df)
    #}), SIMPLIFY = F), recursive = F)
    query_df = map( (y = query_results, n = query_names) -> map( function (x, ds_name = n)
                                                                   x_names = names(x)
                                                                   x_char = string.(x)
    
                                                                   # Make sf polygon filed from spatialFootprint
                                                                   spf_sub = findall(occursin.("spatialFoot", x.names))
                                                                   spf = x[spf_sub]
                                                                   spf = tryparse( Int64, spf[ findall(occursin.("coordinates", names(spf))) ] )
                                                                   spf_sf = make_aoi( hcat( spf[range(1, length(spf), by = 2)], spf[range(2, length(spf), by = 2)] ), type = "sf", quiet = true )

                                                                   df = DataFrame( x_char )
                                                                   rename!( df, x_names )
                                                                   df[spf_sub] = st_as_text(spf_sf)
                                                                   df = DataFrame( df, ds_name )

                                                                   rename!( df[ncol(df)], "product" )
                                                                   return df
                                                                 end, y ), [query_results, query_names] )

    
    ## Read out meta data
    out("Reading meta data of search results from USGS EarthExplorer...")
    meta = gSD_get.( getindex.( query_df, :metadataUrl ) )
    meta_list = map( x -> as_list( xml_contents( xml_contents( content(x) )[1] ) ), meta )
    meta_val = map.( function (x)
                       try
                         z = x[:metadataValue][1]
                       catch exception
                         return nothing
                       end
                       return z
                     end, meta_list )
    meta_name = map.( x -> names(x)[:name], meta_list )


    ## Define meta fields that are usefull for the query output
    if isnothing(meta_fields)
      meta_fields = unique(meta_name)
    end
    meta_subs = map( (mnames, mf = meta_fields) -> map( (x, mn = mnames) -> findall( x .== mn ), mf ), meta_name )
    meta_df = map.( function (v = meta_val, n = meta_name, i = meta_subs )
                      x = v[i]
                      x = map( y -> isnothing(y) ? "" : y, x )
                      x = DataFrame(x)
                      rename!( x, gsub(" ", "", n[i]) )
                      return x
                    end, [meta_val, meta_name, meta_subs] )

    query_df = map( function (q = query_df, m = meta_df)
                      ## apply meaningful order and replace startTime and endTime with meta outputs
                      x = DataFrame( q[:acquisitionDate], m, q[:, 4:end] )
                      reaname!( x[1], names(q[1]) )
                      return x
                    end, [query_df, meta_df] )
    

    return_names = unique(names.(query_df))
    return_df = DataFrame( names!( zeroes( Int64, length(return_names) ) ), return_names )
    return_df = DataFrame( map( function (x, rn = return_names, rdf = return_df)
                                  rdf[1, findall( in(rn), names(x) )] = x
                                  return rdf
                                end, query_df ) )
    return return_df
  else
    return nothing
  end
end



#' preview EE record
#'
#' @param record record
#' @param on_map on_map
#' @param show_aoi show_aoi
#' @param verbose verbose
#'
#' @importFrom getPass getPass
#' @importFrom httr GET write_disk authenticate
#' @importFrom raster stack plotRGB crs crs<- extent extent<- NAvalue
#' @importFrom sf st_as_sfc st_crs as_Spatial
#' @importFrom mapview viewRGB addFeatures
#'
#' @keywords internal
#' @noRd
function EE_preview(record; on_map = true, show_aoi = true, verbose = true)

  if verbose <: Bool
    options[:verbose] = verbose
  end
  
  ## Intercept false inputs and get inputs
  url_icon = record[:browseUrl]
  if isnan(url_icon)
    out("Argument 'record' is invalid or no preview is available.", type = 3)
  end
  if length(url.icon) > 1
    out("Argument 'record' must contain only a single record, represented by a single row data.frame.")
  end
  char_args = Dict( url_icon => url_icon )
  for x in values(char_args) 
    if !all(isa.( x, AbstractString ))
      out( join([ "Argument '", names(x), "' needs to be of type 'character'." ]), type = 3)
    end
  end

  if length(findall("https", url_icon)) == 0
    out("No preview available for this record or product.", msg = true)
  else
    ## Recieve preview
    file_dir = join(tempname(),".jpg")
    gSD_get(url_icon, dir_file = file_dir)
    preview = stack(file_dir)
    #NAvalue(preview) <- 0

    if isTrue(on_map)

      ## create footprint
      footprint = st_as_sfc( Vector(record[:spatialFootprint]) )
      st_crs(footprint) = 4326
      footprint = as_Spatial(footprint)
       
      ## create preview
      crs(preview) = crs(footprint)
      extent(preview) = extent(footprint)
      #preview <- aggregate(preview, 2) # make it faster
       
      ## create map
      map = suppressWarnings(viewRGB(preview, r = 1, g = 2, b = 3) )

      if isTRUE(show_aoi)
        if isFALSE( options[:aoi_set] )
          out("Preview without AOI, since no AOI has been set yet (use 'set_aoi()' to define an AOI).", type = 2)
        else
          aoi_sf = options[:aoi]
          #aoi_sf = make_aoi(aoi_m, type = "sf", quiet = true )
          map = addFeatures(map, aoi_sf)
        end
      end
      map
    else
      plotRGB(preview)
    end
  end
end



#' convert MODIS product names
#'
#' @param names names
#' @keywords internal
#' @noRd
function convMODIS_names(names)
  map( function (x)
         y = split( x, "_", keepempty = false )[1]
         y = y[2:end]
         if length(y) > 1
           y = join( y[1:end-1], "_" )
         end
         return y 
       end, names )
end



#' USGS ESPA ordering functon
#'
#' @param id id
#' @param level level
#' @param username username
#' @param password password
#' @param format format
#' @keywords internal
#' @importFrom httr content
#' @noRd
function ESPA_order( id, username, password, verbose; level = "sr", format = "gtiff" )
  
  ## check query and abort, if not available
  out("Ordering requested items from ESPA...")
  checked = map( function (x, v = verbose)
                   r = gSD_get( join([ options[:gSD_api][:espa], "available-products/", x ]), options[:gSD_usgs_user], options[:gSD_usgs_pass] )
                   if names(content(r)) == "not_implemented"
                     out( join([ "'", x, "': This ID is invalid, as it cannot be found in the ESPA database. Please remove it from input and reexecute." ]), type = 3)
                   end
                   vcat(x, r)
                 end, id )

  ## group request by collection (single or multi order)
  req_data = map( x -> vcat( names( content(x[2]) ), x[1] ), checked )
  coll = map( x -> x[1][1], req.data  )
  coll_uni = unique(coll)
  out(join([ "Collecting from ", length(coll_uni), " collection(s) [", join(coll_uni, ", "), "], resulting in ", length(coll_uni), " order(s)..." ]))
  req_coll = map( (x, c = coll, rd = req_data) -> rd[findall( c .== x )], coll_uni )

  ## build request
  req_body = map( function (x, p = level, f = format)
                    i = join( map( y -> y[2], x ), "\", \"" )
                    join([ "{\"", x[1][1], "\": { \"inputs\": [\"", i, "\"], \"products\": [\"", p, "\"]}, \"format\": \"", f, "\"}" ])
                  end, req_coll )
  
  ## order
  order = map( (x, user = username, pass = password) -> gSD_post(url = join([ options[:gSD_api][:espa], "order/" ]), username = user, password = pass, body = x), req_body )
  order_list = map( x -> content(x)[1], order )
  out(join([ "Products '", join( id, "', '" ), "' have been ordered successfully:" ]))
  out(join([ "[level = '", level, "', format = '", format, "', order ID(s) '", join( order_list, "', '" ), "']." ]))
  return order_list
end



#' USGS ESPA downloading functon
#'
#' @param order.list order.list
#' @param username username
#' @param password password
#' @param file.down file.down
#' @param delay delay
#'
#' @importFrom utils head tail
#'
#' @keywords internal
#' @noRd
## check order(s)
function ESPA_download(order_list, username, password, file_down, dir_out; delay = 10)
  remain_active = true
  ini = true
  show_status = true
  while remain_active
    
    ## get tiems
    items = map( (x, user = username, pass = password) -> content(gSD_get(join([ options[:gSD_api][:espa], "item-status/", x ]), user, pass)),
                 order_list )

    ## get items content
    items = map( x -> map( function (y)
                             rename!( r, names(y) )
                             return r
                           end, x[1] ), items )


    ## make items data.frame containing recieve status
    items = DataFrame( reduce( hcat, reduce.( hcat, items ) ) )
    names_required = map( x -> join( split( split( x, "/" )[end], "_" )[1], collapse = "_" ), file_down )
    items = items[ map( (x, y = items[:name]) -> findall( y .== x ), names_required ), : ]
    items = vcat( items, items[:status] == "complete" )
    items = DataFrame( items, file_down ) #sapply(as.character(items$name), function(x, l = level) paste0(dir_out, "/", x, "_", toupper(level), ".tar.gz"), USE.NAMES = F), stringsAsFactors = F)
    rename!( items[end-1:end], [ Symbol("available"), Symbol("file") ] )

    if ini
      items_df = DataFrame( items, falses(length(items[:status])) )
      rename!( items_df[end], Symbol("recieved") )
      ini = false
    else
      items_df = DataFrame( items, items_df[:recieved] )
    end
    if isTRUE(force)
      emp = map( x -> isfile(x) && rm(x), items_df[:file] )
    end
    items_df[:recieved] = map( isfile, items_df[:file] )

    ## Items to download
    if all(items_df[:available]) && all(items_df[:recieved])
      remain_active = false
    else
  
      ## Download or wait for status
      sub_download = intersect( findall( items_df[:available] .== true ), findall( items_df[:recieved] .== false ) )
      if length(sub_download) > 0
        
        items_get = items_df[sub_download,:]
        out( join([ "Starting download of product(s) '", join( items_get[:name], "', " ), "'." ]), msg = true )
        items_df[Symbol(recieved[sub_download])] = map( function (x, d = dir_out)
                                                          y = DataFrame(x)
                                                          rename!( y, names(x) )
                                                          gSD_download( name=y[:name], url_file=y[:product_dload_url], url_checksum=y[:cksum_download_url], file=y[:file] )
                                                        end, items_get )
        show_status = true
      else
        if isTRUE(show_status)
          out(join([ "Waiting for product(s) '", join( items_df[Symbol(name[findall(items_df[:available] .== false)])], "', " ), "' to be ready for download from ESPA (this may take a while)..." ]))
          out("Note: It is also possible to terminate the function and call it again later by providing the displayed order ID(s) as input to 'espa_order'.", msg = true)
        end
        show_status = false
      end
    end
    sleep(delay) #wait before reconnecting to ESPA to recheck status
  end
end



#' make aoi
#'
#' @param aoi aoi
#' @keywords internal
#' @importFrom sp SpatialPolygons
#' @importFrom sf st_sfc st_polygon st_crs st_as_sf st_coordinates st_transform st_crs<- as_Spatial
#' @noRd
function make_aoi(aoi; type = "matrix", quiet = false)
  ## if not sfc, convert to sfc
  if !inherits(aoi, ["Spatial", "sfc", "matrix"])
    out("Argument 'aoi' needs to be a 'SpatialPolygons' or 'sfc_POLYGON' or 'matrix' object.", type = 3)
  end
  if inherits(aoi, "matrix")
    if !all( aoi[1,:] .== aoi[length(aoi[:,1]),:])
      aoi = hcat( aoi, aoi[1,:] )
    end
    aoi = st_sfc( st_polygon( Vector(aoi) ), crs = 4326 )
    if isFALSE(quiet)
      out(join([ "Argument 'aoi' is a matrix, assuming '", st_crs(aoi)[:proj4string], "' projection." ]), type = 2)
    end
  end
  if inherits(aoi, "Spatial") 
    aoi <- st_as_sf(aoi)
  end

  ## check projection
  if isnan( st_crs(aoi) )
    st_crs(aoi) <- 4326
    if isFALSE(quiet)
      out(join([ "Argument 'aoi' has no projection, assuming '", st_crs(aoi)[:proj4string], "' projection." ]), type = 2)
    end
  end
  if length( grep("WGS84", grep("longlat", st_crs(aoi)[:proj4string], value = true), value = true) ) != 1
    aoi = st_transform(aoi, 4326)
  end

    
  ## get coordinates
  aoi_m = st_coordinates(aoi)[:,[1,2]]
  aoi_sf = st_sfc( st_polygon(aoi_m), crs = 4326 )
  aoi_sp = as_Spatial(aoi_sf)
  
  if type == "matrix" 
    return aoi_m
  end
  if type == "sf"
    return aoi_sf
  end
  if type == "sp" 
    return aoi_sp
  end
end





options = Dict(
  :gSD_api => Dict(
    :dhus => "https://scihub.copernicus.eu/dhus/",
    :s3 => "https://scihub.copernicus.eu/s3/",
    :espa => "https://espa.cr.usgs.gov/api/v0/",
    :ee => "https://earthexplorer.usgs.gov/inventory/json/v/1.4.0/",
    :aws_l8 => "https://landsat-pds.s3.amazonaws.com/c1/L8/",
    :aws_l8_sl => "https://landsat-pds.s3.amazonaws.com/c1/L8/scene_list.gz"
  ),
  :gSD_verbose => false,
  :gSD_cophub_user => false,
  :gSD_cophub_pass => false,
  :gSD_cophub_set => false,
  :gSD_usgs_user => false,
  :gSD_usgs_pass => false,
  :gSD_usgs_set => false,
  :gSD_usgs_apikey => false,
  :gSD_archive => false,
  :gSD_archive_set => false,
  :gSD_aoi => false,
  :gSD_aoi_set => false
)
#toset <- !(names(op.gSD) %in% names(op))
#if(any(toset)) options(op.gSD[toset])

end # module