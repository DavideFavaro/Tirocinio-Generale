module GetSpatialDataDev

#Global variables
#= imports
  library(httr)
  library(xml2)
  source("getSpatialData_dev_internal.R")
=#

using DataFrames
using Dates

include("getSpatialDataDevInternal.jl")

export cophub_api, getSentinel_query2, login_CopHub, set_archive #Functions


function cophub_api(x::AbstractString, p::AbstractString, user, pw) 
  if x == "auto"
    x = ( p == "Sentinel-1" || p == "Sentinel-2" ) ? "operational" : "pre-ops"
  elseif x == "operational"
    x = options[:gSD_api][:dhus]
  elseif x == "pre-ops"
    x = options[:gSD_api][:s3]
    user = "s3guest"
    pw = "s3guest"
  end
  return [ user, pw, x ]
end




function set_archive(dir_archive) 
  if !isa(dir_archive, AbstractString)
    out( "Argument 'dir_archive' needs to be of type 'String'.", type = 3 )
  end
  if !isdir(dir_archive)
    out("The defined directory does not exist.", type = 3) 
  end
  options[:gSD_archive] = dir_archive
  options[:gSD_archive_set] = true
end




function login_CopHub(username; password = nothing) 
  if isnothing(password)
    password = getPass()
  end
  char_args = Dict( :username => username, :password => password )
  for (key, val) in char_args
    if !isa(val, AbstractString)
      out( "Argument '$key' needs to be of type 'String'.", type = 3 )
    end
  end
  options[:gSD_cophub_user] = username
  options[:gSD_cophub_pass] = password
  options[:gSD_cophub_set] = true
end




function getSentinel_query2(time_range, platform; aoi=nothing, username=nothing, password=nothing, hub="auto", verbose=true, ingestiondate=false, extra=nothing)
  if options[:gSD_cophub_set]
    if isnothing(username) 
      username = options[:gSD_cophub_user]
    end
    if isnothing(password) 
      password = options[:gSD_cophub_pass]
    end
  end

  if !isa(username, AbstractString) 
    println("Argument 'username' needs to be of type 'String'. You can use 'login_CopHub()' to define your login credentials globally.")
  end

  if isnothing(password)
    password = getPass()
  end

    if verbose <: Bool
      options[:verbose] = verbose
    #if isnothing(aoi)
    #  if options[:gSD_aoi_set]
    #    aoi = options[:gSD_aoi]
    #  else
    #    out("Argument 'aoi' is undefined and no session AOI could be obtained. Define aoi or use set_aoi() to define a session AOI.", type = 3)
    #  end
    #end
    #aoi = make_aoi(aoi, type = "matrix")
  end


  char_args = Dict( :time_range => time_range, :platform => platform )
  for (key, val) in char_args
    if !isa(val, AbstractString)
      out( "Argument '$key' needs to be of type 'character'.", type = 3)
    end
  end
  if length(time_range) != 2
    out( "Argument 'time_range' must contain two elements (start and stop time).", type = 3 )
  end


  function cop_url( ext_xy, urlRoot, platform, timeRange, rowStart )
    if !ingestiondate
      qs = Dict( :urlRoot => join([urlRoot, "/"]),
                 :search => ( "search?start=", "&rows=100&q=(" ), 
                 :and => "%20AND%20", 
                 :aoi_poly => ( "footprint:%22Intersects(POLYGON((", ")))%22" ), 
                 :platformName => "platformname:", 
                 :time => Dict( "[" => "beginposition:%5b", "to" => "%20TO%20", "]" => "%5d" ) )
    else
      qs = Dict( :urlRoot => join([urlRoot, "/"]),
                 :search => ( "search?start=", "&rows=100&q=(" ), 
                 :and => "%20AND%20", 
                 :aoi_poly => ( "footprint:%22Intersects(POLYGON((", ")))%22" ), 
                 :platformName => "platformname:", 
                 :time => Dict( "[" => "ingestiondate:%5b", "to" => "%20TO%20", "]" => "%5d" ) )
    end
  

    #aoi_str = join( join.( ext_xy, "%20" ), "," )

    timeRange = map( x -> tt = format( DateTime( timeRange ), "%Y-%m-%dT%H:%M:%S.000Z" ), timeRange )
    if isnothing(extra)
      urlstring = join([ qs[:urlRoot], (qs[:search])[1], rowStart,
                         (qs[:search])[2], qs[:platformName], platform,
                         qs[:and], (qs[:time])["["], timeRange[1],
                         (qs[:time])["to"], timeRange[2], (qs[:time])["]"], ")" ])
    else
      urlstring = join([ qs[:urlRoot], (qs[:search])[1], rowStart,
                         (qs[:search])[2], qs[:platformName], platform,
                         qs[:and], extra, qs[:and], (qs[:time])["["],
                         timeRange[1], (qs[:time])["to"], timeRange[2],
                         (qs[:time])["]"], ")" ])  
    end
    println(urlstring)
    println( now() )
    return urlstring
  end

  cred = cophub_api(hub, platform, username, password)
  row_start = -100
  re_query = true
  give_return = true
  query_list = []
  while re_query
    row_start += 100
    query = gSD_get( cop_url( aoi, cred[3], platform, time_range, row_start ), cred[1], cred[2] )
    query_xml = xml_contents(contents(query))
    append!( query_list, map( x -> xml_contents(x), filter(str -> contains("entry", str), query_xml) ) )


    if length(query_list) == 0 && row_start == 0
      out( "No results could be obtained for this request.", msg = true )
      re_query = false
      give_return = false
    end
    if length(query_list) != row_start + 100
      re_query = false
    end
  end






  if give_return
    field_tag = [ "title", "link href", "link rel=\"alternative\"", "link rel=\"icon\"", "summary", "name" ]
    field_names = [ "title", "url", "url.alt", "url.icon", "summary" ]

    query_cont = map( (x, field = field_tag) -> map( (f, y = x) -> filter( line ->  occursin(f, line), y ), field ), query_list )
    query_names = map( (x, field_n = field_names) -> append!( field_n, map( y -> split( y, '\"', keepempty = false )[2], x[length(field_n)+1 : end] ) ), query_cont )
    query_fields = map( x -> map( y -> split( split(y, '<', keepempty = false)[2], '>', keepempty = false )[1], x ), query_cont )
    
    function fun( i, qf, qc, qn )
        x = qf[i]
        x[ ismissing.(x) ] = map( y -> split( split(y, "href=\"")[2], "\">" )[1], qc[i][ ismissing.(x) .== true ] )
        #names(x) <- qn[[i]]
        return x
    end
    query_fields = map( (i, qf = query_fields, qc = query_cont, qn = query_names ) -> fun( i, qf, qc, qn ), 1:length(query_fields) )

    return_names = unique(query_names)
    return_df = rename!( convert( DataFrame, return_names  ), zeros( Int64,  length(return_names) ) )
    return_df = combine( DataFrame( map( (x, rn = return_names, rdf = return_df) -> rdf[1, findall( in(names(x)), rn )] = map( y -> string(y), x ) ) ) )
  end

  if give_return
    return return_df
  end
end

end # module