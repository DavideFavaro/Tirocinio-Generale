module GroundDataV
"""
Module for the download and processing of atmospheric data gathered by measuring stations located in Veneto, Italy
"""



using CombinedParsers
using CombinedParsers.Regexp
using CSV
using DataFrames
using Dates
using HTTP
using JSON
using JSONTables
using Revise


export getData,
       getRegionAttributes, getRegionIds, getRegionStationInfo


@syntax informations = Sequence(
                            :info => Sequence(
                                           Sequence( "<", "PERIODO", ">" ),
                                           "<![CDATA[",
                                           re"[^\[\]]+",
                                           "]]>",
                                           "</PERIODO>"
                                       ),
                            :begin => Sequence(
                                          Sequence( "<", "INIZIO", ">" ),
                                          Numeric(Int64),
                                          "</INIZIO>"
                                      ),
                            :end => Sequence(
                                Sequence( "<", "FINE", ">" ),
                                Numeric(Int64),
                                "</FINE>"
                                ),
                            :proj => Sequence(
                                Sequence( "<", "PROJECTION", ">" ),
                                re"[^<:>]+",
                                ":",
                                Numeric(Int64),
                                "</PROJECTION>"
                                )
                       )

@syntax stations = Repeat( 
                        Sequence(
                            re"<[^<>]+>",
                            Either(
                                Numeric(Int64),
                                Numeric(Float64),
                                Sequence( "<![CDATA[", re"[^\[\]]+", "]]>" ),
                                re"[^<>]+"
                            ),
                            re"</[^<>]+>"
                        )
                    )

@syntax sensor = Sequence(
                    Repeat(
                        re"<(?!D)[^>]+>",
                        Either(
                            Numeric(Int64),
                            Sequence( "<![CDATA[", re"[^\[\]]+", "]]>" ),
                            re"[^>]+"
                        ),
                        re"</[^>]+>"
                    ),
                    Repeat(
                        Sequence( "<DATI ISTANTE=\"", Numeric(Int64), "\">" ),
                        Sequence( "<VM>", Numeric(Float64), "</VM>" ),
                        "</DATI>"
                    ) 
                  )



"""
    getRegionAttributes( [ type::Symbol=:METEO ] )

Obtain the names of the columns of the region's dataframe required by `GroundData.createMap`'s `attributes` parameter to create
`GroundData.standardize`'s `map` parameter
"""
function getRegionAttributes( type::Symbol=:METEO )
    return type == :METEO ?
               [ :paramnm, :unitnm, :value, nothing, :instant, :x, :y, :quota, nothing, nothing, :rmh ] :
               type == :AIRQUALITY ?
                   [ :param, nothing, :value, nothing, :date, :lon, :lat, nothing, nothing, :tipozona ] :
                   throw( DomainError( type, "`type` must be either `:METEO` OR `:AIRQUALITY`" ) )
end



"""
    getRegionAttributes( [ type::Symbol=:METEO ] )

Obtain the names of the columns of the dataframe required for `GroundData.standardize`'s `bridge` parameter
"""
function getRegionIds( type::Symbol=:METEO )
    return type == :METEO ? :idstaz :
               type == :AIRQUALITY ? :codseqst : throw( DomainError( type, "`type` must be either `:METEO` OR `:AIRQUALITY`" ) )
end



"""
    getRegionStationInfo( [ type::Symbol=:METEO  ] )

Obtain the names of the columns of the region's stations dataframe required by `GroundData.createMap`'s `attributes` parameter to be used
in `GroundData.generateUuidsTable`
"""
function getRegionStationsInfo( type::Symbol=:METEO )
    return type == :METEO ?
               [ :idstaz, :nome, :x, :y ] :
               type == :AIRQUALITY ?
                   [ :codseqst, :nome, :lon, :lat ] :
                   throw( DomainError( type, "`type` must be either `:METEO` OR `:AIRQUALITY`" ) )
end



"""
    getStationsInfo()::DataFrames.DataFrame

Obtain the dataframe containing informations on all available ARPAV stations 
"""
function getMeteoStationsData()
#Get the information on the available stations
    # Get the page containing informations on all the available stations
    page = String( HTTP.get( "https://www.arpa.veneto.it/bollettini/meteo/h24/img07/stazioni.xml" ).body )
    # Keep only the useful portion of the page
    useful = split( page, "</PERIODO>" )[2][1:end-15]

    # Split the downloaded string in substrings one for each station
    arr = split( useful, r"</?STAZIONE>", keepempty=false )

    #info = informations(arr[1])
    stats = stations.(arr[2:end])

    #Create the dictionary of the stations, checking the type of `attribute[2]` (contains the value of the attribute)
    dict = [ Dict(
                Symbol( lowercase( String( attribute[1][2] ) ) )
                =>
                isa(attribute[2], Number) ? # If it's a number leave it as is
                   attribute[2] :
                   isa(attribute[2], Array) ? #If it is an array (of strings) convert it to String 
                       lowercase( String( attribute[2] ) ) : # If it is a Tuple ( Es. `("<![CDATA[", "Arabba", "]]>")` ) take the actual value ( `"Arabba"` ) 
                       titlecase( String( attribute[2][2] ) )
                for attribute in station
             ) for station in stats ]

    df = DataFrame(dict)
    insertcols!( df, :rmh => "0m" )
    return df
end

# res = getMeteoStationData()



"""
    getStationsData( stats::AbstractVector{String} )

Obtain the dataframe containing the data of all the stations described by the elements of `stats`
"""
function getMeteoData( ids::AbstractVector{Int64} )
    pages_vect = [ String( HTTP.get("https://www.arpa.veneto.it/bollettini/meteo/h24/img07/$(lpad(id, 4, "0")).xml").body ) for id in ids ]

    # For each of the station remove the first part of the string and the closing tags
    stat_strings_vect = [ split(page, "</ATTIVAZIONE>")[2][1:end-26] for page in pages_vect ]

    # From each station generate the corresponding array of sensors
    sensor_vect = split.( stat_strings_vect, r"</?SENSORE>", keepempty=false )

    # vect[i][j] i-th station, j-th sensor
    # vect[i][j][1] j-th sensor's info (Vector)
    # vect[i][j][2] j-th sensor's data (Vector)
    # vect[i][j][2][l] j-th sensor'data, l-th measurement (Vector)
    # vect[i][j][2][l][4] l-th measurement's values (Vector 1 => mean, 2 => min, 3 => max)
    vect = [ sensor.(sensor_group) for sensor_group in sensor_vect ]

    data_dict_vect = [
                        push!(
                            Dict(
                                Symbol( lowercase( String( attribute[3] ) ) )
                                =>
                                isa( attribute[5], Number ) ? attribute[5] :
                                    isa( attribute[5], Array ) ? String( attribute[5] ) : titlecase( String( attribute[5][2] ) )
                              for attribute in sensor[1]
                            ),
                            :idstaz => id,
                            :instant => entry[2],
                            :value => entry[5]
                        )
                        for (id, station) in zip(ids, vect)
                        for sensor in station   
                        for entry in sensor[2]
                     ]   
    df = DataFrame( data_dict_vect )
    transform!(
        df,
        [:unitnm] => ByRow( x -> x = replace( replace( x, "\xb0" => "°" ), "2" => "²" ) ) => :unitnm,
        [:instant] => ByRow( x -> DateTime( string(x), "yyyymmddHHMM" ) ) => :instant
    )
    return df
end

#   df = getMeteoStationsData()
#   data = getMeteoSensorsInfo( df[1:3, :idstaz] )



"""
    getAqStationsData()

Obtain a dataframe containing informations on the measuring stations in the area
"""
function getAqStationsData()
    page = String( HTTP.get("http://213.217.132.81/aria-json/exported/aria/stats.json").body )
    jst = jsontable( page[14:end] )
    return DataFrame( jst )
end

#   res = getAqStationsData()



"""
    getAqData()

Obtain the data gathered by the stations in the area
"""
function getAqData()
    page = String( HTTP.get("http://213.217.132.81/aria-json/exported/aria/data.json").body )

    # js["stazioni"][i]                                            i-th station
    # js["stazioni"][i]["codseqst"]                                i-th station's id
    # js["stazioni"][i]["misurazioni"][j]                          i-th station's j-th measurement
    # js["stazioni"][i]["misurazioni"][j]["pm10"/"ozono"][l]       j-th measurement's l-th entry 
    js = JSON.parse( page )["stazioni"]

    arr = []
    for station in js
        if !isempty( station["misurazioni"] )
            for measurement in station["misurazioni"]
                for key in collect( keys( measurement ) )
                    for entry in measurement[key]
                        dict = Dict(
                            :codseqst => station["codseqst"],
                            :param => key,
                            :date => entry["data"],
                            :value => entry["mis"]
                        )
                        push!( arr, dict )
                    end
                end
            end
        else
            push!( arr, Dict( :codseqst => station["codseqst"], :param => missing, :date => missing, :value => missing ) )
        end
    end
    return DataFrame(arr)
end

#   df = getAqSensorsData()



"""

Obtain `type` data from `source`
"""
function getData(; type::Symbol=:METEO, source::Symbol=:STATIONS )
    if type == :METEO
        stations = getMeteoStationsData()
        if source == :STATIONS
            return stations
        else
            return getMeteoData( stations[:, :idstaz] )
        end
    else
        if source == :STATIONS
            return getAqStationsData()
        else
            return getAqData()
        end
    end
end

#   ressta = getData()
#   ressen = getData( source=:SENSORS )
#   ressta = getData( type=:AIRQUALITY, source=:STATIONS )
#   ressen = getData( type=:AIRQUALITY, source=:SENSORS )



end # moduled