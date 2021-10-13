module GroundData


using HTTP

using CombinedParsers
using CombinedParsers.Regexp

using DataFrames
using Dates
using CSV

using Revise



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
                                          Numeric(Int),
                                          "</INIZIO>"
                                      ),
                            :end => Sequence(
                                Sequence( "<", "FINE", ">" ),
                                Numeric(Int),
                                "</FINE>"
                                ),
                            :proj => Sequence(
                                Sequence( "<", "PROJECTION", ">" ),
                                re"[^<:>]+",
                                ":",
                                Numeric(Int),
                                "</PROJECTION>"
                                )
                       )


@syntax stations = Repeat( 
                        Sequence(
                            re"<[^<>]+>",
                            Either(
                                Numeric(Int),
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
                            Numeric(Int),
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
    getStationsInfo()::DataFrames.DataFrame

Obtain the dataframe containing informations on all available ARPAV stations 
"""
function getStationsInfo()

#Get the information on the available stations
    # Get the page containing informations on all the available stations
    page = HTTP.get( "https://www.arpa.veneto.it/bollettini/meteo/h24/img07/stazioni.xml" )
    str_page = String(page.body)
    # Keep only the useful portion of the page
    useful = split( str_page, "</PERIODO>" )[2][1:end-15]

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

    return df
end






"""
    getStationsData( stats::AbstractVector{String} )

Obtain the dataframe containing the data of all the stations described by the elements of `stats`
"""
# INCOMPLETO
function getSensorsInfo( stats::AbstractVector{Int64} )
    pages_vect = [ String( HTTP.get("https://www.arpa.veneto.it/bollettini/meteo/h24/img07/$(lpad(stat, 4, "0")).xml").body ) for stat in stats ]

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

    sensor_dict_vect = []
    data_dict_vect = []
    for (id, station) in zip(stats, vect)
        for sensor in station
            sensor_dict = setindex!(
                Dict(
                    Symbol( lowercase( String( attribute[3] ) ) )
                    =>
                    isa( attribute[5], Number ) ? attribute[5] :
                        isa( attribute[5], Array ) ? String( attribute[5] ) : titlecase( String( attribute[5][2] ) )
                    for attribute in sensor[1]
                ),
                id,
                :station_id
            )
            
            for entry in sensor[2]
                data_dict = Dict(
                    :station_id => id,
                    :sensor_id => sensor_dict[:id],
                    :instant => entry[2],
                    :value => entry[5]
                )
                push!( data_dict_vect, data_dict )
            end

            push!( sensor_dict_vect, sensor_dict )   
        end
    end

    dfs = ( DataFrame( sensor_dict_vect ), DataFrame( data_dict_vect ) )

    return dfs
end




df = getStationsInfo()
sen_dict, data_dict = getSensorsInfo( df[1:3, :idstaz] )






end # moduled