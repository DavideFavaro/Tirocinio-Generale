module GroundDataAA
"""
Module for the download and processing of atmospheric data gathered by measuring stations located in Alto Adige, Italy
"""

#=
Alto Adige:
    Link:
        Stazioni Meteo (JSON):
            http://dati.retecivica.bz.it/services/meteo/v1/stations
        Dati Sensori Meteo (JSON):
            http://dati.retecivica.bz.it/services/meteo/v1/sensors
        Stazioni QA (JSON):
            http://dati.retecivica.bz.it/services/airquality/stations
        Dati sensori QA (JSON):
            http://dati.retecivica.bz.it/services/airquality/timeseries
=#



using CSV
using DataFrames
using HTTP
using JSONTables


export getData,
       attributes, ids


const attributes = Dict(
                      :METEO      => [ :DESC_I, :UNIT, :VALUE, nothing, :DATE, :LONG, :LAT, :ALT, :ALT, nothing, nothing ],
                      :AIRQUALITY => [ :MCODE, nothing, :VALUE, nothing, :DATE, :LONG, :LAT, nothing, :FLAGS, nothing ]
                   )

const ids = Dict(
                :METEO      => :SCODE,
                :AIRQUALITY => :SCODE
            )

const stat_info = Dict(
                          :METEO      => [ :SCODE, :NAME_I, :LONG, :LAT ],
                          :AIRQUALITY => [ :SCODE, :NAME_I, :LONG, :LAT ]
                      )



"""
    getData(; <keyword arguments> )

Obtain data of category `type` from `source`

# Arguments
 - `type::Symbol=:METEO`: defines the type of data to be downloaded may either be `:METEO` or `:AIRQUALITY`
 - `source::Symbol=:STATIONS`: defines if the data to be downloaded has to regard information on the stations or their actual measurements, as such may either be `:STATIONS` or `:SENSORS`
"""
function getData(; type::Symbol=:METEO, source::Symbol=:STATIONS )
    opt1 = type == :METEO ? "meteo/v1" : "airquality"
    opt2 = source == :STATIONS ? "stations" : type == :METEO ? "sensors" : "timeseries"

    page = String( HTTP.get( "http://dati.retecivica.bz.it/services/$opt1/$opt2" ).body )

    if source == :STATIONS
        if type == :METEO
            chars = " : "
            div = "\r\n\t\t},\r\n\t\t"
            lim = 11
        else
            chars = ":"
            div = "}\r\n,"
            lim = 5
        end
        features = split( page , "\"features\"$chars"  )[2]
        stations = split( features, "\"properties\"$chars" )
        stations = [ split( station, div )[1] for station in stations ]
        page = "[" * join( stations[2:end], "," )[1:end-lim] * "]"
    end 

    data = jsontable(page)

    return DataFrame(data)
end

#   resAA = getDataAA( type=:METEO, source=:STATIONS )
#   resAA = getDataAA( type=:METEO, source=:SENSROS )
#   resAA = getDataAA( type=:AIRQUALITY, source=:STATIONS )
#   resAA = getDataAA( type=:AIRQUALITY, source=:SENSROS )


end # module