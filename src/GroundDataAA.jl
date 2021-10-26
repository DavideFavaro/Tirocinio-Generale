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


str = occursin( "GroundDataAA.jl", @__FILE__ ) ? "" : "src\\"
include("$(@__DIR__)\\$(str)Global.jl")


export getData



"""
    getDataAA(; type::Symbol=:METEO, source::Symbol=STATIONS )

Obtain the data of the specified `type` regarding the stations themselves or the sensor's measureents
"""
function getData(; type::Data_Type=METEO, source::Data_Source=STATIONS )
    opt1 = type == METEO ? "meteo/v1" : "airquality"
    opt2 = source == STATIONS ? "stations" : type == METEO ? "sensors" : "timeseries"

    page = String( HTTP.get( "http://dati.retecivica.bz.it/services/$opt1/$opt2" ).body )

    if source == STATIONS
        if type == METEO
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

#   resAA = getDataAA( type=METEO, source=STATIONS )
#   resAA = getDataAA( type=type[1], source=source[2])
#   resAA = getDataAA( type=type[2], source=source[1] )
#   resAA = getDataAA( type=type[2], source=source[2])


end # module