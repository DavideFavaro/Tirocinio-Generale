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



using HTTP

using CSV
using DataFrames
using JSONTables


export getDataAA, getDataT


@enum Data_Type METEO=1 AIRQUALITY=2 
@enum Data_Source STATIONS=1 SENSORS=2




"""
    getDataAA(; type::Data_Type=METEO, source::Data_Source=STATIONS )

Obtain the data of the specified `type` regarding the stations themselves or the sensor's measureents
"""
function getDataAA(; type::Data_Type=METEO, source::Data_Source=STATIONS )
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



end # module