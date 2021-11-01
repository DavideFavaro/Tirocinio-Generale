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
using Dates
using HTTP
using JSONTables


export getData,
       getRegionAttributes, getRegionIds, getRegionStationsInfo



"""
    getRegionAttributes( [ type::Symbol=:METEO ] )

Obtain the names of the columns of the region's dataframe required by `GroundData.createMap`'s `attributes` parameter to create
`GroundData.standardize`'s `map` parameter
"""
function getRegionAttributes( type::Symbol=:METEO )
    return type == :METEO ?
               [ :DESC_I, :UNIT, :VALUE, nothing, :DATE, :LONG, :LAT, :ALT, nothing, nothing, :rmh ] :
               type == :AIRQUALITY ?
                   [ :MCODE, nothing, :VALUE, nothing, :DATE, :LONG, :LAT, nothing, :FLAGS, nothing ] :
                   throw( DomainError( type, "`type` must be either `:METEO` OR `:AIRQUALITY`" ) )
end



"""
    getRegionAttributes( [ type::Symbol=:METEO ] )

Obtain the names of the columns of the dataframe required for `GroundData.standardize`'s `bridge` parameter
"""
function getRegionIds( type::Symbol=:METEO )
    if type != :METEO && type != :AIRQUALITY
        throw( DomainError( type, "`type` must be either `:METEO` OR `:AIRQUALITY`" ) )
    end
    return :SCODE
end



"""
    getRegionStationInfo( [ type::Symbol=:METEO  ] )

Obtain the names of the columns of the region's stations dataframe required by `GroundData.createMap`'s `attributes` parameter to be used
in `GroundData.generateUuidsTable`
"""
function getRegionStationsInfo( type::Symbol=:METEO )
    if type != :METEO && type != :AIRQUALITY
        throw( DomainError( type, "`type` must be either `:METEO` OR `:AIRQUALITY`" ) )
    end
    return [ :SCODE, :NAME_I, :LONG, :LAT ]
end



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

    data = DataFrame( jsontable(page) )
    
    if type == :METEO && source == :SENSORS
        insertcols!( data, :rmh => "0m" )
        transform!( data, [:DATE] => ByRow( x -> ismissing(x) ? missing : DateTime( x[1:19], "yyyy-mm-ddTH:M:S" ) ) => :DATE ) 
    end

    return data
end

#   ressta = getData( type=:METEO, source=:STATIONS )
#   ressen = getData( type=:METEO, source=:SENSORS )
#   ressta = getData( type=:AIRQUALITY, source=:STATIONS )
#   ressen = getData( type=:AIRQUALITY, source=:SENSORS )



end # module