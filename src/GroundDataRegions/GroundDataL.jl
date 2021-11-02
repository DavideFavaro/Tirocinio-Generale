module GroundDataL
"""
Module for the download and processing of atmospheric data gathered by measuring stations located in Lombardia, Italy
"""

#=
Si possono ottenere dati dal link (ALL CSV):
    Stazioni:
        https://www.dati.lombardia.it/resource/wkyn-7szs.csv
    Stazioni meteo provincia di Milano:
        https://www.dati.lombardia.it/resource/my9y-8ykb.csv

    Stazioni meteo:
        https://www.dati.lombardia.it/resource/nf78-nj6b.csv
    Dati sensori meteo:
        https://www.dati.lombardia.it/resource/647i-nhxk.csv
    Stazioni aria:
        https://www.dati.lombardia.it/resource/ib47-atvt.csv
    Dati sensori aria:
        https://www.dati.lombardia.it/resource/nicp-bhqi.csv
Si possono ottenere direttamente CVS dei dati.
Dati sensori meteo sembra fare al caso nostro ma ha un campo che indica cosa rappresenta il dato che è encoded
=#


using CSV
using DataFrames
using HTTP


export getData,
       getRegionAttributes, getRegionIds, getRegionStationInfo




"""
    getRegionAttributes( [ type::Symbol=:METEO ] )

Obtain the names of the columns of the region's dataframe required by `GroundData.createMap`'s `attributes` parameter to create
`GroundData.standardize`'s `map` parameter
"""
function getRegionAttributes( type::Symbol=:METEO )
    return type == :METEO ?
               [ :tipologia, :unit_dimisura, :valore, nothing, :data, :lng, :lat, :quota, :stato, nothing, :rmh ] :
               type == :AIRQUALITY ?
                   [ :nometiposensore, :unitamisura, :valore, nothing, :data, :lng, :lat, :quota, :stato, nothing ] :
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
    return :idsensore
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
    return [ :idstazione, :nomestazione, :lng, :lat ]
end



"""
    getData(; <keyword arguments> )

Obtain data of category `type` from `source`

# Arguments
 - `type::Symbol=:METEO`: defines the type of data to be downloaded may either be `:METEO` or `:AIRQUALITY`
 - `source::Symbol=:STATIONS`: defines if the data to be downloaded has to regard information on the stations or their actual measurements, as such may either be `:STATIONS` or `:SENSORS`
"""
function getData(; type::Symbol=:METEO, source::Symbol=:STATIONS )
    str = type == :METEO ? ( source == :STATIONS ? "nf78-nj6b" : "647i-nhxk" ) : ( source == :STATIONS ? "ib47-atvt" : "nicp-bhqi" )
    data = HTTP.get( "https://www.dati.lombardia.it/resource/$str.csv" )
    df = CSV.read( data.body, DataFrame )

    if type ==:METEO && source == :SENSORS
        insertcols!( df, :rmh => "0m" ) 
    end

    return df
end

#   ressta = getData( type=:METEO, source=:STATIONS )
#   ressen = getData( type=:METEO, source=:SENSORS )
#   ressta = getData( type=:AIRQUALITY, source=:STATIONS )
#   ressen = getData( type=:AIRQUALITY, source=:SENSORS )



end # module