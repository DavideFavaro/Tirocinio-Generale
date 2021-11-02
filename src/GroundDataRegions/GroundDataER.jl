module GroundDataER
"""
Module for the download and processing of atmospheric data gathered by measuring stations located in Emilia Romagna, Italy
"""

#=
Link ai dati:
    Stazioni qualità aria - Link ottenuto dal bottone per il download (CSV):
        https://docs.google.com/spreadsheets/d/1UFdHVhzULdE-F4o1CFZ7JJ4hC7xJB1IBt_LUCB4K1Uk/export?format=csv
    Dati AQ:
        https://sdati-test.datamb.it/arex/


Per gli altri dati non ci sono link espliciti e nemmeno bottoni o opzioni per il download 
=#



export getData,
       getRegionAttributes, getRegionIds, getRegionStationInfo



"""
    getRegionAttributes( [ type::Symbol=:METEO ] )

Obtain the names of the columns of the region's dataframe required by `GroundData.createMap`'s `attributes` parameter to create
`GroundData.standardize`'s `map` parameter
"""
function getRegionAttributes( type::Symbol=:METEO )
    return type == :METEO ?
               [ nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing ] :
               type == :AIRQUALITY ?
                   [ nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing ] :
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
    return nothing
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
    return [ nothing, nothing, nothing, nothing ]
end



"""
"""
function getData(; type::Symbol=:METEO, source::Symbol=:STATIONS )
    println("Non ancora implementato")
    return nothing
end



end # module