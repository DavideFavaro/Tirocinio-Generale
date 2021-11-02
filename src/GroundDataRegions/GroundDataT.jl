module GroundDataT
"""
Module for the download and processing of atmospheric data gathered by measuring stations located in Trentino, Italy
"""

#=
Trentino:
    Altri link utili:
        https://bollettino.appa.tn.it/aria/scarica
        https://dati.trentino.it/dataset/anagrafica-stazioni-meteo-stazioni-automatiche
        https://dati.trentino.it/dataset/dati-recenti-delle-stazioni-meteo


    Link:
        Dati qualità dell'aria (CSV):
            https://bollettino.appa.tn.it/aria/opendata/csv/last/
        Stazioni meteo automatiche(XML):
            http://dati.meteotrentino.it/service.asmx/listaStazioni
        Dati di una stazione automatica (XML):
            http://dati.meteotrentino.it/service.asmx/ultimiDatiStazione?codice=T0409
            ( Codice è l'id di una stazione )
=#



using CombinedParsers
using CombinedParsers.Regexp
using CSV
using DataFrames
using EzXML
using HTTP
using Revise


export getData,
       getRegionAttributes, getRegionIds, getRegionStationInfo


@syntax station = Repeat(
                      "<", re"[^>]+", ">",
                      Either( Numeric(Int64), Numeric(Float64), re"[^>]+" ),
                      Sequence( "</", re"[^>]+", ">" )
                  )



"""
    getRegionAttributes( [ type::Symbol=:METEO ] )

Obtain the names of the columns of the region's dataframe required by `GroundData.createMap`'s `attributes` parameter to create
`GroundData.standardize`'s `map` parameter
"""
function getRegionAttributes( type::Symbol=:METEO )
    return type == :METEO ?
               [ :info, :unit, :value, nothing, :date, :longitudine, :latitudine, :quota, nothing, nothing, :rmh ] :
               type == :AIRQUALITY ?
                   [ :Inquinante, :Unita_di_misura, :Valore, nothing, :Data, :Longitude, :Latitude, nothing, nothing, nothing ] :
                   throw( DomainError( type, "`type` must be either `:METEO` OR `:AIRQUALITY`" ) )
end



"""
    getRegionAttributes( [ type::Symbol=:METEO ] )

Obtain the names of the columns of the dataframe required for `GroundData.standardize`'s `bridge` parameter
"""
function getRegionIds( type::Symbol=:METEO )
    return type == :METEO ? :codice :
               type == :AIRQUALITY ? :NOME_STA : throw( DomainError( type, "`type` must be either `:METEO` OR `:AIRQUALITY`" ) )
end



"""
    getRegionStationInfo( [ type::Symbol=:METEO  ] )

Obtain the names of the columns of the region's stations dataframe required by `GroundData.createMap`'s `attributes` parameter to be used
in `GroundData.generateUuidsTable`
"""
function getRegionStationsInfo( type::Symbol=:METEO )
    return type == :METEO ?
               [ :codice, :nome, :longitudine, :latitudine ] :
               type == :AIRQUALITY ?
                   [ nothing, :NOME_STA, :Longitude, :Latitude  ] :
                   throw( DomainError( type, "`type` must be either `:METEO` OR `:AIRQUALITY`" ) )
end



"""
    getMeteoStationsData()

Obtain informations regarding the meteorological stations in Trentino
"""
function getMeteoStationsData()
    
    str = replace(
            replace(
                String( HTTP.get("http://dati.meteotrentino.it/service.asmx/listaStazioni").body ),
                r"\r\n *" => ""
            ),
            r"<fine */>" => "<fine>nothing</fine>"
          )

    stations_strs = split( str, r"</?anagrafica>", keepempty=false )[2:end-1]

    stations = station.(stations_strs)

    stations_dicts = [ Dict(
                        Symbol( String( attribute[2] ) ) => attribute[4] isa Number ? attribute[4] : String( attribute[4] )
                        for attribute in station 
                      ) for station in stations ]

    return DataFrame( stations_dicts )
end



"""
    getMeteoData( ids::AbstractVector{String} )

Obtain the meteorological data from the stations in `ids`
"""
function getMeteoData( ids::AbstractVector{String} )
    pages = [ String( HTTP.get("http://dati.meteotrentino.it/service.asmx/ultimiDatiStazione?codice=$id").body ) for id in ids ]

    units = Dict(
                "temperatura" => "°C",
                "pioggia" => "mm",
                "v" => "m/s",
                "vmax" => "m/s",
                "d" => "gN",
                "rsg" => "W/m²",
                "rh" => "%"
            )

    xmlpages = EzXML.parsexml.(pages)

    data_vect = []
    for (id, station) in zip(ids, xmlpages)
        stn = collect( eachelement( root(station) ) )[5:end]
        for attribute in stn
            for measurement in eachelement(attribute)
                msrmt = collect( eachelement(measurement) )
                for entry in msrmt[2:end]
                    dict = Dict(
                        :value => parse( Float64, entry.content ),
                        :codice => id,
                        :attribute => attribute.name,
                        :info => entry.name,
                        :unit => units[ entry.name ],
                        :date => msrmt[1].content
                    )
                    push!( data_vect, dict )
                end
            end
        end
    end
    df = DataFrame(data_vect)
    insertcols!( df, :rmh => "0m" )
    return df
end

# sdf = getMeteoStationsDataT()
# dct = getMeteoDataT( sdf[1:10,:codice] )



"""
"""
#=
Nel CSV si ha:
    - TRENTO PSC al posto di PARCO S. CHIARA
    - TRENTO VBZ al posto di VIA BOLZANO
    - ROVERETO LGP al poso di ROVERETO
    - RIVA GAR al psoto di RIVA DEL GARDA
    - BORGO VALSUGANA al posto di BORGO VAL
manca:
    - stazione A22 (AVIO)
=#
function getAQStationsData()
    return CSV.read( ".\\Dati stazioni\\rete_qual_air_TN.csv", DataFrame )
end



"""
    getAQData()
    
Return a `CSV.File` containing the data on air quality collected from measuring stations in Trentino Alto Adige
"""
function getAQData()
    data = HTTP.get( "https://bollettino.appa.tn.it/aria/opendata/csv/last/" )
    df = CSV.read( data.body, DataFrame )
    transform!( df, [:Data, :Ora] => ByRow( (date, time) -> date + Time( time == 24 ? 0 : time ) ) => :Data )
    select!( df, Not(4) )
    rename!( df, Symbol("Unit\xe0 di misura") => :Unita_di_misura )
    transform!( df, [:Unita_di_misura] => ByRow( x -> x = replace( replace( x, "\xb5" => "μ" ), "c" => "³" ) ) => :Unita_di_misura )
    transform!( df, [:Stazione] => ByRow( x -> x = uppercase(x) ) => :NOME_STA )
    return df
end

#   c = getAQData()



"""
    getData(; <keyword arguments> )

Obtain data of category `type` from `source`

# Arguments
 - `type::Symbol=:METEO`: defines the type of data to be downloaded may either be `:METEO` or `:AIRQUALITY`
 - `source::Symbol=:STATIONS`: defines if the data to be downloaded has to regard information on the stations or their actual measurements, as such may either be `:STATIONS` or `:SENSORS`
"""
function getData(; type::Symbol=:METEO, source::Symbol=:STATIONS )
    if type == :METEO
        stations = getMeteoStationsData()
        if source == :STATIONS
            return stations
        else
            return getMeteoData( stations[:, :codice] )
        end
    else
        if source == :STATIONS
            return getAQStationsData()
        else
            return getAQData()
        end
    end
end

#   ressta = getData( type=:METEO, source=:STATIONS )
#   ressen = getData( type=:METEO, source=:SENSORS )
#   ressta = getData( type=:AIRQUALITY, source=:STATIONS )
#   ressen = getData( type=:AIRQUALITY, source=:SENSORS )



end # module