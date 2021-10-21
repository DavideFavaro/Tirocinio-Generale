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



using HTTP

using CombinedParsers
using CombinedParsers.Regexp

using CSV
using DataFrames
using Dates

using EzXML

using Revise

export getDataT


@enum Data_Type METEO=1 AIRQUALITY=2 
@enum Data_Source STATIONS=1 SENSORS=2


@syntax station = Repeat(
                      "<", re"[^>]+", ">",
                      Either( Numeric(Int64), Numeric(Float64), re"[^>]+" ),
                      Sequence( "</", re"[^>]+", ">" )
                  )



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
                        :value => entry.content,
                        :station_id => id,
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
    return DataFrame(data_vect)
end

# sdf = getMeteoStationsDataT()
# dct = getMeteoDataT( sdf[1:10,:codice] )



"""
    getAQData()
    
Return a `CSV.File` containing the data on air quality collected from measuring stations in Trentino Alto Adige
"""
function getAQData()
    data = HTTP.get( "https://bollettino.appa.tn.it/aria/opendata/csv/last/" )

    return CSV.File( data.body )
end

#   c = getAQData()



"""
    getDataT(; type::Data_Type=METEO, source::Data_Source=STATIONS )

Obtain informations on the `type` stations or their sensor's data
"""
function getDataT(; type::Data_Type=METEO, source::Data_Source=STATIONS )
    if type == METEO
        stations = getMeteoStationsData()
        if source == STATIONS
            return stations
        else
            return getMeteoData( stations[:, :codice] )
        end
    else
        return  DataFrame( getAQData() )
    end
end

#   df = getDataT( type=METEO, source=STATIONS )
#   df = getDataT( type=METEO, source=SENSORS )
#   df = getDataT( type=AIRQUALITY, source=STATIONS )



end # module