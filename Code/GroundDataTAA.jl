module GroundDataTAA
#=
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
        ( Codie è l'id di una stazione )
=#

"""
Module for the download and processing of atmospheric data gathered by measuring stations located in Trentino Alto Adige, Italy
"""




using CombinedParsers
using CombinedParsers.Regexp

using CSV
using DataFrames

using Downloads
using HTTP



@syntax station = Repeat(
    "<", re"[^>]+", ">",
    Either( Numeric(Int64), Numeric(Float64), re"[^>]+" ),
    Sequence( "</", re"[^>]+", ">" )
)


@syntax data = Repeat(
    "<", re"[^>]+", ">",
    Repeat(
        "<", re"[^>]+", ">",
        Repeat(
            "<", re"[^>]+", ">",
            Either(
                Numeric(Int64),
                Numeric(Float64),
                re"[^>]+"
            ),
            "</", re"[^>]+", ">"
        ),
        "</", re"[^>]+", ">"
    ),
    "</", re"[^>]+", ">"
)



"""
    getAQData()
    
Return a `CSV.File` containing the data on air quality collected from measuring stations in Trentino Alto Adige
"""
function getAQData()
    data = HTTP.get( "https://bollettino.appa.tn.it/aria/opendata/csv/last/" )
    data_csv = CSV.File( data.body )
    return data_csv
end


#   c = getAQData()


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
















function getMeteoData( ids::AbstractVector{String} )
    
    strs = [ replace(
                replace(
                   String( HTTP.get("http://dati.meteotrentino.it/service.asmx/ultimiDatiStazione?codice=$id").body ),
                   r"\r\n *" => ""
                ),
                r" */>" => ">TBD</ >" 
             ) for id in ids ]

    data_strs = [ split(
                    split( str, r"</?datiOggi[^>]*>", keepempty=false )[2],
                    "</rain>",
                    keepempty=false
                  )[2] for str in strs ]


    # data_vect[i]                      i-th station
    # data_vect[i][j]                   i-th station's j-th attribute
    # data_vect[i][j][4]                j-th attribute's measurements
    # data_vect[i][j][4][l]             j-th attribute's l-th measurement
    
    # data_vect[i][j][4][l][2]          l-th measurement's name
    # data_vect[i][j][4][l][4]          l-th measurement's attributes  
    # data_vect[i][j][4][l][4][k][2]    l-th measurement's k-th attribute's name
    # data_vect[i][j][4][l][4][k][4]    k-th attribute's value
    data_vect = data.(data_strs)
    

    #   attributes = Dict(
    #                   Symbol( String( attribute[2] ) )
    #                   =>
    #                   DataFrame( [
    #                       push!(
    #                           Dict(
    #                               Symbol( String( value[2] ) ) => value[4] isa Number ? value[4] : String( value[4] )
    #                               for value in measurement[4]
    #                           ),
    #                           :station_id => id,
    #                           :attribute => String( attribute[2] )
    #                       )
    #                       for measurement in attribute[4]
    #                   ] )
    #                   for (id, station) in zip(ids, data_vect)
    #                   for attribute in station
    #                )








    
    attributes = Dict(
                     Symbol( String( attribute[2] ) ) => DataFrame()
                     for station in data_vect
                     for attribute in station
                 )

    for (id, station) in zip(ids, data_vect)
        for attribute in station
            append!(
                attributes[ Symbol( String( attribute[2] ) ) ],
                DataFrame( [
                    push!(
                        Dict(
                            Symbol( String( value[2] ) ) => value[4] isa Number ? value[4] : String( value[4] )
                            for value in measurement[4]
                        ),
                        :station_id => id,
                        :attribute => String( attribute[2] )
                    )
                    for measurement in attribute[4]
                ] )
            )
        end
    end

    #   for key in keys(attributes)
    #       attributes[key] = vcat( DataFrame.( attributes[key] ) )
    #   end

    return attributes
end


sdf = getMeteoStationsData()
data = getMeteoData( sdf[1:3,:codice] )


end # module