module GroundDataTAA
"""
Module for the download and processing of atmospheric data gathered by measuring stations located in Trentino Alto Adige, Italy
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
            ( Codie è l'id di una stazione )

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



using CombinedParsers
using CombinedParsers.Regexp

using CSV
using DataFrames
using JSONTables

using HTTP


export getDataAA, getDataT


@enum Data_Type METEO=1 AIRQUALITY=2 
@enum Data_Source STATIONS=1 SENSORS=2
@enum Data_Region T=1 AA=2


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



#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#                                                               DATI DELL'ALTO ADIGE
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

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



#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#                                                                   DATI DEL TRENTINO
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

"""
    getMeteoStationsData()

Obtain informations regarding the meteorological stations in Trentino
"""
function getMeteoStationsDataT()
    
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
function getMeteoDataT( ids::AbstractVector{String} )
    
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

    return attributes
end

# sdf = getMeteoStationsData()
# data = getMeteoData( sdf[1:3,:codice] )



"""
    getAQData()
    
Return a `CSV.File` containing the data on air quality collected from measuring stations in Trentino Alto Adige
"""
function getAQDataT()
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
        stations = getMeteoStationsDataT()
        if source == STATIONS
            return stations
        else
            return getMeteoDataT( stations[:, :codice] )
        end
    else
        return  DataFrame( getAQDataT() )
    end
end







end # module