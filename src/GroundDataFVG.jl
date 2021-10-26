module GroundDataFVG
"""
Module for the download and processing of atmospheric data gathered by measuring stations located in Friuli Venezia Giulia, Italy 
"""

#=
Accesso risorse Friuli:
    https://www.dati.friuliveneziagiulia.it/resource/qp5k-6pvm.csv

    Qualità dell'aria:
        PM10:                       qp5k-6pvm
        PM2.5:                      d63p-pqpr
        Ozono:                      7vnx-28uy
        Monossido di carbonio:      t274-vki6
        Biossido di zolfo:          2zdv-x7g2
        Biossido di azoto:          ke9b-p6z2
        Pollini:                    svph-8w2g
    
    Meteo:
        Elenco sensori:             498i-2j88
        Previsioni e dati:          j654-ykm6 (Probabilmente non funzionante, il link è esterno e non ha il download)
        Dati stazioni:              4wxn-35av (Probabilmente non funzionante, il link è esterno e non ha il download)


Link:
    Dati stazioni:
        https://dev.meteo.fvg.it/xml/stazioni/GRA.xml
=#



using CombinedParsers
using CombinedParsers.Regexp
using CSV
using DataFrames
using HTTP
using JSONTables


export getDataFVG


@syntax meteo_data = Repeat(
                        "<",
                        re"[^ >]+",
                        Optional( re" [^ =]+=\"[^ \"]+\"" ),
                        Optional( re" [^ =]+=\"[^\"]+\"" ),
                        ">",
                        Either( Numeric(Int64), Numeric(Float64), re"[^>]+" ),
                        "</", re"[^>]+", ">"
                     ) 


const attributes = Dict(
                      :METEO      => [ :param, :unit, :value, nothing, :observation_time, :longitude, :latitude, :station_altitude, :rel_measure_height#=, nothing, nothing=#],
                      :AIRQUALITY => [ :parametro, :unita_misura, :value, nothing, :data_misura, :longitudine, :latitudine, nothing, :dati_insuff, nothing ]
                   )

const ids = Dict(
                :METEO      => :nome,
                :AIRQUALITY => nothing
            )



"""
"""
function getMeteoStationsData()
    return CSV.read( ".\\Dati stazioni\\stazioni_meteoclimatiche-FVG.csv", DataFrame )
end



#   prec_type
#       0:nulla; 1:pioggia; 2:pioggia e neve; 3:neve
#   cloudiness
#       0:n.d.; 1:sereno; 2:poco nuvoloso; 3:variabile; 4:nuvoloso; 5:coperto
"""
"""
function getMeteoData()
    resources = [ "ARI", "BAR", "BGG", "BIC", "BOA", "BOR", "BRU", "CAP", "CDP", "CER", "CHI", "CIV", "CMT", "COD", "COR",
                  "ENE", "FAG", "FOS", "FSP", "GEM", "GRA", "GRG", "GRM", "LAU", "LIG", "LSR", "MAT", "MGG", "MNF", "MUS",
                  "PAL", "PDA", "PIA", "213200", "POR", "PRD", "RIV", "SAN", "SGO", "SPN", "TAL", "TAR", "TOL", "TRI",
                  "UDI", "VIV", "ZON" ]
    #   prec_type = [ "nulla", "pioggia", "pioggia e neve", "neve" ]
    #   cloudiness = [ "n.d.", "sereno", "poco nuvoloso", "variabile", "nuvoloso", "coperto" ]
    
    data_str = []
    for res in resources
        try
            page = HTTP.get("https://dev.meteo.fvg.it/xml/stazioni/$res.xml")
            push!( data_str, String(page.body) )
        catch e
            if !isa( e, HTTP.ExceptionRequest.StatusError )
                throw(e)
            end
        end
    end
    
    data_split = @. replace( replace( getindex( split( data_str, r"</?meteo_data>" ), 2 ), r"\n *" => "" ), r"<!--[^-]+-->" => "" )
    data_parse = meteo_data.(data_split)

    vect = []
    for data in data_parse
        # Attributes of a single station, they are shared between all the parameters measured by the station
        others = [
                     Symbol( String( attribute[2] ) ) => attribute[6] isa Number ? attribute[6] : String( attribute[6] )
                     for attribute in data[1:6]
                 ]
        for attribute in data[7:end]
            # Parameters measured by a single station
            dict = push!(
                       Dict(
                          :param => !ismissing( attribute[4] ) ? String( attribute[4][5] ) : String( attribute[2] ),
                          :value => attribute[6] isa Number ? attribute[6] : String( attribute[6] ),
                          :unit => !ismissing(attribute[3]) ? String( attribute[3][5] ) : missing
                       ),
                       others...
                   )
            push!( vect, dict )
        end
    end

    df = select!( DataFrame(vect), Not(2, 7) )

    rel_heights = [ length( split( param, " a " ) ) == 2 ? split( param, " a " )[2] : "0m" for param in df[:, :param] ]

    transform!( df, [:station_name] => ByRow( x -> x = uppercase(x) ) => :nome )

    insertcols!( df, :rel_measure_height => rel_heights )

    return df
end



"""
media_giornaliera, media_giornaliera, media_oraria_max, media_mobile_8h_max, media_giornaliera, media_oraria_max
"""
function getAQData()
    codes = [ "qp5k-6pvm" , "d63p-pqpr", "7vnx-28uy", "t274-vki6", "2zdv-x7g2", "ke9b-p6z2" ]
    params = [ "PM10", "PM2.5", "Ozono", "Monossido di carbonio", "Biossido di zolfo", "Biossido di Azoto" ]
    val_types = [ "media_giornaliera", "media_giornaliera", "media_oraria_max", "media_mobile_8h_max", "media_giornaliera", "media_oraria_max" ]
    data = [ DataFrame( CSV.File( HTTP.get( "https://www.dati.friuliveneziagiulia.it/resource/$code.csv" ).body ) ) for code in codes ]
    for (i, (df, param)) in enumerate(zip(data, params))
        rename!( df, val_types[i] => "value" )
        insertcols!(df, :type => val_types[i] )
        !in( "parametro", names(df) ) && insertcols!( df, :parametro => param )
    end
    dataframe = reduce( (x, y) -> vcat( x, y, cols=:intersect ), data )
    return dataframe
end



"""
"""
function getData(; type::Symbol=:METEO, source::Symbol=:STATIONS )
    if type == :METEO
        if source == :STATIONS
            return getMeteoStationsData()
        else
            return getMeteoData()
        end
    else
        return getAQData()
    end
end

#   resFVG = getData( type=:METEO, source=:STATIONS )
#   resFVG = getData( type=:METEO, source=:SENSORS )
#   resFVG = getData( type=:AIRQUALITY )



end # module