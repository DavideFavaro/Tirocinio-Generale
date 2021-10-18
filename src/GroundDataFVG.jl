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
using JSONTables

using HTTP


@enum Data_Type METEO=1 AIRQUALITY=2 
@enum Data_Source STATIONS=1 SENSORS=2


@syntax meteo_data = Repeat(
                        "<",
                        re"[^ >]+",
                        Optional( re" [^ =]+=\"[^ \"]+\"" ),
                        Optional( re" [^ =]+=\"[^\"]+\"" ),
                        ">",
                        Either( Numeric(Int64), Numeric(Float64), re"[^>]+" ),
                        "</", re"[^>]+", ">"
                     ) 



#   function getAQData()
#       codes = [ "qp5k-6pvm" , "d63p-pqpr", "7vnx-28uy", "t274-vki6", "2zdv-x7g2", "ke9b-p6z2" ]
#       params = [ "PM10", "PM2.5", "Ozono", "Monossido di carbonio", "Biossido di zolfo", "Biossido di Azoto" ]
#       data = [ DataFrame( CSV.File( HTTP.get( "https://www.dati.friuliveneziagiulia.it/resource/$code.csv" ).body ) ) for code in codes ]
#       for (df, param) in zip(data, params)
#           !in( "parametro", names(df) ) && insertcols!( df, "parametro" => param )
#   
#       end
#       dataframe = reduce( (x, y) -> vcat( x, y, cols=:intersect ), data )
#       return dataframe
#   end


function getMeteoData()
    resources = [ "ARI", "BAR", "BGG", "BIC", "BOA", "BOR", "BRU", "CAP", "CDP", "CER", "CHI", "CIV", "CMT", "COD", "COR",
                  "ENE", "FAG", "FOS", "FSP", "GEM", "GRA", "GRG", "GRM", "LAU", "LIG", "LSR", "MAT", "MGG", "MNF", "MUS",
                  "PAL", "PDA", "PIA", "213200", "POR", "PRD", "RIV", "SAN", "SGO", "SPN", "TAL", "TAR", "TOL", "TRI",
                  "UDI", "VIV", "ZON" ]
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
    vect = [
               Dict(
                    (
                        !ismissing( attribute[4] ) ?
                            !ismissing( attribute[3] ) ?
                                Symbol( String( attribute[4][5] ) * " ($( String( attribute[3][5] ) ))" ) :
                                Symbol( String( attribute[4][5] ) ) :
                            Symbol( String( attribute[2] ) )
                    ) =>
                    attribute[6] isa Number ?
                        attribute[6] :
                        String( attribute[6] )
                    for attribute in data
               )
               for data in data_parse
           ]

    keys_groups = unique( keys.(vect) )
    grouped_vect = [ filter( x -> keys(x) == ks, vect ) for ks in keys_groups ]
    dfs_vect = DataFrame.(grouped_vect)
    data = dfs_vect[1]
    for df in dfs_vect[2:end]
        append!( data, df, cols=:union )
    end

    return data
end

#   res = getMeteoData()




function getDataFVG(; type::Data_Type=METEO )
    if type == METEO
        return getMeteoData()
    else
        return getAQData()
    end
end


end # module