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


include("./src/Global.jl")


export getDataL



"""
    getDataL(; type::Data_Type, source::Data_Source )

Download the data specified by `type` and `source`, returning a CSV file
"""
function getDataL(; type::Data_Type=METEO, source::Data_Source=STATIONS )

    str = type == METEO ? ( source == STATIONS ? "nf78-nj6b" : "647i-nhxk" ) : ( source == STATIONS ? "ib47-atvt" : "nicp-bhqi" )
    data = HTTP.get( "https://www.dati.lombardia.it/resource/$str.csv" )
    data_csv = CSV.File( data.body )

    return DataFrame(data_csv)
end

#   resLst = getDataL( type=METEO, source=STATIONS )
#   resLse = getDataL( type=METEO, source=SENSORS )
#   resL = getDataL( type=AIRQUALITY, source=STATIONS )
#   resL = getDataL( type=AIRQUALITY, source=SENSORS )



end # module