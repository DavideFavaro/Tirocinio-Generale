module GroundDataLombardia

"""
Module for the download and processing of atmospheric data gathered by measuring stations located in Lombardia, Italy
"""

using CSV
using Downloads
using HTTP

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


@enum Data_Type meteo=1 airquality=2 
@enum Data_Source stations=1 sensors=2


"""
    getDataL(; type::Data_Type, source::Data_Source )

Download the data specified by `type` and `source`, returning a CSV file
"""
function getDataL(; type::Data_Type=meteo, source::Data_Source=stations )

    str = type == meteo ? ( source == stations ? "nf78-nj6b" : "647i-nhxk" ) : ( source == stations ? "ib47-atvt" : "nicp-bhqi" )
    data = HTTP.get( "https://www.dati.lombardia.it/resource/$str.csv" )
    data_csv = CSV.File( data.body )

    return data_csv
end


data = getDataL()

end # module