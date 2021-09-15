module Test

using HTTP
using Downloads


# https://scihub.copernicus.eu/dhus/odata/v1
# https://scihub.copernicus.eu/dhus/search?q=*&rows=25
# https://scihub.copernicus.eu/apihub/odata/v1/
# https://scihub.copernicus.eu/dhus/odata/v1/Products('2b17b57d-fff4-4645-b539-91f305c27c69')/$value

# QUESTO EFFETTIVAMENTE SCARICA MA MANCA L'AUTENTICAZIONE
#https://scihub.copernicus.eu/dhus/odata/v1/Products('2b17b57d-fff4-4645-b539-91f305c27c69')/$value
# S3B_SL_2_LST____20210915T101845_20210915T102145_20210915T124701_0179_057_065_2160_LN2_O_NR_004  id di un file, tentare di usarlo in Product(...) ritorna not found
# 2e9876f9-d10b-4063-819b-def25ccaee79  codice ottenuto guardando l'XML su https://scihub.copernicus.eu/dhus/odata/v1/Products

out = "D:\\Vario\\Stage\\dt"

auth = HTTP.Base64.base64encode("davidefavaro:Tirocinio") #Trasformazione dei dati nel formato valido per la verifica

d = Downloads.download( "https://scihub.copernicus.eu/dhus/odata/v1/Products('2b17b57d-fff4-4645-b539-91f305c27c69')/$value",
                    out,
                    headers = [ "Authorization" => "Basic $auth" ], #info di autorizzazioone
                    progress = ( total, now ) -> println("$now / $total"), #Funzione che permette di controllare lo stato del download
                    verbose = true ) #Più info


end #module