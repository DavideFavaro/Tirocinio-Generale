module Test

using HTTP
using Downloads



"""
    authenticate( username::AbstractString, password::AbstractString[, type::AbstractString ] )

Create an authentication token of type "type" for the user
"""
function authenticate( username::AbstractString, password::AbstractString, type::AbstractString = "Base" )
    if type == "Base"
        return HTTP.Base64.base64encode("$username:$password") #Trasformazione dei dati nel formato valido per la verifica
    end
end


out = "D:\\Vario\\Stage\\data.zip"
out2 = "C:\\Users\\DAVIDE-FAVARO\\Desktop\\data"
out3 = "C:\\Users\\DAVIDE-FAVARO\\Desktop\\XML"
auth = authenticate("davidefavaro","Tirocinio") 
address = "https://scihub.copernicus.eu/dhus"




#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#                                                    DOWNLOAD DEL .ZIP CORRISPONDENTE AD UN DATO ID
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||


# https://scihub.copernicus.eu/dhus/odata/v1
# https://scihub.copernicus.eu/dhus/search?q=*&rows=25
# https://scihub.copernicus.eu/dhus/odata/v1/Products('2b17b57d-fff4-4645-b539-91f305c27c69')/$value

# S3B_SL_2_LST____20210915T101845_20210915T102145_20210915T124701_0179_057_065_2160_LN2_O_NR_004  Nome di un file, tentare di usarlo in Product(...) ritorna not found
# 1def7d25-ecd6-49ef-8e3a-243f35857951 uuid file 240 Mb
# b0526ddf-c3f3-4065-8dc5-1891ffc8e326 uuid file 165 Mb

# https://scihub.copernicus.eu/dhus/odata/v1/Products('2b17b57d-fff4-4645-b539-91f305c27c69')/$value
# https://scihub.copernicus.eu/dhus/odata/v1/Products?filter=Name eq 'S3B_SL_2_LST____20210915T101845_20210915T102145_20210915T124701_0179_057_065_2160_LN2_O_NR_004'

"""
    getData( fileId::AbstractString, targetDirectory::AbstractString, authToken::AbstractString )

Download the archive containing the data identified by "fileId", into the chosen directory using the user authentication token
"""
function getData( fileId::AbstractString, targetDirectory::AbstractString, authToken::AbstractString )

    Downloads.download( "https://scihub.copernicus.eu/dhus/odata/v1/Products('$fileId')",
                        targetDirectory*"\\data.zip", #Destinazione
                        headers = [ "Authorization" => "Basic $authToken" ], #info di autorizzazioone
                        progress = ( total, now ) -> println("$(now/total*100)% ( $now / $total )"), #Funzione che permette di controllare lo stato del download
                        verbose = true ) #Più info
end

#getData( "68fc21ec-d9e1-44f1-be45-f487031423a1", "C:\\Users\\DAVIDE-FAVARO\\Desktop", authenticate("davidefavaro", "Tirocinio")  )








#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#                                                            DOWNLOAD DEI FILE XML COLLEGATI AI PRODOTTI
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||


# https://scihub.copernicus.eu/dhus/search?start=0&rows=1&q=ingestiondate:
#   [NOW-6MONTHS%20TO%20NOW]%20AND%20footprint:%22Intersects(POLYGON((9.5000%2047.0000,%2014.0000%2047.0000,%2014.0000%2044.0000,%209.5000%2044.0000,%209.5000%2047.0000)))%22


#PER IL MOMENTO RITORNA SOLO IL TOTAL COUNT, INOLTRE DEVO TESTARE QUNATO E' VELOCE COMBINEDPERSER RISPETTO A COME HO FATTO IO

function getProducts( targetDirectory::AbstractString, authToken::AbstractString, maxNumber::Integer )
    start = 0
    rows = 0
    address = "https://scihub.copernicus.eu/dhus"
    aoi = "POLYGON((9.5000%2047.0000,%2014.0000%2047.0000,%2014.0000%2044.0000,%209.5000%2044.0000,%209.5000%2047.0000))"
    query2 = "[NOW-6MONTHS%20TO%20NOW]%20AND%20footprint:\"Intersects($aoi)\""
    query = "/search?start=$start&rows=$rows&q=ingestiondate:$query2"

    # Get number of total products 
    #Downloads.download( "$address$query",
    #                    targetDirectory*"\\data.xml",
    #                    headers = [ "Authorization" => "Basic $(authenticate("davidefavaro","Tirocinio"))" ], #info di autorizzazioone
    #                    progress = ( total, now ) -> println("$(now/total*100)% ( $now / $total )"), #Funzione che permette di controllare lo stato del download
    #                    verbose = true ) #Più info
    Downloads.download( "$address$query", targetDirectory*"\\data.xml", headers = [ "Authorization" => "Basic $(authenticate("davidefavaro","Tirocinio"))" ] )
    lines = readlines( targetDirectory*"\\data.xml" )
    count = tryparse( Int64, split( split( filter( line -> occursin( "totalResults", line ), lines )[1], ">" )[2], "<" )[1] )

    # DA AGGIUNGERE LA PARTE CHE ITERA SU TUTTI GLI XML 

    return count
end

    


getProducts( out3, authenticate("davidefavaro","Tirocinio"), 10 )
                        
                        


                        
                            
























end #module