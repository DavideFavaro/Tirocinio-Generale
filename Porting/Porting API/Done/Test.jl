module Test

using HTTP
using Downloads

using CombinedParsers
using CombinedParsers.Regexp


# TEST ROW:
#   test[1]
#   "C:\\Users\\DAVIDE-FAVARO\\Desktop\\XML\\Row_Test.xml"
#   "<opensearch:totalResults>13265</opensearch:totalResults>"

# TEST HEADER:
#   test[2]"
#   "C:\\Users\\DAVIDE-FAVARO\\Desktop\\XML\\Head_Test.xml"

# TEST PRODUCT:
#   test[3]
#   "C:\\Users\\DAVIDE-FAVARO\\Desktop\\XML\\Prod_Test.xml"

# TEST DOCUMENT:
#   "C:\\Users\\DAVIDE-FAVARO\\Desktop\\XML\\1.xml"

source = [ "C:\\Users\\DAVIDE-FAVARO\\Desktop\\XML\\",
           "C:\\Users\\Lenovo\\Desktop\\XML\\" ]

test = [ "Row_Test.xml",
         "Head_Test.xml",
         "Prod_Test.xml",
         "1.xml" ]



# Credo che i Repeat() senza indicazione numerica causino loop infiniti
# In 1.xml il minimo di attributi era 28 il massimo 48

# Match su tutte le stringhe contenute tra due tag dal contenuto qualunque e tag singoli ma non su "<entry>" e "</entry>"
@syntax row = Either( Sequence( re"<.+>",
                                Either( Numeric(Int), re".+" ),
                                re"</.+>\n" ),
                      re"<(?!(/?entry)).*>\n" )

# Match sull eprime 17 "row"
#@syntax header = Repeat( 17, row )


@syntax product = Sequence( "<entry>\n",
                            Repeat(row),
                            "<str name=\"gmlfootprint\">",
                            re".+",
                            "</str>",
                            Repeat( row ),
                            "</entry>\n" )


a = readline( source[2]*test[1] )
b = join( readlines( source[2]*test[2] ), '\n' ) * '\n'
c = join( readlines( source[2]*test[3] ), '\n' ) * '\n'
d = readlines( source[2]*test[4] )

row(a)
header( join(b) )
product( join(c) )
document( joni(d) )


for line in b
    println( ( row(line) ) )
end 






"""
    authenticate( username::AbstractString, password::AbstractString[, type::AbstractString ] )

Create an authentication token of type "type" for the user
"""
function authenticate( username::AbstractString, password::AbstractString, type::AbstractString = "Base" )
    if type == "Base"
        return HTTP.Base64.base64encode("$username:$password") #Trasformazione dei dati nel formato valido per la verifica
    end
end


out = [ "D:\\Vario\\Stage",
        "C:\\Users\\Lenovo\\Desktop\\XML",
        "C:\\Users\\DAVIDE-FAVARO\\Desktop",
        "C:\\Users\\DAVIDE-FAVARO\\Desktop\\XML" ]






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

#getData( "68fc21ec-d9e1-44f1-be45-f487031423a1", out[3], authenticate("davidefavaro", "Tirocinio")  )








#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#                                                            DOWNLOAD DEI FILE XML COLLEGATI AI PRODOTTI
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||


"""
    getProducts( targetDirectory::AbstractString, authToken::AbstractString[, maxNumber::Union{Integer, Nothing} ] )

Obtain "maxNumber" pages, eachone containing the XML description of 100 products, through "authToken" and store them in "targetDirectory"

#Arguments
- `targetDirectory::AbstractString`: path to the directory of destination
- `authToken::AbstractString`: Authentication token, obtained through "authenticate()" function
- `maxNumber::Integer`: maximum number of pages to download
"""
function getProducts( targetDirectory::AbstractString, authToken::AbstractString, maxNumber::Union{Integer, Nothing} = nothing )
# Definition of the components of The URL
    aoi = "POLYGON((9.5000%2047.0000,%2014.0000%2047.0000,%2014.0000%2044.0000,%209.5000%2044.0000,%209.5000%2047.0000))"
    query2 = "[NOW-6MONTHS%20TO%20NOW]%20AND%20footprint:\"Intersects($aoi)\""
    query = "search?start=0&rows=0&q=ingestiondate:$query2"
    file0 = targetDirectory*"\\0.xml"

# Get number of total products
    #Download a short XML file
    Downloads.download( "https://scihub.copernicus.eu/dhus/$query", file0, headers = [ "Authorization" => "Basic $authToken" ] )
    lines = readlines( file0 )
    #Filter its contents to obtain the total number ofproducts
    count = tryparse( Int64, split( split( filter( line -> occursin( "totalResults", line ), lines )[1], ">" )[2], "<" )[1] )
    #After retrieving the total count, remove the XML file
    #rm( file0 )


# Get the XML files of the existing products
    #Check if maxNumber has an actual value, if so use it to limit the number of collected pages
        #One page holds the XML specifications of 100 products
    if !isnothing(maxNumber) && maxNumber > 0
        maxNumber *= 100
        if maxNumber < count
            count = maxNumber
        end
    end
    #Get the desired number of pages
    for i in range( 0, count, step=100 )
        query = "search?start=$i&rows=100&q=ingestiondate:$query2"
        Downloads.download( "https://scihub.copernicus.eu/dhus/$query", targetDirectory*"\\$((i ÷ 100) + 1).xml", headers = [ "Authorization" => "Basic $authToken" ] )
    end
end


getProducts( out[4], authenticate("davidefavaro","Tirocinio"), 1 )
                        
                        


                        
                            
























end #module