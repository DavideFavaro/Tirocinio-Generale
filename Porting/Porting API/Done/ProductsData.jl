module Test

using HTTP
using Downloads

using CombinedParsers
using CombinedParsers.Regexp

using DataFrames
using Dates
#using Unitful

using CSV



test = [ "C:\\Users\\DAVIDE-FAVARO\\Desktop\\XML\\1.xml",
         "C:\\Users\\DAVIDE-FAVARO\\Desktop\\XML\\2.xml",
         "C:\\Users\\Lenovo\\Desktop\\XML\\1.xml",
         "C:\\Users\\Lenovo\\Desktop\\XML\\2.xml",
         "C:\\Users\\Lenovo\\Desktop\\XML\\Prod_Test.xml" ]







@syntax product = Repeat( Either( Sequence( "<",
                                            :type    => re"[^< >]+ ",
                                            :opening => re"name=\"[^<\">]+\">",
                                            :content => re"[^<>]+",
                                            :closing => re"</[^<>]+>" ),
                                  Sequence( re"<[^<>]+>", re"[^<>]+", re"</[^<>]+>" ),
                                  re"<[^<>]+>" ) )


#=

"""
    parseConvert( xmlType::AbstractString, value::AbstractString )

Given the xml type assigned to "value", convert the latter to its correct Julia type 
"""
function parseConvert( xmlType::AbstractString, value::AbstractString )
    if xmlType == "str"
        if !isnothing( match( re"^[0-9]+[.,][0-9]+$", value ) )
            return tryparse( Float64, value )
        elseif !isnothing( match( re"^[0-9]+$", value ) )
            return tryparse( Int64, value )
        else
            return value
        end
    end
    if xmlType == "int"
        return tryparse( Int64, value )
    end
    if xmlType == "double"
        return tryparse( Float64, value )
    end
    if xmlType == "date"
        return DateTime( value[1:end-1], "y-m-dTH:M:S.s" )
    end
    throw( DomainError( (xmlType, value), "Undefined type of $value" )  )
end





=#

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
        "C:\\Users\\Lenovo\\Desktop\\XML\\Test",
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


#=


"""
    getData( fileId::AbstractString, targetDirectory::AbstractString, authToken::AbstractString )

Download the archive containing the data identified by "fileId", into the chosen directory using the user authentication token
"""
function getData( fileId::AbstractString, targetDirectory::AbstractString, authToken::AbstractString )
    # Downloads.download( "https://scihub.copernicus.eu/dhus/odata/v1/Products('$fileId')/\$value",
    #                     targetDirectory*"\\$fileId.zip", #Destinazione
    #                     headers = [ "Authorization" => "Basic $authToken" ], #info di autorizzazioone
    #                     progress = ( total, now ) -> println("$(now/total*100)% ( $now / $total )"), #Funzione che permette di controllare lo stato del download
    #                     verbose = true ) #Più info

    Downloads.download( "https://scihub.copernicus.eu/dhus/odata/v1/Products('$fileId')/\$value",
                        targetDirectory*"\\$fileId.zip", #Destinazione
                        headers = [ "Authorization" => "Basic $authToken" ] ) #info di autorizzazioone
end

getData( "b57f225e-d288-4e4a-bb35-2a7eb75d60e4", out[3], authenticate("davidefavaro", "Tirocinio")  )


=#



#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#                                                            DOWNLOAD DEI FILE XML COLLEGATI AI PRODOTTI
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||



function getProductsPages( authToken::AbstractString, maxNumber::Union{Integer, Nothing} = nothing )
# Definition of the components of The URL
    aoi = "POLYGON((9.5000%2047.0000,%2014.0000%2047.0000,%2014.0000%2044.0000,%209.5000%2044.0000,%209.5000%2047.0000))"
    query2 = "[NOW-6MONTHS%20TO%20NOW]%20AND%20footprint:\"Intersects($aoi)\""
    query = "search?start=0&rows=0&q=ingestiondate:$query2"
    iob = IOBuffer()

# Get number of total products
    #Download the firs page (This way we can get the number we are interested in and also the first page that we would download anyway)
    Downloads.download( "https://scihub.copernicus.eu/dhus/$query", iob, headers = [ "Authorization" => "Basic $authToken" ] )
    lines = String(take!(iob))

    #count = tryparse( Int64, split( split( filter( line -> occursin( "totalResults", line ), lines )[1], ">" )[2], "<" )[1] )
    val = tryparse( Int64, split( split( lines, "totalResults>" )[2], "<" )[1] )
    count = [ val ÷ 100, val % 100 ]

    if !isnothing( maxNumber ) && maxNumber > 0 && maxNumber < val
        count[1] = maxNumber ÷ 100
        count[2] = maxNumber % 100
    end
    if count[1] > 0
        for i in range( 0, count[1] - 1, step=1 )
            query = "search?start=$(i*100)&rows=100&q=ingestiondate:$query2"
            Downloads.download( "https://scihub.copernicus.eu/dhus/$query", iob, headers = [ "Authorization" => "Basic $authToken" ] )
        end
    end
    if count[2] > 0
        query = "search?start=$(count[1] * 100)&rows=$(count[2])&q=ingestiondate:$query2"
        Downloads.download( "https://scihub.copernicus.eu/dhus/$query", iob, headers = [ "Authorization" => "Basic $authToken" ] )
    end
    return iob
end



#NON FUNZIONA, PER QUALCHE MOTIVO IL SECONDO SPLIT NON DIVIDE LE PAGINE IN PRODOTTI

function getPageProducts( fileIO::IO )
# Get the downloaded pages containing the XML representations of the products
    original = String(take!(fileIO))
    println(typeof(original))
    println(length(original))
    println()

# Split the downloaded files into pages   
    pages = split( original, r"<\?xml [^<>]+><[^<>]+>", keepempty=false )
    println(typeof(pages))
    println(length(pages))
    println()
    return pages

# Split the result in a vector of ready-to-be-parsed strings eliminating "<entry>" and "</entry>" tags in the process 
    vector = split.( pages, r"<.?entry>", keepempty=false )
    println(typeof(vector))
    println(length(vector))
    println()

# Parse the strings
    products = [ product(x)[8:end] for x in vector ]


    #for i in 1:length(products)
    #    println()
    #    println("PRODOTTO $i:\n")
    #    for p in products[i]
    #        println( "$( join(p[:opening][10]) ) : $( join(p[:content]) )" )
    #    end
    #    println()
    #end

# Generate a vector of dictionaries containing the details for each product of the original page
    return [ Dict( Symbol(join(p[:opening][7])) => parseConvert( join(p[:type][1]), join(p[:content]) ) for p in products[i] ) for i in 1:length(products) ]
end


io = getProductsPages( authenticate("davidefavaro","Tirocinio"), 200 )
dict = getPageProducts( io )






















function getProductsDF( targetDirectory::AbstractString, authToken::AbstractString, maxNumber::Union{Integer, Nothing} = nothing )
# Download "maxNumber" pages in "targetDirectory"
    io = getProductsPages(  authToken, maxNumber )

# Create a vector of dictionaries of the products
    dict_vect = getPageProducts( io )

# Obtain the existing subsets of attributes of the products
    keys_groups = unique( keys.(dict_vect) )
# Divide the dictionaries in groups homogeneus on their attributes
    grouped_vect = [ filter( x -> keys(x) == ks, dict_vect ) for ks in keys_groups ]
# Turn each group in a DataFrame 
    dfs_vect = DataFrame.(grouped_vect)
# Merge all Dataframes using append to create a Dataframe with the union of all the possible columns and the right "missing" values 
    data = dfs_vect[1]
    for df in dfs_vect[2:end]
        append!( data, df, cols=:union )
    end

    return data
end


























"""
    getProducts( targetDirectory::AbstractString, authToken::AbstractString[, maxNumber::Union{Integer, Nothing} ] )

Obtain "maxNumber" pages, eachone containing the XML description of 100 products, through "authToken" and store them in "targetDirectory"

#Arguments
- `targetDirectory::AbstractString`: path to the directory of destination
- `authToken::AbstractString`: Authentication token, obtained through "authenticate()" function
- `maxNumber::Integer`: maximum number of pages to download
"""
function getProductsPages( targetDirectory::AbstractString, authToken::AbstractString, maxNumber::Union{Integer, Nothing} = nothing )
# Definition of the components of The URL
    aoi = "POLYGON((9.5000%2047.0000,%2014.0000%2047.0000,%2014.0000%2044.0000,%209.5000%2044.0000,%209.5000%2047.0000))"
    query2 = "[NOW-6MONTHS%20TO%20NOW]%20AND%20footprint:\"Intersects($aoi)\""
    query = "search?start=0&rows=0&q=ingestiondate:$query2"
    file0 = targetDirectory*"\\0.xml"

# Get number of total products
    #Download the firs page
    Downloads.download( "https://scihub.copernicus.eu/dhus/$query", file0, headers = [ "Authorization" => "Basic $authToken" ] )
    line = readlines( file0 )[9]
    rm( file0 )

    #Filter its contents to obtain the total number ofproducts
    #count = tryparse( Int64, split( split( filter( line -> occursin( "totalResults", line ), lines )[1], ">" )[2], "<" )[1] )
    val = tryparse( Int64, split( split( line, ">" )[2], "<" )[1] )
    count = [ val ÷ 100, val % 100 ]


    if !isnothing( maxNumber ) && maxNumber > 0 && maxNumber < val
        count[1] = maxNumber ÷ 100
        count[2] = maxNumber % 100
    end

    if count[1] > 0
        for i in range( 0, count[1] - 1, step=1 )
            query = "search?start=$(i*100)&rows=100&q=ingestiondate:$query2"
            Downloads.download( "https://scihub.copernicus.eu/dhus/$query", targetDirectory*"\\$(i + 1).xml", headers = [ "Authorization" => "Basic $authToken" ] )
        end
    end
    if count[2] > 0
        query = "search?start=$(count[1] * 100)&rows=$(count[2])&q=ingestiondate:$query2"
        Downloads.download( "https://scihub.copernicus.eu/dhus/$query", targetDirectory*"\\$(count[1] + 1).xml", headers = [ "Authorization" => "Basic $authToken" ] )
    end
end


# getProductsPages( out[5], authenticate("davidefavaro","Tirocinio"), 300 )



"""
    getPageProducts( path::AbstractString )

Given the path to a file downloaded with "getProductsPages()", return an array of the dictionaries containing the informatiions on each of the products of the page 
"""
function getPageProducts( path::AbstractString )
# Get the downloaded page containing the XML representations of the products
    original = readlines( path )
# Join all the rows in a single string
    string = join( original[19:end-1] )
# Split the result in a vector of ready-to-be-parsed strings eliminating "<entry>" and "</entry>" tags in the process 
    vector = split( string, r"<.?entry>", keepempty=false )
# Parse the strings
    products = [ product(x)[8:end] for x in vector ]

# For each product print its details
    #for i in 1:length(products)
    #    println()
    #    println("PRODOTTO $i:\n")
    #    for p in products[i]
    #        println( "$( join(p[:opening][10]) ) : $( join(p[:content]) )" )
    #    end
    #    println()
    #end

# Generate a vector of dictionaries containing the details for each product of the original page
    return [ Dict( Symbol(join(p[:opening][7])) => parseConvert( join(p[:type][1]), join(p[:content]) ) for p in products[i] ) for i in 1:length(products) ]
end

# getPageProducts( out[5] )


"""
    getProductsDF( targetDirectory::AbstractString, authToken::AbstractString, maxNumber::Union{Integer, Nothing} = nothing )

Generate in "targetDirectory" the DataFrame containing the data of "maxNumber" undreds of products using "authToken"
"""
function getProductsDF( targetDirectory::AbstractString, authToken::AbstractString, maxNumber::Union{Integer, Nothing} = nothing )
# Download "maxNumber" pages in "targetDirectory"
    getProductsPages( targetDirectory, authToken, maxNumber )
    
# Check how many pages have been downloaded
    len = length( readdir( targetDirectory ) )

# For each of the pages obtain the dictionaries of its products and add them to "res"
    dict_vect = []
    for i in 1:len
        target = targetDirectory*"\\$i.xml"
        dict_vect = vcat( dict_vect, getPageProducts( target ) )
        #rm( target )
    end

# Obtain the existing subsets of attributes of the products
    keys_groups = unique( keys.(dict_vect) )
# Divide the dictionaries in groups homogeneus on their attributes
    grouped_vect = [ filter( x -> keys(x) == ks, dict_vect ) for ks in keys_groups ]
# Turn each group in a DataFrame 
    dfs_vect = DataFrame.(grouped_vect)
# Merge all Dataframes using append to create a Dataframe with the union of all the possible columns and the right "missing" values 
    data = dfs_vect[1]
    for df in dfs_vect[2:end]
        append!( data, df, cols=:union )
    end

    return data
end

df = getProductsDF( out[3], authenticate("davidefavaro","Tirocinio"), 20 )                   



"""
    saveProductsDF( targetDirectory::AbstractString, data::DataFrame; overwrite::Bool=false )

Save "data" in "targetDirectory" if not already existing or if "overwrite" is true, otherwise append its content to "data.csv"
"""
function saveProductsDF( targetDirectory::AbstractString, data::DataFrame; overwrite::Bool=false )
    CSV.write( targetDirectory*"\\data.csv", data, append = !overwrite && in("data.csv", readdir(targetDirectory)) )
end

saveProductsDF( out[2], df )



"""
    createUnitsDF( columns::AbstractVector{AbstractString}, units::AbstractVector{AbstractString} )

Given the array of column names and the array of their respective units of measure, create a df associating both
"""
function createUnitsDF( columns::AbstractVector{AbstractString}, units::AbstractVector{AbstractString} )
    return Dataframe( Dict( columns .=> units ) )
end












end #module