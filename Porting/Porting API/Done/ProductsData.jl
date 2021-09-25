module Test

using HTTP
using Downloads

using CombinedParsers
using CombinedParsers.Regexp

using DataFrames
using Dates
#using Unitful
using CSV
using ZipFile
using NCDatasets



test = [ "C:\\Users\\DAVIDE-FAVARO\\Desktop\\XML\\1.xml",
         "C:\\Users\\DAVIDE-FAVARO\\Desktop\\XML\\2.xml",
         "C:\\Users\\Lenovo\\Desktop\\XML\\1.xml",
         "C:\\Users\\Lenovo\\Desktop\\XML\\2.xml",
         "C:\\Users\\Lenovo\\Desktop\\XML\\Prod_Test.xml" ]

out = [ "D:\\Vario\\Stage",
        "C:\\Users\\Lenovo\\Desktop\\XML",
        "C:\\Users\\Lenovo\\Desktop\\XML\\Test",
        "C:\\Users\\DAVIDE-FAVARO\\Desktop",
        "C:\\Users\\DAVIDE-FAVARO\\Desktop\\XML" ]



@syntax product = Repeat( Either( Sequence( "<",
                                            :type    => re"[^< >]+ ",
                                            :opening => re"name=\"[^<\">]+\">",
                                            :content => re"[^<>]+",
                                            :closing => re"</[^<>]+>" ),
                                  Sequence( re"<[^<>]+>", re"[^<>]+", re"</[^<>]+>" ),
                                  re"<[^<>]+>" ) )



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



"""
    authenticate( username::AbstractString, password::AbstractString[, type::AbstractString ] )

Create an authentication token of type "type" for the user
"""
function authenticate( username::AbstractString, password::AbstractString, type::AbstractString = "Base" )
    if type == "Base"
        return HTTP.Base64.base64encode("$username:$password") #Trasformazione dei dati nel formato valido per la verifica
    end
end




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






#dir = "C:\\Users\\Lenovo\\Desktop\\XML\\Test\\dirnc"
#
#nc = NetCDF.read( dir*"\\cartesian_fn.nc" )
#nc1 = NetCDF.read( dir*"\\geometry_tn.nc")
#
#plot(nc)
#plot!(nc1)



using NCDatasets
using ZipFile
using Plots


test = [ "C:\\Users\\DAVIDE-FAVARO\\Desktop\\XML\\1.xml",
         "C:\\Users\\DAVIDE-FAVARO\\Desktop\\XML\\2.xml",
         "C:\\Users\\Lenovo\\Desktop\\XML\\1.xml",
         "C:\\Users\\Lenovo\\Desktop\\XML\\2.xml",
         "C:\\Users\\Lenovo\\Desktop\\XML\\Prod_Test.xml" ]

out = [ "D:\\Vario\\Stage",
        "C:\\Users\\Lenovo\\Desktop\\XML",
        "C:\\Users\\Lenovo\\Desktop\\XML\\Test",
        "C:\\Users\\DAVIDE-FAVARO\\Desktop",
        "C:\\Users\\DAVIDE-FAVARO\\Desktop\\XML" ]


"""
    unzipIn( dir::AbstractString )

Unzip a ".zip" file stored in "dir", creating a new directory within "dir" containing the file contents 
"""
function unzip( dir::AbstractString )
    zipPath = "$dir\\$(readdir(dir)[1])"
    println(zipPath)
    zip = ZipFile.Reader(zipPath)
    new = mkdir(zipPath[1:end-4])
    newdir = zipPath[1:end-4]*"\\"
    start = length(zip.files[1].name) + 1
    for i in 2:length(zip.files)
        write( newdir*zip.files[i].name[start:end], read(zip.files[i]) )
    end
    return new
end


# Per ora ottiene le informazioni sui file .nc contenuti nello zip scaricato
function plotNC(  dir::AbstractString )
    cur = pwd()
    cd(dir)
    files = readdir()
    info = [ NCDataset( file ) for file in files[1:end-1] ]
    cd(cur)
    return info
end


dir = unzip( out[3] )
info = plotNC( dir )



#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#                                                            DOWNLOAD DEI FILE XML COLLEGATI AI PRODOTTI
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||


"""
    getProducts( authToken::AbstractString[, maxNumber::Union{Integer, Nothing} ] )

Obtain "maxNumber" pages, eachone containing the XML description of 100 products, through "authToken", returning the IOBuffer that contains them
"""
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
    val = tryparse( Int64, split( split( lines, "totalResults>" )[2], "<" )[1] )
    # "count" has two values: the first is the number of pages (undreds of products) and the second is a number of products smaller than 100 to be downloaded in an additional page
    count = [ val ÷ 100, val % 100 ]

# Check if maxNumber has a sensible value
    if !isnothing( maxNumber ) && maxNumber > 0 && maxNumber < val
        count[1] = maxNumber ÷ 100
        count[2] = maxNumber % 100
    end

# Download the desired number pages and sore the in an IOBuffer

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

# io = getProductsPages( authenticate("davidefavaro","Tirocinio"), 300 )



"""
    getPageProducts( fileIO::IO )

Given the IOBuffer obtained through "getProductsPages()", return an array of the dictionaries containing the informations on each of the products of the page 
"""
function getPageProducts( fileIO::IO )
# Get the downloaded pages containing the XML representations of the products
    original = replace( String( take!(fileIO) ), "\n" => "" )
# Split the downloaded files into pages
    pages = split( original, r"<\?xml [^<>]+><[^<>]+>", keepempty=false )
# Split the result in a vector of ready-to-be-parsed strings 
    vector = reduce( vcat, [ split( page, r"</?entry>", keepempty=false )[2:end-1] for page in pages ] )
# Parse the strings
    products = [ product(x)[8:end] for x in vector ]
# Generate a vector of dictionaries containing the details for each product of the original page
    return [ Dict( Symbol(join(p[:opening][7])) => parseConvert( join(p[:type][1]), join(p[:content]) ) for p in products[i] ) for i in 1:length(products) ]
end

# dict = getPageProducts(io)



"""
    getProductsDF( authToken::AbstractString, maxNumber::Union{Integer, Nothing} = nothing )

Generate the DataFrame containing the data of "maxNumber" products using "authToken", if "maxNumber" is not specified, the df will include all available products 
"""
function getProductsDF( authToken::AbstractString, maxNumber::Union{Integer, Nothing} = nothing )
# Download "maxNumber" pages and return the buffer containing them
    io = getProductsPages(  authToken, maxNumber )
# Create a vector of dictionaries of the products
    dict_vect = getPageProducts( io )
# Obtain the existing subsets of attributes of the products
    keys_groups = unique( keys.(dict_vect) )
# Divide the dictionaries in groups homogeneus on their attributes
    grouped_vect = [ filter( x -> keys(x) == ks, dict_vect ) for ks in keys_groups ]
# Turn each group in a DataFrame 
    dfs_vect = DataFrame.(grouped_vect)
# Merge all Dataframes using "DataFrames.append" to create a Dataframe with the union of all the possible columns and the right "missing" values 
    data = dfs_vect[1]
    for df in dfs_vect[2:end]
        append!( data, df, cols=:union )
    end
    return data
end

df = getProductsDF( authenticate("davidefavaro","Tirocinio"), 20 )                   



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