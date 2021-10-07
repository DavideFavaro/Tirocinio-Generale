module Test

using HTTP
using Downloads

using CombinedParsers
using CombinedParsers.Regexp

using ArchGDAL
using DataFrames
using Dates
using CSV

using Revise



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


# Syntax to parse the XMl reppresentation of a product
@syntax product = Repeat(
                    Either(
                        Sequence( "<",
                                  :type    => re"[^< >]+ ",
                                  :opening => re"name=\"[^<\">]+\">",
                                  :content => re"[^<>]+",
                                  :closing => re"</[^<>]+>"
                        ),
                        Sequence(
                            re"<[^<>]+>",
                            re"[^<>]+",
                            re"</[^<>]+>"
                        ),
                        re"<[^<>]+>"
                    )
                  )









# Funzione presa da Stackoverflow https://stackoverflow.com/questions/48104390/julias-most-efficient-way-to-choose-longest-array-in-array-of-arrays
function maxLenIndex( vect::AbstractVector )
    len = 0
    index = 0
    @inbounds for i in 1:length( vect )
        l = length( vect[i] )
        l > len && ( index = i; len=l )
    end
    return index
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



#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#                                                            DOWNLOAD DEI FILE XML COLLEGATI AI PRODOTTI
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||


"""
    getProductsBuffer( authToken::AbstractString[, maxNumber::Union{Integer, Nothing} ] )

Obtain "maxNumber" XML description of products, through "authToken", returning the IOBuffer that contains them all
"""
function getProductsBuffer( authToken::AbstractString, maxNumber::Union{Integer, Nothing} = nothing )
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

# io = getProductsBuffer( authenticate("davidefavaro","Tirocinio"), 300 )














"""
    getProductsDicts( fileIO::IO )

Given the IOBuffer obtained through "getProductsBuffer()", return an array of the dictionaries containing the informations on each of the products of the page 
"""
function getProductsDicts( fileIO::IO )
# Get the downloaded XML representations of the products
    original = replace( String( take!(fileIO) ), "\n" => "" )
# Split the downloaded files into pages
    pages = split( original, r"<\?xml [^<>]+><[^<>]+>", keepempty=false )
# Split the result in a vector of ready-to-be-parsed strings representing single products
    vector = reduce( vcat, [ split( page, r"</?entry>", keepempty=false )[2:end-1] for page in pages ] )
# Parse the strings
    products = [ product(x)[8:end] for x in vector ]
# Generate a vector of dictionaries containing the details for each product of the original page adding to each of them an additional value to account for the
    # original order of the data
    return [setindex!(
                Dict( Symbol( join( prod[:opening][7] ) ) => parseConvert( join( prod[:type][1] ), join( prod[:content] ) ) for prod in products[i] ),
                i,
                :orderNumber
            )
            for i in 1:length(products)]
end

#   dict = getProductsDicts(io)

#io = getProductsBuffer( authenticate("davidefavaro","Tirocinio"), 300 )
#res = getProductsDicts(io)





"""
    getProductsDF( authToken::AbstractString, maxNumber::Union{Integer, Nothing} = nothing )

Generate the DataFrame containing the data of "maxNumber" products using "authToken", if "maxNumber" is not specified, the df will include all available products 
"""
function getProductsDF( authToken::AbstractString, maxNumber::Union{Integer, Nothing} = nothing )
# Download "maxNumber" pages and return the buffer containing them
    io = getProductsBuffer(  authToken, maxNumber )
# Create a vector of dictionaries of the products
    dict_vect = getProductsDicts( io )
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

# Convert `footprint` e `gmlfootprint` columns in geometries
    data[!, :footprint] = ArchGDAL.fromWKT.( data[:, :footprint] )
    data[!, :gmlfootprint] = ArchGDAL.fromGML.( replace.( replace.( data[:, :gmlfootprint], "&lt;" => "<" ), "&gt;" => ">" ) )

# Order the rows based on `orderNumber`, then remove said column
    sort!(data, [:orderNumber] )
    select!( data, Not(:orderNumber) )

# Create the column indicating wether a product has been downloaded 
    insertcols!( data, ( :available => fill( true, nrow(data) ) ) ),( :downloaded => fill( false, nrow(data) ) )


    return data
end

















"""
    saveProductsDF( targetDirectory::AbstractString, data::DataFrame; overwrite::Bool=false )

Save "data" in "targetDirectory" if not already existing or if "overwrite" is true, otherwise append its content to "data.csv"
"""
function saveProductsDF( targetDirectory::AbstractString, data::DataFrame; overwrite::Bool=false )
    condition = !overwrite && in("data.csv", readdir(targetDirectory))
    if condition
        # Preexisting data
        old_data = CSV.read( targetDirectory*"\\data.csv", DataFrame )
        
# DA TESTARE QUANDO CI SARANNO EFFETTIVAMENTE DATI NUOVI
        # Index of the last element of the preexisting data
            # It will be used to find the duplicated data in "data"
        old_last_id = odata[end, :uuid]
        #Index of the first element of the new data
            # It will be used to find the now-unavailable products in "old_data"
        new_first_id = data[1, :uuid]

        # Index of the last duplicate in "data"
        duplicated_index = findfirst( ==(old_last_id), data[:, :uuid] )
        # Index of the last available product in "old_data"
        unavailable_index = findfirst( ==(new_first_id), old_data )

        #filter!( prod -> rownumber(prod) > duplicated_index, data )
        data = data[1:duplicated_index, :]
        old_data[ 1:unavailable_index, :available ] = false
    end
    CSV.write( targetDirectory*"\\data.csv", data, append = condition )
end



"""
    createUnitsDF( columns::AbstractVector{AbstractString}, units::AbstractVector{AbstractString} )

Given the array of column names and the array of their respective units of measure, create a df associating both
"""
function createUnitsDF( columns::AbstractVector{AbstractString}, units::AbstractVector{AbstractString} )
    return Dataframe( Dict( columns .=> units ) )
end



"""
    setDownloaded( fileIDs::Union{ AbstractString, AbstractVector{AbstractString} } )

Mark all products in "fileIDs" as already downloaded
"""
function setDownloaded( fileIDs::Union{ AbstractString, AbstractVector{AbstractString} } )
    updateProductsVal( fileIDs, "products", "downloaded", true )
end



"""
    setUnavailable( fileID::Union{ AbstractString, AbstractVector{AbstractString} } )

Mark all products in "fileIDs" as unavailable
"""
function setUnavailable( fileIDs::Union{ AbstractString, AbstractVector{AbstractString} } )
    updateProductsVal( fileIDs, "products", "available", false )
end





#df = getProductsDF( authenticate("davidefavaro","Tirocinio"), 1000 )
#saveProductsDF( out[2], df )



"C:\\Users\\Lenovo\\Documents\\GitHub\\Tirocinio\\Dati di prova"

df = getProductsDF( authenticate("davidefavaro","Tirocinio") )
saveProductsDF( "C:\\Users\\Lenovo\\Documents\\GitHub\\Tirocinio\\Dati di prova", df )





end #module