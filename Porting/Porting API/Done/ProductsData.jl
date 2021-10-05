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

    data[!, :footprint] = ArchGDAL.fromWKT.( data[:, :footprint] )
    data[!, :gmlfootprint] = ArchGDAL.fromGML.( replace.( replace.( data[:, :gmlfootprint], "&lt;" => "<" ), "&gt;" => ">" ) )

    insertcols!( )

    return data
end



"""
    saveProductsDF( targetDirectory::AbstractString, data::DataFrame; overwrite::Bool=false )

Save "data" in "targetDirectory" if not already existing or if "overwrite" is true, otherwise append its content to "data.csv"
"""
function saveProductsDF( targetDirectory::AbstractString, data::DataFrame; overwrite::Bool=false )
    CSV.write( targetDirectory*"\\data.csv", data, append = !overwrite && in("data.csv", readdir(targetDirectory)) )
end



"""
    createUnitsDF( columns::AbstractVector{AbstractString}, units::AbstractVector{AbstractString} )

Given the array of column names and the array of their respective units of measure, create a df associating both
"""
function createUnitsDF( columns::AbstractVector{AbstractString}, units::AbstractVector{AbstractString} )
    return Dataframe( Dict( columns .=> units ) )
end


df = getProductsDF( authenticate("davidefavaro","Tirocinio"), 1000 )
saveProductsDF( out[2], df )









df = CSV.read( split( @__DIR__, "Porting")[1] * "\\Dati di Prova\\data.csv", DataFrame )

df[!, :footprint] = ArchGDAL.fromWKT.( df[:, :footprint] )
df[!, :gmlfootprint] = ArchGDAL.fromGML.( replace.( replace.( df[:, :gmlfootprint], "&lt;" => "<" ), "&gt;" => ">" ) )





end #module



#   #=
#   gmls[i][:geometry][:content][:type]
#       contiene il poligono (Polygon)
#   
#   gmls[i][:geometry][:content][:ref]
#       contiene il riferimento (4362)
#   
#   gmls[i][:body][4]
#       contiene le coordinate come vettore di caratteri
#   =#
#   gmls = gml.(gmlfootprint)
#   
#   
#   joined = join( join.([ "&lt;gml:",
#                          gmls[1][:geometry][:content][:type],      # Polygon
#                          "&gt;",
#                          gmls[1][:body][4],                        # 'coordinates'
#                          "&lt;/gml:",
#                          gmls[1][:geometry][:content][:type],      # Polygon
#                          "&lt;"
#                        ]))




#   # Syntax to parse and decompose a GML format string
#   @syntax gml = Sequence(
#                       :geometry => Sequence(
#                                       :ltgml => "&lt;gml:",
#                                       :content => Sequence(
#                                                       :type => re"[^ &;]+",
#                                                       :name => re" [^ &;#]+#",
#                                                       :ref => Numeric(Int),
#                                                       :context => re"\" [^ &;]+"
#                                                   ),
#                                       :gt => "&gt;",
#                                       :spacing => re" *"
#                                   ),    
#                       :body => Repeat(
#                                   Either(
#                                       Sequence(
#                                           :ltgml => "&lt;gml:",
#                                           :content => re"[^ &;]+",
#                                           :gt => "&gt;",
#                                           :spacing => re" *"
#                                       ),
#                                       re"[0-9 .,-]+",
#                                       Sequence(
#                                           :ltgml => "&lt;/gml:",
#                                           :content => re"[^&;]+",
#                                           :gt => "&gt;",
#                                           :spacing => re" *"
#                                       ),
#                                   )
#                               )
#                     )