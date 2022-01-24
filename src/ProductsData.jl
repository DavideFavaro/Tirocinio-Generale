module ProductsData
"""
Module for the Download and processing of descriptors of satellitar data from the Copernicus project databases
"""



using CombinedParsers
using CombinedParsers.Regexp
using CSV
using Dates
using DataFrames
using Downloads
using HTTP

import ArchGDAL as agd



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
        Sequence(
            "<",
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




# Funzione presa da: https://stackoverflow.com/questions/48104390/julias-most-efficient-way-to-choose-longest-array-in-array-of-arrays
function maxLenIndex( vect::AbstractVector )
    len = 0
    index = 0
    @inbounds for i in 1:length(vect)
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



"""
"""
function getAoi( path::AbstractString )::AbstractString
    data = agd.read(path)
    layer = agd.getlayer(data, 0)
    geometry = agd.getgeom(collect(layer)[1])
    src_crs = agd.getspatialref(geometry)
    trg_crs = agd.importEPSG(4326)
    agd.createcoordtrans(src_crs, trg_crs) do transform
        agd.transform!(geometry, transform)
    end 
    geom = agd.toWKT(geometry)
    return geom[1:7] * replace( geom[9:end], " " => "%20" )
end



#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#                                                            DOWNLOAD DEI FILE XML COLLEGATI AI PRODOTTI
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

"""
    getProductsBuffer( authToken::AbstractString[, start::Integer, maxNumber::Union{Integer, Nothing} ] )

Obtain `maxNumber` XML description of products, starting from `start`, through `authToken`, returning the IOBuffer that contains them all
"""
function getProductsBuffer( authToken::AbstractString; aoi::AbstractString="POLYGON((9.5000%2047.0000,%2014.0000%2047.0000,%2014.0000%2044.0000,%209.5000%2044.0000,%209.5000%2047.0000))",
                            numMonths::Integer=6, start::Integer=0, max::Union{Integer, Nothing} = nothing, last::Bool=false )
 # Definition of the components of The URL
    query2 = "[NOW-$(numMonths)MONTHS%20TO%20NOW]%20AND%20footprint:\"Intersects($aoi)\""
    query = "search?start=0&rows=0&q=ingestiondate:$query2"
    iob = IOBuffer()

 # Get number of total products
    #Download the firs page (This way we can get the number we are interested in and also the first page that we would download anyway)
    Downloads.download( "https://scihub.copernicus.eu/dhus/$query", iob, headers = [ "Authorization" => "Basic $authToken" ] )
    lines = String(take!(iob))
    val = tryparse( Int64, split( split( lines, "totalResults>" )[2], "<" )[1] )
    # "count" has two values: the first is the number of pages (undreds of products) and the second is a number of products smaller than 100 to be downloaded in an additional page
    count = val

    #Check if `start` has a sensible value
    if ( last && start == 0 ) || start > val
        start = val -1
    end

    # Check if `maxNumber` has a sensible value
    if !isnothing(max) && max > 0 && max < val
        count = max
    end

    # Calculate the value of count in terms of pages of 100 products each accounting for the starting point
    if start < count
        count -= start
    end
    count = [ count ÷ 100, count % 100 ]

 # Download the desired number pages and store the in an IOBuffer
    if count[1] > 0
        for i in 0:count[1]-1
            query = "search?start=$(start + ( i * 100 ))&rows=100&q=ingestiondate:$query2"
            Downloads.download( "https://scihub.copernicus.eu/dhus/$query", iob, headers = [ "Authorization" => "Basic $authToken" ] )
        end
    end
    if count[2] > 0
        query = "search?start=$(start + (count[1] * 100))&rows=$(count[2])&q=ingestiondate:$query2"
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
    pages = split(original, r"<\?xml [^<>]+><[^<>]+>", keepempty=false)

    # Split the result in a vector of ready-to-be-parsed strings representing single products
    vector = reduce( vcat, [ split( page, r"</?entry>", keepempty=false )[2:end-1] for page in pages ] )

    # Parse the strings
    products = [ product(x)[8:end] for x in vector ]

    # Generate a vector of dictionaries containing the details for each product of the original page adding to each of them an additional value to account for the
        # original order of the data
    return [ push!(
                Dict( Symbol( join( prod[:opening][7] ) ) => parseConvert( join( prod[:type][1] ), join( prod[:content] ) ) for prod in products[i] ),
                :orderNumber => i, 
            )
            for i in 1:length(products) ]
end

#   dict = getProductsDicts(io)

#   io = getProductsBuffer( authenticate("davidefavaro","Tirocinio"), 5115, 50 )
#   res = getProductsDicts(io)



"""
    getProductsDF( authToken::AbstractString[, start::Integer, maxNumber::Union{Integer, Nothing} ] )

Generate the DataFrame containing the data of `maxNumber` products using `authToken`, if `maxNumber` is not specified, the df will include all available products,
if `start` is specified the downloaded data will begin from that index
"""
function getProductsDF( authToken::AbstractString; aoi::AbstractString="", numMonths::Integer=6, start::Integer=0, max::Union{Integer, Nothing}=nothing, last::Bool=false )::DataFrame
 # Download "maxNumber" pages and return the buffer containing them
    io = getProductsBuffer(authToken, aoi=aoi, numMonths=numMonths, start=start, max=max, last=last)
 # Create a vector of dictionaries of the products
    dict_vect = getProductsDicts(io)
 # Obtain the existing subsets of attributes of the products
    keys_groups = unique( keys.(dict_vect) )
 # Divide the dictionaries in groups homogeneus on their attributes
    grouped_vect = [ filter( x -> keys(x) == ks, dict_vect ) for ks in keys_groups ]
 # Turn each group in a DataFrame 
    dfs_vect = DataFrame.(grouped_vect)
 # Merge all Dataframes using "DataFrames.append" to create a Dataframe with the union of all the possible columns and the right "missing" values 
    data = dfs_vect[1]
    for df in dfs_vect[2:end]
        append!(data, df, cols=:union)
    end

 # Convert `footprint` e `gmlfootprint` columns in geometries
    data[!, :footprint] = agd.fromWKT.( data[:, :footprint] )
    data[!, :gmlfootprint] = agd.fromGML.( replace.( replace.( data[:, :gmlfootprint], "&lt;" => "<" ), "&gt;" => ">" ) )

 # Order the rows based on `orderNumber`, then remove said column
    sort!( data, [:orderNumber] )
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

        # Index of the first entry of the preexisting data
            # It will be used to find the duplicates in "data"
        old_first_id = old_data[1, :uuid]

        # Index of the first duplicate in "data"
        first_dup_idx = findfirst( ==(old_first_id), data[:, :uuid] )

        # Remove all the duplicates from "data"
        filter!( prod -> rownumber(prod) > first_dup_idx, data )
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



function checkAvailabile()
    # Get the last available product as a DataFrame row
    data = getProductsDF( authenticate("davidefavaro","Tirocinio"), max=1, last=true  )

    # Load the preexisting data
    old_data = CSV.read( *( @__DIR__, "\\Dati di prova\\data.csv"), DataFrame )

    # Obtain the `uuid` of the last available product
    last_available_id = data[end, :uuid]

    # Obtain the index of the last available product of "old_data"
    last_available_idx = findfirst( ==(last_available_id), old_data[:, :uuid] )

    # Set all the products not of level 0, from the last available to the end, as unavailable 
    for i in last_available_idx:size(old_data)[1]
        if ismissing(old_data[i, :productlevel]) || old_data[i, :productlevel] != "L0"
            old_data[i, :available] = false
        end
    end

    CSV.write( *( @__DIR__, "\\Dati di prova\\data.csv" ), old_data )
end

#   df = getProductsDF( authenticate("davidefavaro","Tirocinio"), 1000 )
#   saveProductsDF( out[2], df )


#   df = getProductsDF( authenticate("davidefavaro","Tirocinio"), max=1, last=true )
#   odf = CSV.read( "D:\\Documents and Settings\\DAVIDE-FAVARO\\My Documents\\GitHub\\Tirocinio\\Dati di prova\\data.csv", DataFrame )
#   odf = CSV.read( "C:\\Users\\Lenovo\\Documents\\GitHub\\Tirocinio\\Dati di prova\\data.csv", DataFrame )


#   saveProductsDF( "C:\\Users\\Lenovo\\Documents\\GitHub\\Tirocinio\\Dati di prova", df )
#   saveProductsDF( "D:\\Documents and Settings\\DAVIDE-FAVARO\\My Documents\\GitHub\\Tirocinio\\Dati di prova", df )


#   odf = CSV.read( "C:\\Users\\Lenovo\\Documents\\GitHub\\Tirocinio\\Dati di prova\\data.csv", DataFrame )
#   checkAvailabile()
#   ndf = CSV.read( "C:\\Users\\Lenovo\\Documents\\GitHub\\Tirocinio\\Dati di prova\\data.csv", DataFrame )




# :footprint o :gmlfootprint ci dovrebbero dare il poligono
#= La selezione dei prodotti potrebbe avvenire inbase a:
    - :processinglevel / :productlevel ( livello 2A )
    - :instrumentname / :instrumentshortname
    - :size
    - :cloudcoverpercentage + :platformname ( cc minima per Sentinel-2A )
    - :status ?
    - :available
=#

#   sat_file = *( @__DIR__, "\\..\\Mappe\\sat\\sette_sorelle.shp" )
#   aoi = getAoi(sat_file)
#   df = getProductsDF( authenticate("davidefavaro","Tirocinio"), aoi=aoi )
#
# PRENDERE IL PRIMO PRODOTTO PER OGNI MESE POTREBBE NON ESSERE UNA BUON IDEA
#   idxs = [ findfirst( date -> month(date) == m, df[:, :beginposition] ) for m in month( now() - Month(6) ) : month( now() ) ]
#
#   res = df[ idxs, : ]
#
#   aoi_geom = agd.fromWKT( replace(aoi, "%20" => " ") )
#   
#   intersections = agd.intersection.()

function downloadProduct( authentication::AbstractString, aoi_path::AbstractString, num_per_month::Integer, from::Integer, to::Integer )
    aoi = getAoi(aoi_path)
 # DA AGGIUNGERE LA POSSIBILITA' DI SPECIFICARE IL MESE DA CUI INIZIARE A PRENDERE I PRODOTTI
    df::DataFrame = getProductsDF( authentication, aoi=aoi, numMonths=to )
    idxs = Vector{Int64}()
    condition( date, platform, clouds, level, m, sat ) = month(date) == m && platform == sat && (platform != "Sentinel-2" || ( !ismissing(clouds) && clouds < 30.0 )) && !ismissing(level) && level == "L2"
    for m in month.(collect( now()-Month(6):Month(1):now() ))
        ind::Int64 = 1
        for i in 1:num_per_month
            for sat in ["Sentinel-1", "Sentinel-2", "Sentinel-3"]
                first::Union{Int64, Nothing} = findfirst( row -> condition(row..., m, sat), eachrow( df[ ind:end, [:beginposition, :platformname, :cloudcoverpercentage, :productlevel] ] ) )
                if !isnothing(first)
                    ind = first 
                    push!(idxs, first)
                    ind += 1
                end
            end
        end
    end
    return df[idxs, :]
end

# SI OTTENGONO SOLO Sentinel-3 PERCHE' TUTTI GLI ALTRI PRODOTTI HANNO :productlevel missing
# PER MESE 6 SI HANNO SOLO Sentinel-3

#   dir = "C:\\Users\\DAVIDE-FAVARO\\Desktop"
#   dir = "D:\\Z_Tirocinio_Dati\\Copernicus Data"
#   authToken = authenticate("davidefavaro","Tirocinio")
#   sat_file = *( @__DIR__, "\\..\\Mappe\\sat\\sette_sorelle.shp" )
#   res = downloadProduct( authToken, sat_file, 1, 12, 6 )
#= 
    for id in res[:, :uuid]
        Downloads.download( "https://scihub.copernicus.eu/dhus/odata/v1/Products('$id')/\$value", dir*"\\$id.zip", headers = [ "Authorization" => "Basic $authToken" ] )
    end


    uuid = "de9494a3-fdd9-46ce-98cb-42a2632e8b87"
    Downloads.download( "https://scihub.copernicus.eu/dhus/odata/v1/Products('$uuid')/\$value", dir*"\\$uuid", headers = [ "Authorization" => "Basic $authToken" ] )
=#




# @code_warntype res = downloadProduct( authToken, sat_file, 1, 12, 6 )



using Rasters
using Shapefile

sat_file = *( @__DIR__, "\\..\\Mappe\\sat\\sette_sorelle.shp" )
ncdf_files = "C:\\Users\\DAVIDE-FAVARO\\Desktop\\Dati Copernicus\\1609_S3B_OL_1_EFR"
#   ncdf_files = "D:\\Z_Tirocinio_Dati\\Copernicus Data\\S3A_SL_2_LST____20191122T183837"
#   files = NamedTuple( Symbol(f) = ncdf_files * "\\" * f for f in readdir(ncdf_files)[1:end-1] )
files = ncdf_files * "\\" .* readdir(ncdf_files)[1:end-1]

#   data = [ Rasters.Raster(f) for f in files[1:end-1] ];
data = RasterSeries(files, :n)

#   shape = Shapefile.Handle(sat_file).shapes[1]
shape = Shapefile.Table(sat_file)
#   shp = Shapefile.Polygon( shape[1].MBR, shape[1].parts, shape[1].points )

poly = Shapefile.Polygon(shape.geometry[1].MBR, [1,2,3,4], shape.geometry[1].points[1:end-1])

mask(data[1]; to=poly)

crop(data, to=shape)

end #module




RasterSeries{
    Raster{
        Union{Missing, Float32},
        2,
        Tuple{
            Dim{
                :columns,
                DimensionalData.Dimensions.LookupArrays.NoLookup{
                    Base.OneTo{Int64}
                }
            },
            Dim{
                :rows,
                DimensionalData.Dimensions.LookupArrays.NoLookup{
                    Base.OneTo{Int64}
                }
            }
        },
        Tuple{},
        Rasters.FileArray{
            Rasters.NCDfile,
            Union{Missing, Float32},
            2,
            Symbol,
            DiskArrays.GridChunks{2},
            DiskArrays.Unchunked
        },
        Symbol,
        DimensionalData.Dimensions.LookupArrays.Metadata{
            Rasters.NCDfile,
            Dict{Symbol, Any}
        },
        Missing
    },
    1,
    Tuple{
        Dim{
            :n,
            DimensionalData.Dimensions.LookupArrays.NoLookup{
                Base.OneTo{Int64}
            }
        }
    },
    Tuple{},
    Vector{
        Raster{
            Union{Missing, Float32},
            2,
            Tuple{
                Dim{
                    :columns,
                    DimensionalData.Dimensions.LookupArrays.NoLookup{
                        Base.OneTo{Int64}
                    }
                },
                Dim{
                    :rows,
                    DimensionalData.Dimensions.LookupArrays.NoLookup{
                        Base.OneTo{Int64}
                        }
                    }
            },
            Tuple{},
            Rasters.FileArray{
                Rasters.NCDfile,
                Union{Missing, Float32},
                2,
                Symbol,
                DiskArrays.GridChunks{2},
                DiskArrays.Unchunked
            },
            Symbol,
            DimensionalData.Dimensions.LookupArrays.Metadata{
                Rasters.NCDfile,
                Dict{Symbol, Any}
            },
            Missing
        }
    }
}



