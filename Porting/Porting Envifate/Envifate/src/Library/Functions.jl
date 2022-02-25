module Functions
"""
Module containing auxiliary functions
"""



using ArchGDAL
using CombinedParsers
using CombinedParsers.Regexp
using Rasters



export substance_extract, texture_extract, air_extract, cn_extract, cn_list_extract, array2raster!, writeRaster!, applystyle,
       getindex, setindex!, convert, -,
       getOrigin, getCellDims, getSidesDistances, toCoords, toIndexes,
       compute_position, expand!



const agd = ArchGDAL



@syntax dims = Sequence( "Pixel Size = (", Numeric(Float64), ",", Numeric(Float64), ")" )
@syntax points = Sequence( re"[^(]+", "(", re" *", Numeric(Float64), ",", re" *", Numeric(Float64), re".+" )



Base.getindex( collection::Raster, x::Float64, y::Float64 ) = collection[X(Near(x)), Y(Near(y))][1]
Base.setindex!( collection::Raster, v, x::Float64, y::Float64, ) = collection[X(Near(x)), Y(Near(y))] .= v
Base.convert(::Type{Int64}, n::AbstractFloat) = round(Int64, n)
Base.:-( x::Tuple{Vararg{T1}}, y::Vector{T2} ) where {T1<:Number, T2<:Number} = length(x) == length(y) ? Tuple( xi-yi for (xi, yi) in zip(x, y) ) : throw(ArgumentError("`x` and `y` must have the same size"))


#=
function substance_extract( substance_id, fields, dbloc = "" )
    # estrazione valori sostanze
    db = sql.DB(dbloc*"substance.db")
    sql_fields = join( fields, "," )
    query_substance = sql.Stmt( db, "SELECT ? FROM substance WHERE id = ?AAA" )
    sql.bind!( query_substance, [ sql_fields, substance_id ] )
    results = dbi.execute( query_substance )
    res_fields = [ x for x in results ]
    return res_fields
end

function texture_extract( texture_name, fields, dbloc = "" )
    # estrazione valori sostanze
    db = sql.DB(dbloc*"substance.db")
    sql_fields = join( fields, "," )
    query_texture = sql.Stmt( db, "SELECT ? FROM texture WHERE nome LIKE ?" )
    sql.bind!( query_texture, [ sql_fields, texture_name ] )
    results = dbi.execute( query_texture ) 
    res_fields = [ x  for x in results ]
    return res_fields
end

function air_extract( stability_class, outdoor, dbloc::AbstractString=*( @__DIR__, "\\") )
    db = sql.DB(dbloc*"substance.db")
    query_texture = sql.Stmt( db, "SELECT sigmay1, sigmay2, sigmayexp, sigmaz1, sigmaz2, sigmazexp FROM air_stability WHERE class LIKE ?NNN AND outdoor LIKE ?NNN" )
    sql.bind!( query_texture, [ stability_class, outdoor ] )
    results = dbi.execute( query_texture )
    res_fields = [ x for x in results ]
    return res_fields
end

function cn_extract( cnl, soil, dbloc::AbstractString=*( @__DIR__, "\\") )
    db = sql.DB(dbloc*"substance.db")
    classecn = "cn_"*String(cnl)
    query_cn = sql.Stmt( db, "SELECT ? FROM cn WHERE id = ?AAA" )
    sql.bind!( query_cn, [ classecn, soil ] )
    results = dbi.execute(query_cn)
    res_fields = [ x for x in results ]
	return res_fields
end

function cn_list_extract( dbloc::AbstractString=*( @__DIR__, "\\") )
	db = sql.DB(dbloc*"substance.db")
    query_cn = sql.Stmt( db, "SELECT * FROM cn" )
    results = dbi.execute(query_cn)
 #    listaclc = Dict()
 #    for row in results
 #        lista_soil = [ x for x in row ]
 #        listaclc[ row[5] ] = lista_soil
 #    end
    listaclc = Dict( row[5] => [ x for x  in row ] for row in results )
    return listaclc
end
=#


#= SOSTITUITI DALLA NUOVA writeRaster
function array2raster!( newRasterfn, xmin, ymin, pixelWidth, pixelHeight, xsize, ysize, array )
	# vedi https://pcjericks.github.io/py-gdalogr-cookbook/raster_layers.html#create-raster-from-array
    # cols = array.shape[1]
    # rows = array.shape[0]
    cols = xsize
    rows = ysize
    originX = xmin
    originY = ymin

    driver = agd.getdriver("GTiff")
    outRaster = agd.create( newRasterfn, driver, cols, rows, 1, ArchGDAL.GDT_Byte )
    agd.setgeotransform!( outRaster, [ originX, pixelWidth, 0, originY, 0, pixelHeight ] )
    outband = agd.getband(1)
    agd.write!( outband, array )
    #outRasterSRS = agd.importEPSG(4326)
    #agd.setproj!( outRaster,  agd.toWKT(outRasterSRS) )
 """
     outband.FlushCache()
 """
end

function writeRaster!( newRasterfn, xmin, ymin, pixelWidth, pixelHeight, xsize, ysize, array )
    reversed_arr = reverse(array) # reverse array so the tif looks like the array
    array2raster!( newRasterfn, xmin, ymin, pixelWidth, pixelHeight, xsize, ysize, reversed_arr ) # convert array to raster
end
=#
"""
    writeRaster( data::Array{Float32}, driver::agd.Driver, geotransform::Vector{Float64}, refsys::AbstractString, noData::Real, output_path::AbstractString=".\\raster.tiff", output::Bool=false )

Given a NxMxH dimensional matrix `data`, create a raster file with H NxM bands as `output_path` file, with `refsys` and `geotransfrom` as spatial references,
using `driver` to define the format.
If `output` is set to true return the new raster, otherwise return nothing.  
"""
function writeRaster( data::Array{Float32}, driver::ArchGDAL.Driver, geotransform::Vector{Float64}, resolution::Real, refsys::AbstractString, noDataValue::Real, output_file_path::AbstractString=".\\raster.tiff", output::Bool=false )
    rows, cols, bands = length(size(data)) < 3 ? (size(data)..., 1) : size(data) 
    res_raster = agd.create(output_file_path, driver=driver, width=rows, height=cols, nbands=bands, dtype=Float32)
 #= NON SO QUALE SIA IL COMANDO PER SETTARE I METADATI CON `ArchGDAL`
    target_ds.SetMetadata(
        Dict(
            "credits" => "Envifate - Francesco Geri, Oscar Cainelli, Paolo Zatelli, Gianluca Salogni, Marco Ciolli - DICAM Università degli Studi di Trento - Regione Veneto",
            "modulo" => "Analisi sedimentazione marina",
            "descrizione" => "Simulazione di sedimento disperso in ambiente marino",
            "srs" => refsys,
            "data" => today()
        )
    )
 =#
    for i in 1:bands
        agd.setnodatavalue!(agd.getband(res_raster, i), noDataValue)
        agd.write!(res_raster, data[:, :, i], i)
    end
    agd.setgeotransform!(res_raster, geotransform)
    agd.setproj!(res_raster, refsys)
    # NON SO QUANTO SIA NECESSARIO QUESTO IF
    if !output
        res_raster = nothing
        GC.gc()
    end
    return res_raster
end



#=
function applystyle( layer, colore, opacity )
    return nothing
end
=#


# Additional functions

"""
    getCellDims( dtm::ArchGDAL.AbstractDataset )

Return the cells' dimentions in meters for an ArchGDAL raster dataset
"""
function getCellDims( dtm::ArchGDAL.AbstractDataset )
    info = split( agd.gdalinfo( dtm ), "\n" , keepempty=false )
    pos = findfirst(occursin.("Pixel Size", info))
    size_pars = dims(info[pos])
    return size_pars[2], size_pars[4]
end



"""
    getOrigin( dtm::ArchGDAL.AbstractDataset )

Return the coordinates of the origin point of the raster
"""
function getOrigin( dtm::ArchGDAL.AbstractDataset )
    info = split( agd.gdalinfo( dtm ), "\n" , keepempty=false )
    pos = findfirst(occursin.("Origin", info))
    origin_pars = points(info[pos])
    return origin_pars[4], origin_pars[7]
end



"""
    getSidesDistances( dtm::ArchGDAL.AbstractDataset )

Return distances from left, upper, right and lower sides of an ArchGDAL raster dataset
"""
function getSidesDistances( dtm::ArchGDAL.AbstractDataset )
    info = split( agd.gdalinfo( dtm ), "\n" , keepempty=false )
    pos = findfirst(occursin.("Upper Left", info))
    dists_pars = points.( getindex.(Ref(info), [pos, pos+3] ) )
    return dists_pars[1][4], dists_pars[1][7], dists_pars[2][4], dists_pars[2][7]
end



"""
    toCoords( dtm::ArchGDAL.AbstractDataset, r::Integer, c::Integer )

Convert convert the indexes `r` and `c` to the coordinates of the respective cell in `dtm` raster
"""
function toCoords( dtm::ArchGDAL.AbstractDataset, r::Int64, c::Int64 )
    gtf = agd.getgeotransform(dtm)
    return gtf[[1,4]] .+ ( (r + 1/2) .* gtf[[2,5]] ) .+ ( (c + 1/2) .* gtf[[3,6]] )
end

function toCoords( geotransform::Vector{Float64}, r::Int64, c::Int64 )
    return geotransform[[1,4]] .+ ( (r + 1/2) .* geotransform[[2,5]] ) .+ ( (c + 1/2) .* geotransform[[3,6]] )
end

function toCoords( dtm::Rasters.Raster{T}, r::Int64, c::Int64 ) where {T}
    return dtm.dims[1][r], dtm.dims[2][c]
end


"""
    toIndexes( dtm::ArchGDAL.AbstractDataset, x::Real, y::Real )

Convert coordinates to the indexes of the respective cell in `dtm` raster
"""
function toIndexes( dtm::ArchGDAL.AbstractDataset, x::Real, y::Real )
    gtf = agd.getgeotransform(dtm)
    return round.( Int64, ( ( (x, y) .- gtf[[1,4]] ) ./ gtf[[2,6]] ) )
end

function toIndexes( geotransform::Vector{Float64}, x::Real, y::Real )
    return round.( Int64, ( ( (x, y) .- geotransform[[1,4]] ) ./ geotransform[[2,6]] ) )
end



"""
    compute_position( r0::Integer, c0::Integer, ri::Integer, ci::Integer, direction::Real )

Given the indexes of the source cell (`r0`, `c0`), those of the current one (`ri`, `ci`) and the angular direction of flow, compute the (`x`, `y`) coordinates
"""
function compute_position( dtm::AbstractArray, r0::Integer, c0::Integer, ri::Integer, ci::Integer, direction::Real )
    Δx, Δy = toCoords(dtm, r0, c0) - toCoords(dtm, ri, ci)
    dir = deg2rad(direction)
    sindir, cosdir = sin(dir), cos(dir)
    return Δx * cosdir - Δy * sindir, Δx * sindir + Δy * cosdir
end



"""
    expand!( condition::Function, positions::AbstractVector, results::AbstractVector, dtm::AbstractArray, indx_x::Integer, indx_y::Integer, object )::Nothing

Recursively compute the concentration of a substance spreading at cell (`indx_x`, `indx_y`) of `dtm` and in the adjacent cells,
adding all cells touched by the substance in `positions` and the relative concentration in `results` and accounting for the specificity of the
substance and the physical context through `object`.

The function, for each cell starting from the one at (`indx_x`, `indx_y`), checks whether it has been already visited, if so skips it, otherwise computes the concentration
of substance in it, keeping track of the total distance from the source through the first element of `positions`, which is always the source point.
If the concentration is sufficient according to `condition`, the function continues on the four cardinal adjacent cells.

The parameter `object` accounts for the properties of the substance and the enviroment surrounding it at each instant ad is modified every time the concentration is computed
allowing to keep track of the changes.

The function calls to method `compute_result!` of which exists a version specific for each module that employs `expand!`.
"""
function expand!( condition::Function, positions::AbstractVector, results::AbstractVector, dtm::AbstractArray, indx_x::Integer, indx_y::Integer, object )
    if indx_x < 1 || indx_x > size(dtm, 1) || indx_y < 1 || indx_y > size(dtm, 1)
        return nothing
    end
    if (indx_x, indx_y) in positions
        xs = [ indx_x, indx_x-1, indx_x ]
        ys = [ indx_y+1, indx_y, indx_y-1 ]
        expand!( condition, positions, results, dtm, indx_x+1, indx_y, object )
        expand!.( condition, Ref(positions), Ref(results), Ref(dtm), xs, ys, deepcopy(object) )
        return nothing
    else
        result = compute_result!(dtm, positions[1]..., indx_x, indx_y, object)
        if condition(result)
            push!( positions, (ind_x, ind_y) )
            push!( results, result )
            xs = [ indx_x, indx_x-1, indx_x ]
            ys = [ indx_y+1, indx_y, indx_y-1 ]
            expand!( condition, positions, results, dtm, indx_x+1, indx_y, object )
            expand!.( condition, Ref(positions), Ref(results), Ref(dtm), xs, ys, deppcopy(object) )
        end
        return nothing
    end
end







#= NON UTILIZZATE
function rasterFromPointsValues( points::Vector{Tuple{Int64, Int64}}, values::Vector{T}, driver::agd.Driver, resolution::Real, geotransform::Vector{Float64}, refsys::AbstractString, noDataValue::T, output_file_path::AbstractString ) where {T <: Real}
    if length(points) != length(values)
        throw(DomainError("There must be the same number of points in `points` and values in `values`"))
    end

    #  <sort points in increasing order>

    minR = minimum(x -> x[1], points) 
    minC = minimum(x -> x[2], points)
    maxR = maximum(x -> x[1], points)
    maxC = maximum(x -> x[2], points)

    #   data = [ isnothing( findfirst(p -> p == (r, c), points) ) ? noData : values[findfirst(p -> p == (r, c), points)] for r in minR:maxR, c in minC:maxC ]
    data = fill(noDataValue, maxR-minR, maxC-minC)
    for r in minR:maxR, c in minC:maxC
        match = findfirst(p -> p == (r, c), points)
        if !isnothing(match)
            data[r-minR+1, c-minC+1] = values[match]
        end
    end

    writeRaster(data, driver, geotransform, refsys, noDataValue, output_file_path, false )
end

function rasterFromPointsValues( points::Vector{Tuple{Int64, Int64}}, value::T, driver::agd.Driver, resolution::Real, geotransform::Vector{Float64}, refsys::AbstractString, noDataValue::T, output_file_path::AbstractString ) where {T <: Real}
    # <sort points in increasing order>
    
    minR = minimum(x -> x[1], points) 
    minC = minimum(x -> x[2], points)
    maxR = maximum(x -> x[1], points)
    maxC = maximum(x -> x[2], points)

    data = [ isnothing(findfirst(p -> p == (r, c), points)) ? noDataValue : value for r in minR:maxR, c in minC:maxC ]

    writeRaster( data, driver, geotransform, refsys, noDataValue, output_file_path, false )
end
=#




end # module
#=
""" Main problems with raster data handling:
- There seems to be no direct way to convert spatial measurements in degrees to meters.
- There is no one package that seems to be complete, the two more promising for the sake of the project seems to be `ArchGDAL.jl` and `Rasters.jl`
    but both seem to be lacking in some department that is central for the development of the project, while also being incompatible with eachother.
- Currently, in `ArchGDAL.jl`, there seems to be no way of obtaining the exact (X, Y) coordinates of a cell at indexes (row, column).
    It's possible to obtain the minimum X coordinate value, the maximum Y value and the value of the displacement between a cell
    and the next one on both axis as the geotransformation vector with:
        `ArchGDAL.getgeotransform( raster )`
    which returns a vector of six values including the aforementioned ones.
    There is also the function:
        `ArchGDAL.applygeotransform( geotransform, row, column )`
    that would theoretically serve the purpose, but in practice returns an inexact value, albeit not by much.
    The method seems also to be the main way of obtaining indexes from coordinates, using the inverse of the geotransform as input
    along with the X and Y values, obtaining the inverse transformationvector with:
        `ArchGDAL.invgeotransform!( geotransform, destinationvector )`
    having `geotransform` the transformation to invert and `destinationvector` the vector  that will be updated with the new values
    (it needs to be declared and intitialized beforehand as a vector of `Float64` values).
- The package `Rasters.jl`, for every `Raster` object created mantains two arrays of the coordinates corresponding to each cell,
    thus allowing to find the coordinates of a cell by accessing these arrays with the row and column of the cell.
    The problem can also be partially bypassed with the use of `Selectors` defined in the `DimensionalData.jl` package.
    These objects allow to index an Array through the use of coordinates, for example:
        `raster[Near(coordinateX), Near(coordinateY)]`
    this line accesses the raster `raster` at the cell with coordinates values closest to (`coordinateX`, `coordinateY`).
    While solving the problem of `ArchGDAL` the package lacks the same functionalities, not being able to operate directly or
    even read `Shapefiles` (however, it does support the use of shapefiles, read through the `Shapefile.jl` package, as input for
    some of its methods).
"""
=#
#= TESTING 


import ArchGDAL as agd
using Rasters

# Testing conversione coordinate -> indici e indici -> coordinate


dtm_file = split( @__DIR__ , "\\Porting\\")[1] * "\\Mappe\\DTM_wgs84.tiff"
dtm = agd.read(dtm_file)
rdtm = Raster(dtm_file)


sat = agd.read("D:\\Documents and Settings\\DAVIDE-FAVARO\\My Documents\\GitHub\\Tirocinio\\Mappe\\sat WGS84\\sette_sorelle.shp")
feature = collect(agd.getlayer(sat, 0))[1]
poly = agd.getgeom(feature, 0)
line = agd.getgeom(poly, 0)
point = agd.getpoint(line, 0)
x, y = point[1:2]

gtf = agd.getgeotransform(dtm)
igtf = deepcopy(gtf)
agd.invgeotransform!(gtf, igtf)

r1, c1 = agd.applygeotransform(igtf, x, y)
x1, y1 = agd.applygeotransform(gtf, r1, c1)



r2 = round( Int64, ((x - gtf[1]) / gtf[2]) )
c2 = round( Int64, ((y - gtf[4]) / gtf[6]) )
x2 = gtf[1] + ((r2 - 1) * gtf[2])
y2 = gtf[4] + ((c2 - 1) * gtf[6])


x3 = gtf[1] + ( (r2 + 1/2) * gtf[2] ) + ( (c2 + 1/2) * gtf[3] )
y3 = gtf[4] + ( (r2 + 1/2) * gtf[5] ) + ( (c2 + 1/2) * gtf[6] )


x, y = (10.5998, 44.7666)





# Test creazione raster


dtm_file = split( @__DIR__ , "\\Porting\\")[1] * "\\Mappe\\DTM_wgs84.tiff"
dtm = agd.read(dtm_file)

points = Vector{Tuple{Int64, Int64}}()
values = Vector{Float32}()
for r in 1000:1010, c in 1000:1010
    if r == 1010 && c == 1002
        break
    end
    push!( points, (r, c) )
    push!( values, 156.f0 + (r-1000)*(c-1000) )
end

maxR = maximum( point -> point[1], points )
minR = minimum( point -> point[1], points )
maxC = maximum( point -> point[2], points )
minC = minimum( point -> point[2], points )
rows = maxR - minR + 1
cols = maxC - minC + 1

gtf = agd.getgeotransform(dtm)
gtf[1] += (minR - 1) * gtf[2]
gtf[4] += (maxC - 1) * gtf[6]

rfs = agd.getproj(dtm)

noData = -9999.f0

data = [ isnothing( findfirst(p -> p == (r, c), points) ) ? noData : values[findfirst(p -> p == (r, c), points)] for r in minR:maxR, c in minC:maxC ]

writeRaster( data, agd.getdriver("GTiff"), gtf, rfs, noData, "C:\\Users\\DAVIDE-FAVARO\\Desktop\\test.tiff", false )








target_ds = agd.create( "C:\\Users\\DAVIDE-FAVARO\\Desktop\\test.tiff", driver=agd.getdriver("GTiff"), width=rows, height=cols, nbands=1, dtype=Float32)

agd.setnodatavalue!(agd.getband(target_ds, 1), noData)

agd.write!(target_ds, data, 1)

agd.setgeotransform!(target_ds, gtf)
agd.setproj!(target_ds, rfs)

=#