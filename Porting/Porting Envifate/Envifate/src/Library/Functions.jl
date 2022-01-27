module Functions
"""
Module containing auxiliary functions
"""



import ArchGDAL as agd


using CombinedParsers
using CombinedParsers.Regexp
#   import DBInterface as dbi
#   import SQLite as sql



export substance_extract, texture_extract, air_extract, cn_extract, cn_list_extract, array2raster!, writeRaster!, applystyle,
       +, -, *, /, ^,
       getCellDims, getSidesDistances, toCoords, toIndexes


@syntax dims = Sequence( "Pixel Size = (", Numeric(Float64), ",", Numeric(Float64), ")" )
@syntax points = Sequence( re"[^(]+", "(", re" *", Numeric(Float64), ",", re" *", Numeric(Float64), re".+" )



#= Old version
    Base.convert(::Type{Int64}, n::Float64) = round(Int64, n)
    Base.:-( x::Tuple{Number, Number}, y::Tuple{Number, Number} ) = ( x[1] - y[1], x[2] - y[2] )
    Base.:-( x::Vector{T}, y::Tuple{T, T} ) where {T <: Number} = length(x) == length(y) ? [ e1 - e2 for (e1, e2) in zip(x, y) ] : throw(ArgumentError("`x` and `y` must have the same size"))
    Base.:-( x::Tuple{T, T}, y::Vector{T} ) where {T <: Number} = length(x) == length(y) ? Tuple( e1 - e2 for (e1, e2) in zip(x, y) ) : throw(ArgumentError("`x` and `y` must have the same size"))
    Base.:+( x::Tuple{Number, Number}, y::Tuple{Number, Number} ) = ( x[1] + y[1], x[2] + y[2] )
    Base.:*( x::Tuple{Number, Number}, y::Number ) = ( x[1] * y, x[2] * y )
    Base.:*( x::Number, y::Tuple{Number, Number} ) = y * x
    Base.:*( x::Tuple{Number, Number}, y::Tuple{Number, Number} ) = ( x[1] * y[1], y[1] * y[2] )
    Base.:/( x::Tuple{Number, Number}, y::Number ) = ( x[1] / y, x[2] / y )
    Base.:/( x::Number, y::Tuple{Number, Number} ) = y / x
    Base.:/( x::Tuple{Number, Number}, y::Tuple{Number, Number} ) = ( x[1] / y[1], x[2] / y[2] )
    Base.:^( x::Tuple{Number, Number}, y::Number ) = ( x[1]^y, x[2]^y )
=#
#= Non sono necessarie basta usare il broadcasting
    Base.:+( x::Tuple{T1, T2}, y::Tuple{T1, T2} ) where {T1<:Number, T2<:Number} = length(x) == length(y) ? Tuple( xi+yi for (xi, yi) in zip(x,y) ) : throw(ArgumentError("`x` and `y` must have the same size"))
    Base.:+( x::Tuple{Vararg{T1}}, y::Tuple{Vararg{T2}} ) where {T1<:Number, T2<:Number} = length(x) == length(y) ? Tuple( xi+yi for (xi, yi) in zip(x,y) ) : throw(ArgumentError("`x` and `y` must have the same size"))
    Base.:-( x::Tuple{T1, T2}, y::Tuple{T1, T2} ) where {T1<:Number, T2<:Number} = length(x) == length(y) ? Tuple( xi-yi for (xi, yi) in zip(x,y) ) : throw(ArgumentError("`x` and `y` must have the same size"))
    Base.:-( x::Tuple{Vararg{T1}}, y::Tuple{Vararg{T2}} ) where {T1<:Number, T2<:Number} = length(x) == length(y) ? Tuple( xi-yi for (xi, yi) in zip(x,y) ) : throw(ArgumentError("`x` and `y` must have the same size"))
    Base.:-( x::Vector{T1}, y::Tuple{Vararg{T2}} ) where {T1<:Number, T2<:Number} = length(x) == length(y) ? [ xi-yi for (xi, yi) in zip(x, y) ] : throw(ArgumentError("`x` and `y` must have the same size"))
    Base.:*( x::Tuple{Vararg{T1}}, y::Tuple{Vararg{T2}} ) where {T1<:Number, T2<:Number} = length(x) == length(y) ? Tuple( xi*yi for (xi, yi) in zip(x,y) ) : throw(ArgumentError("`x` and `y` must have the same size"))
    Base.:*( x::Tuple{Vararg{T1}}, y::T2 ) where {T1<:Number, T2<:Number} = Tuple( xi*y for xi in x )
    Base.:*( x::T1, y::Tuple{Vararg{T2}} ) where {T1<:Number, T2<:Number} = y * x
    Base.:/( x::Tuple{Vararg{T1}}, y::Tuple{Vararg{T2}} ) where {T1<:Number, T2<:Number} = length(x) == length(y) ? Tuple( xi/yi for (xi, yi) in zip(x,y) ) : throw(ArgumentError("`x` and `y` must have the same size"))
    Base.:/( x::Tuple{Vararg{T1}}, y::T2 ) where {T1<:Number, T2<:Number} = Tuple( xi/y for xi in x )
    Base.:/( x::T1, y::Tuple{Vararg{T2}} ) where {T1<:Number, T2<:Number} = Tuple( x/yi for yi in y )
    Base.:^( x::Tuple{Vararg{T1}}, y::T2 ) where {T1<:Number, T2<:Number} = Tuple( xi^y for xi in x )
=#
"""
Convert a floating point number to integer by rounding it
"""
Base.convert(::Type{Int64}, n::Real) = round(Int64, n)
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

Given a N dimensional matrix `data`, create a raster file as `output_path` with `refsys` as spatial reference and `geotransfrom``, using `driver` to define the format.
If `output` is set to true return the new raster.  
"""
function writeRaster( data::Array{Float32}, driver::agd.Driver, geotransform::Vector{Float64}, refsys::AbstractString, noData::Real, output_path::AbstractString=".\\raster.tiff", output::Bool=false )
    rows, cols, bands = length(size(data)) < 3 ? (size(data)..., 1) : size(data) 
    res_raster = agd.create(output_path, driver=driver, width=rows, height=cols, nbands=bands, dtype=Float32)
    for i in 1:bands
        agd.setnodatavalue!(agd.getband(res_raster, i), noData)
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
    getCellDims( dtm::AbstractDataset )

Return the cells' dimentions in meters for an ArchGDAL raster dataset
"""
function getCellDims( dtm::AbstractDataset )
    info = split( agd.gdalinfo( dtm ), "\n" , keepempty=false )
    pos = findfirst(occursin.("Pixel Size", info))
    size_pars = dims(info[pos])
    return size_pars[2], size_pars[4]
end



"""
    getOrigin( dtm::AbstractDataset )

Return the coordinates of the origin point of the raster
"""
function getOrigin( dtm::AbstractDataset )
    info = split( agd.gdalinfo( dtm ), "\n" , keepempty=false )
    pos = findfirst(occursin.("Origin", info))
    origin_pars = points(info[pos])
    return origin_pars[4], origin_pars[7]
end



"""
    getSidesDistances( dtm::AbstractDataset )

Return distances from left, upper, right and lower sides of an ArchGDAL raster dataset
"""
function getSidesDistances( dtm::AbstractDataset )
    info = split( agd.gdalinfo( dtm ), "\n" , keepempty=false )
    pos = findfirst(occursin.("Upper Left", info))
    dists_pars = points.( getindex.(Ref(info), [pos, pos+3] ) )
    return dists_pars[1][4], dists_pars[1][7], dists_pars[2][4], dists_pars[2][7]
end



"""
    toCoords( dtm::AbstractDataset, r::Integer, c::Integer )

Convert convert the indexes `r` and `c` to the coordinates of the respective cell in `dtm` raster
"""
function toCoords( dtm::AbstractDataset, r::Integer, c::Integer )
    Δx, Δy = getCellDims(dtm)
    left, up = getOrigin(dtm)
    x = r * Δx + left
    y = c * Δy + up
    return x, y
end



"""
    toIndexes( dtm::AbstractDataset, x::Real, y::Real )

Convert coordinates to the indexes of the respective cell in `dtm` raster
"""
function toIndexes( dtm::AbstractDataset, x::Real, y::Real )
    Δx, Δy = getCellDims(dtm)
    left, up = getSidesDistances(dtm)[1:2]
    r = round( Int64, ( x - left ) / Δx  )
    c = round( Int64, ( y - up ) / Δy )
    return r, c
end



"""
    compute_position( r0::Integer, c0::Integer, ri::Integer, ci::Integer, direction::Real )

Given the indexes of the source cell (`r0`, `c0`), those of the current one (`ri`, `ci`) and the angular direction of flow, compute the (`x`, `y`) coordinates
"""
function compute_position( dtm::AbstractArray, r0::Integer, c0::Integer, ri::Integer, ci::Integer, direction::Real )
    Δx, Δy = Functions.toCoords(dtm, r0, c0) - Functions.toCoords(dtm, ri, ci)
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
    if (indx_x, indx_y) in positions
        xs = [ indx_x, indx_x-1, indx_x ]
        ys = [ indx_y+1, indx_y, indx_y-1 ]
        expand!( condition, positions, concentrations, dtm, indx_x+1, indx_y, object )
        expand!.( condition, Ref(positions), Ref(concentrations), Ref(dtm), xs, ys, deepcopy(object) )
        return nothing
    else
        result = compute_result!(dtm, positions[1]..., ind_x, ind_y, object)
        if condition(result)
            push!( positions, (ind_x, ind_y) )
            push!( results, result )
            xs = [ indx_x, indx_x-1, indx_x ]
            ys = [ indx_y+1, indx_y, indx_y-1 ]
            expand!( condition, positions, concentrations, dtm, indx_x+1, indx_y, object )
            expand!.( condition, Ref(positions), Ref(results), Ref(dtm), xs, ys, deppcopy(object) )
        end
        return nothing
    end
end



end # module


















import ArchGDAL as agd



dtm_file = split( @__DIR__ , "\\Porting\\")[1] * "\\Mappe\\DTM_wgs84.tiff"
dtm = agd.readraster(dtm_file)

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

