module Lakes



using ArchGDAL
using ArgParse
using Dates
using Sys


include("../Library/Functions.jl")


export run_lake



const agd = ArchGDAL



mutable struct Lake
    concentration::Float64
    time::Time
    distance_x::Float64
    distance_y::Float64
    fickian_x::Float64
    fickian_y::Float64
    velocity_x::Float64
    velocity_y::Float64
    wind_direction::Float64
    λk::Float64

    C::Float64

    Lake(concentration, time, distance_x, distance_y, fickian_x, fickian_y, velocity_x, velocity_y, wiind_direction, λk) = new(concentration, time, distance_x, distance_y, fickian_x, fickian_y, velocity_x, velocity_y, wiind_direction, λk)
end


function calc_concentration!( l::Lake )
 #=
     c1_1 = l.distance_x - ( l.velocity_x * l.time )
     c1_2 = c1_1^2
     c1_3 = c1_2 / ( 4 * l.fickian_x * l.time )
   
     c2_1 = l.distance_y - ( l.velocity_y * l.time )
     c2_2 = c2_1^2
     c2_3 = c2_2 / ( 4 * l.fickian_y * l.time )
   
     c3_2 = exp( -(c1_3 + c2_3) )
     c3_1 = exp( -l.λk * l.time )
   
     c4 = c3_2 * c3_1
   
     c5 = l.ma / ( 4π * l.time * √(l.fickian_x * l.fickian_y) ) 
     l.C = c4 * c5
 =#
    c1, c2 = ( (l.distance_x, l.distance_y) - ( (l.velocity_x, l.velocity_y) * l.time ) )^2 / ( 4 * (l.fickian_x, l.fickian_y) * l.time )
    c3 = ℯ^( -(c1 + c2) - (l.λk * l.time) )
    c4 = l.concentration / ( 4π * l.time * √(l.fickian_x * l.fickian_y) )
    l.C = c4 * c3
    return l.C
end



"""
    compute_result!( dtm::AbstractArray, r0::Integer, c0::Integer, ri::Integer, ci::Integer, lake::Lake )

Given the raster `dtm` and the indexes (`r0`, `c0`) of the source, modify the postion values of object `lake` and return the concentration at indexes (`ri`, `ci`)
"""
function compute_result!( dtm::AbstractArray, r0::Integer, c0::Integer, ri::Integer, ci::Integer, lake::Lake )
    lake.distance_x, lake.distance_y = Functions.compute_position(dtm, r0, c0, ri, ci, lake.wind_direction)
    return calc_concentration!(lake)
end



"""
    run_lake( source, wind_direction, pollutant_mass, flow_mean_speed, resolution::In64, hours::Int64,
              fickian_x::Real=0.05, fickian_y::Real=0.05, λk::Real=0.0, output_path::AbstractString=".\\lake_otput_model.tiff" )


#Arguments
- `source`: source point of the contaminants.
- `wind_direction`: direction of the wind as an angle in degrees.
- `pollutant_mass`: initial mass of contaminants.
- `flow_mean_speed`:  mean flow speed of the water.
- `resolution::In64`: size of the cell for the analysis.
- `hours::Int64`: time span of the analysis in hours.
- `fickian_x::Real=0.05`: X
- `fickian_y::Real=0.05`: X
- `λk::Real=0.0`: X
- `output_path::AbstractString=".\\lake_otput_model.tiff"`: output file path.
"""
function run_lake( source, wind_direction, pollutant_mass, flow_mean_speed, resolution::In64, hours::Int64,
                   fickian_x::Real=0.05, fickian_y::Real=0.05, λk::Real=0.0, output_path::AbstractString=".\\lake_otput_model.tiff" )

    geom = agd.getgeom(collect(agd.getlayer(source, 0))[1])
    refsys = agd.toWKT(agd.getspatialref(geom))
    x_source = agd.getx(geom, 0)
    y_source = agd.gety(geom, 0)
    r_source, c_source = toIndexes(dem, x_source, y_source)

    hours *= 3600
    velocity_x = √( round( flow_mean_speed * cos(deg2rad(wind_direction)), digits=3 )^2 )
    velocity_y = √( round( flow_mean_speed * sin(deg2rad(wind_direction)), digits=3 )^2 )

    points = [(r_source, c_source)]
    values = [pollutant_mass]
    lake = Lake(pollutant_mass, hours, 0, 0, fickian_x, fickian_y, velocity_x, velocity_y, wind_direction, λk)
    expand!(points, values, dem,  r_source, c_source, lake)

    maxR = maximum( point -> point[1], points )
    minR = minimum( point -> point[1], points )
    maxC = maximum( point -> point[2], points )
    minC = minimum( point -> point[2], points )

 #= SENZA FUNZIONI
    rows = maxR - minR
    cols = maxC - minC
    minX, maxY = toCoords(dtm, minR, maxC)

    valNoData = -9999

    gtiff_driver = agd.getdriver("GTiff")
    target_ds = agd.create( output_path, gtiff_driver, rows, cols, 1, agd.GDAL.GDT_Float32 )
    agd.setgeotransform!( target_ds, [ minX, resolution, 0.0, maxY, 0.0, -resolution ] )
    agd.setproj!(target_ds, refsys)
    band1 = agd.getband(target_ds, 1)
    agd.setnodatavalue!( band1, convert(Float64, valNoData) )
    band = agd.read(band1)
    agd.fillraster!(band, valNoData)

    for (point, value) in zip(points, values)
        r, c = point - (minR, minC)
        band[r, c] = value
    end

    agd.write!( target_ds, band )
 =#

    geotransform = agd.getgeotransform(dem)
    geotransform[[1, 4]] .+= (minR - 1, maxC - 1) .* geotransform[[2, 6]]

    #   data = [ isnothing( findfirst(p -> p == (r, c), points) ) ? noData : values[findfirst(p -> p == (r, c), points)] for r in minR:maxR, c in minC:maxC ]
    data = fill(noDataValue, maxR-minR, maxC-minC)
    for r in minR:maxR, c in minC:maxC
        match = findfirst(p -> p == (r, c), points)
        if !isnothing(match)
            data[r-minR+1, c-minC+1] = values[match]
        end
    end
    Functions.writeRaster(data, agd.getdriver("GTiff"), geotransform, resolution, refsys, noData, output_path, false)
end



end # module












using ArchGDAL
const agd = ArchGDAL


using Rasters


mutable struct Thing
    val::Float64
    a::Float64
    b::Float64
    x::Float64
    y::Float64
    dir::Float64
    res::Float64

    Thing(val, a, b, x, y, dir) = new(val, a, b, x, y, dir, 0)
end



function toCoords( dtm, r::Integer, c::Integer )
    return dtm.dims[1][r], dtm.dims[2][c]
end


function compute_position( dtm, r0::Integer, c0::Integer, ri::Integer, ci::Integer, direction::Real )
    Δx, Δy = toCoords(dtm, r0, c0) .- toCoords(dtm, ri, ci)
    dir = deg2rad(direction)
    sindir, cosdir = sin(dir), cos(dir)
    return Δx * cosdir - Δy * sindir, Δx * sindir + Δy * cosdir
end


function calc_val!( t::Thing )
    t.res = √(abs(t.a * t.x)) + √(abs(t.b * t.y)) / t.val
    return t.res
end


function compute_result!(mat, r0::Int64, c0::Int64, ri::Int64, ci::Int64, object::Thing)
    object.x, object.y = compute_position(mat, r0, c0, ri, ci, object.dir)
    return calc_val!(object)
end


function expand!( condition::Function, positions::AbstractVector, results::AbstractVector, dtm, indx_x::Integer, indx_y::Integer, object )
    if indx_x < 1 || indx_x > size(dtm, 1) || indx_y < 1 || indx_y > size(dtm, 1)
        return nothing
    end
    if (indx_x, indx_y) in positions
        xs = [ indx_x, indx_x-1, indx_x ]
        ys = [ indx_y+1, indx_y, indx_y-1 ]
        expand!( condition, positions, results, dtm, indx_x+1, indx_y, object )
        expand!.( condition, Ref(positions), Ref(results), Ref(dtm), xs, ys, Ref(deepcopy(object)) )
        return nothing
    else
        result = compute_result!(dtm, positions[1]..., indx_x, indx_y, object)
        println(result)
        if condition(result)
            push!( positions, (ind_x, ind_y) )
            push!( results, result )
            xs = [ indx_x, indx_x-1, indx_x ]
            ys = [ indx_y+1, indx_y, indx_y-1 ]
            expand!( condition, positions, concentrations, dtm, indx_x+1, indx_y, object )
            expand!.( condition, Ref(positions), Ref(results), Ref(dtm), xs, ys, Ref(deppcopy(object)) )
        end
        return nothing
    end
end

#=
dtm = agd.read(split( @__DIR__ , "\\Porting\\")[1] * "\\Mappe\\DTM_wgs84.tiff")
band = agd.getband(dtm, 1)
mat = band[3951:3975, 5951:5975]
=#
dtm = Raster(split( @__DIR__ , "\\Porting\\")[1] * "\\Mappe\\DTM_wgs84.tiff")

src_r = 3963
src_c = 5963
thing = Thing(20, 4, 3, src_r, src_c, 40)
pos = [ (src_r, src_c) ]
res = [ calc_val!(thing) ]
expand!( x -> x > 20 && x < 200, pos, res, dtm, src_r, src_c, thing )