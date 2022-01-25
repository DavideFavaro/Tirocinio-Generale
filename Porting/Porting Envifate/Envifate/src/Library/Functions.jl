module Functions

export substance_extract, texture_extract, air_extract, cn_extract, cn_list_extract, array2raster!, writeRaster!, applystyle,
       +, -, *, /, ^,
       getCellDims, getSidesDistances, toCoords, toIndexes

import ArchGDAL as agd
using CombinedParsers
using CombinedParsers.Regexp
#   import DBInterface as dbi
#   import SQLite as sql



@syntax dims = Sequence( "Pixel Size = (", Numeric(Float64), ",", Numeric(Float64), ")" )
@syntax points = Sequence( re"[^(]+", "(", re" *", Numeric(Float64), ",", re" *", Numeric(Float64), re".+" )

#= Old version =#
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

#=
Base.convert(::Type{Int64}, n::Float64) = round(Int64, n)
Base.:+( x::Tuple{T1, T2}, y::Tuple{T1, T2} ) where {T1<:Number, T2<:Number} = length(x) == length(y) ? Tuple( xi+yi for (xi, yi) in zip(x,y) ) : throw(ArgumentError("`x` and `y` must have the same size"))
Base.:+( x::Tuple{Vararg{T1}}, y::Tuple{Vararg{T2}} ) where {T1<:Number, T2<:Number} = length(x) == length(y) ? Tuple( xi+yi for (xi, yi) in zip(x,y) ) : throw(ArgumentError("`x` and `y` must have the same size"))
Base.:+( x::Tuple{T1, T2}, y::Tuple{T1, T2} ) where {T1<:Number, T2<:Number} = length(x) == length(y) ? Tuple( xi-yi for (xi, yi) in zip(x,y) ) : throw(ArgumentError("`x` and `y` must have the same size"))
Base.:-( x::Tuple{Vararg{T1}}, y::Tuple{Vararg{T2}} ) where {T1<:Number, T2<:Number} = length(x) == length(y) ? Tuple( xi-yi for (xi, yi) in zip(x,y) ) : throw(ArgumentError("`x` and `y` must have the same size"))
Base.:-( x::Tuple{Vararg{T1}}, y::Vector{T2} ) where {T1<:Number, T2<:Number} = length(x) == length(y) ? Tuple( xi-yi for (xi, yi) in zip(x, y) ) : throw(ArgumentError("`x` and `y` must have the same size"))
Base.:-( x::Vector{T1}, y::Tuple{Vararg{T2}} ) where {T1<:Number, T2<:Number} = length(x) == length(y) ? [ xi-yi for (xi, yi) in zip(x, y) ] : throw(ArgumentError("`x` and `y` must have the same size"))
Base.:*( x::Tuple{Vararg{T1}}, y::Tuple{Vararg{T2}} ) where {T1<:Number, T2<:Number} = length(x) == length(y) ? Tuple( xi*yi for (xi, yi) in zip(x,y) ) : throw(ArgumentError("`x` and `y` must have the same size"))
Base.:*( x::Tuple{Vararg{T1}}, y::T2 ) where {T1<:Number, T2<:Number} = Tuple( xi*y for xi in x )
Base.:*( x::T1, y::Tuple{Vararg{T2}} ) where {T1<:Number, T2<:Number} = y * x
Base.:/( x::Tuple{Vararg{T1}}, y::Tuple{Vararg{T2}} ) where {T1<:Number, T2<:Number} = length(x) == length(y) ? Tuple( xi/yi for (xi, yi) in zip(x,y) ) : throw(ArgumentError("`x` and `y` must have the same size"))
Base.:/( x::Tuple{Vararg{T1}}, y::T2 ) where {T1<:Number, T2<:Number} = Tuple( xi/y for xi in x )
Base.:/( x::T1, y::Tuple{Vararg{T2}} ) where {T1<:Number, T2<:Number} = Tuple( x/yi for yi in y )
Base.:^( x::Tuple{Vararg{T1}}, y::T2 ) where {T1<:Number, T2<:Number} = Tuple( xi^y for xi in x )
=#

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


function applystyle( layer, colore, opacity )
    return nothing
end
=#


# Additional functions

"""
Return the dimentions of the cell of an ArchGDAL raster dataset
"""
function getCellDims( dtm )
    info = split( agd.gdalinfo( dtm ), "\n" , keepempty=false )
    pos = findfirst(occursin.("Pixel Size", info))
    size_pars = dims(info[pos])
    return size_pars[2], size_pars[4]
end


"""
Return distances from left upper right and lower sides of an ArchGDAL raster dataset
"""
function getOrigin( dtm )
    info = split( agd.gdalinfo( dtm ), "\n" , keepempty=false )
    pos = findfirst(occursin.("Origin", info))
    origin_pars = points(info[pos])
    return origin_pars[4], origin_pars[7]
end


"""
Return distances from left, upper, right and lower sides of an ArchGDAL raster dataset
"""
function getSidesDistances( dtm )
    info = split( agd.gdalinfo( dtm ), "\n" , keepempty=false )
    pos = findfirst(occursin.("Upper Left", info))
    dists_pars = points.( getindex.(Ref(info), [pos, pos+3] ) )
    return dists_pars[1][4], dists_pars[1][7], dists_pars[2][4], dists_pars[2][7]
end


"""
Convert indexes to the coordinates of the respective cell in `dtm` raster
"""
function toCoords( dtm, r::Integer, c::Integer )
    Δx, Δy = getCellDims(dtm)
    left, up = getOrigin(dtm)
    x = r * Δx + left
    y = c * Δy + up
    return x, y
end


"""
Convert coordinates to the indexes of the respective cell in `dtm` raster
"""
function toIndexes( dtm, x::Real, y::Real )
    Δx, Δy = getCellDims(dtm)
    left, up = getSidesDistances(dtm)[1:2]
    r = round( Int64, ( x - left ) / Δx  )
    c = round( Int64, ( y - up ) / Δy )
    return r, c
end



#= DA MODIFICARE PER RENDERLA GENERICA
"""
Recursively compute the concentration of each point and add the value and its indexes to positions
"""
function expand!( positions::AbstractVector, results::AbstractVector, dtm, indx_x::Integer, indx_y::Integer, data_struct )
  if (indx_x, indx_y) in positions
    xs = [ indx_x+1, indx_x, indx_x-1, indx_x ]
    ys = [ indx_y, indx_y+1, indx_y, indx_y-1 ]
    expand!.( Ref(positions), Ref(concentrations), Ref(dtm), xs, ys, data_struct )
    return nothing
  else

 # SOSTITUIRE CON UNA O UN'INSIEME DI FUNZIONI CHE MODIFICHINO I VALORI NECESSARI NELLO STRUCT
    Δx, Δy = Functions.toCoords(dtm, positions[1][1], positions[1][2]) - Functions.toCoords(dtm, indx_x, indx_y)
    dir = deg2rad(plume.wind_direction)
    sindir = sin(dir)
    cosdir = cos(dir)
    data_struct.x = Δy * sindir + Δy * cosdir
    data_struct.y = Δx * cosdir - Δy * sindir

    result = calcSediment!(sediment)


 # SOSTITUIRE CON UN CONTROLLO CHE DIPENDA DA PARAMETRI ESTERNI
    if round(result, digits=5) > 0
        push!( positions, (ind_x, ind_y) )
        push!( results, result )
        xs = [ indx_x+1, indx_x, indx_x-1, indx_x ]
        ys = [ indx_y, indx_y+1, indx_y, indx_y-1 ]
        expand!.( Ref(positions), Ref(results), Ref(dtm), xs, ys, data_struct )
    end
    return nothing
  end
end
=#






#------------------------------------------------ TESTING------------------------------------------------------------------------------
#=
import ArchGDAL as agd
import GeoArrays as ga


dtm_file = split( @__DIR__ , "\\Porting\\")[1] * "\\Mappe\\DTM_32.tiff"
dtm = agd.read(dtm_file)
gdtm = ga.read(dtm_file)

x, y = ga.coords(gdtm, [4000, 6000])
x1, y1 = toCoords(dtm, 4000, 6000)

ga.indices(gdtm, [x, y])
toIndexes(dtm, x1, y1)
=#























mutable struct DAF
    secondary_source_concentration::Float64
    x::Float64
    y::Float64
    α_x::Float64
    α_y::Float64
    decay_coeff::Float64
    darcy_velocity::Float64
    kd
    soil_density::Float64
    tera_e
    orthogonal_extension::Float64
    time
 # AGGIUNTE ALL'ORIGINALE PER SEMPLIFICARE I CALCOLI
    acquifer_flow_direction::Float64
    algorithm::Symbol
    option::Symbol
  
    R
    DAF
    DAF_tot
  
    DAF(secondary_source_concentration,x,y,α_x,α_y,decay_coeff,darcy_velocity,kd,soil_density,tera_e,orthogonal_extension,time,acquifer_flow_direction,algorithm,option) = new(secondary_source_concentration,x,y,α_x,α_y,decay_coeff,darcy_velocity,kd,soil_density,tera_e,orthogonal_extension,time,acquifer_flow_direction,algorithm,option)
end

mutable struct Plume
    """docstring for element"""
    # concentration: concentrazione inquinante m3/sec
    # d: distanza (coordinata x)
    # y,z: coordinate in metri
    # stability: classe di stabilità
    # wind_speed: velocità del vento nella direzione x
    # stack_height: altezza camino
    # gas_speed: velocità gas
    # stack_diameter: diametro camino
    # smoke_temperature: temperatura gas all'uscita del camino
    # temperature: temperatura ambiente
  
    concentration::Real
    d::Real
    y::Real
    z::Real
    stability::AbstractString
    outdoor::AbstractString
    stack_height::Real
    stack_diameter::Real
    wind_direction::Integer
    wind_speed::Real
    gas_speed::Real
    smoke_temperature::Real
    temperature::Real
    max_domain::Real
  
    H
    σy
    σz
    g1
    g2
    
    Plume(concentration,d,y,z,stability,outdoor,stack_height,stack_diameter,wind_direction,wind_speed,gas_speed,smoke_temperature,temperature,max_domain)=new(concentration,d,y,z,stability,outdoor,stack_height,stack_diameter,wind_direction,wind_speed,gas_speed,smoke_temperature,temperature,max_domain)
end

mutable struct Lake
    concentration::Float64
    time::Time
    distance_x::Float64
    distance_y::Float64
    fickian_x::Float64
    fickian_y::Float64
    velocity_x::Float64
    velocity_y::Float64
    λk::Float64

    C::Float64

    Lake( concentration, time, distance_x, distance_y, fickian_x, fickian_y, velocity_x, velocity_y, λk ) = new( concentration, time, distance_x, distance_y, fickian_x, fickian_y, velocity_x, velocity_y, λk )
end

@with_kw mutable struct Sediment
    """docstring for element"""
  
    # dredged_mass: input sedimento (kg/sec)
    # t: tempo finale
    # h: profondità metri
    # Dx,Dy: coefficienti di diffusione
    # x0,y0: coordinate sorgente
    # x,y: coordinate target point
    # V: velocità media corrente
    # w: velocità di sedimentazione
    # dt: delta t per la discretizzazione dell'integrale
    # U: amplitude of the oscillatory current
    # tide: tidal cycle es. 12 (ore)
  
    dredged_mass::Float64
    time
    mean_depth::Int
    x_dispersion_coeff::Float64
    y_dispersion_coeff::Float64
    x::Int
    y::Int
    mean_flow_speed::Float64
    mean_sedimentation_velocity::Float64
    time_intreval::Int
    current_oscillatory_amplitude::Float64 = 0.0
    tide::Int = 0
  
    ω = 0
    ew


    function Sediment(dredged_mass, time, mean_depth, x_dispersion_coeff, y_dispersion_coeff, x, y, mean_flow_speed, mean_sedimentation_velocity, time_intreval, current_oscillatory_amplitude, tide)
        if current_oscillatory_amplitude > 0 && tide > 0
            ω = 2π/tide
            return new(dredged_mass, time, mean_depth, x_dispersion_coeff, y_dispersion_coeff, x, y, mean_flow_speed, mean_sedimentation_velocity, time_intreval, current_oscillatory_amplitude, tide, ω)
        end
    end
end












function expand!( positions::AbstractVector, results::AbstractVector, dtm::AbstractArray, indx_x::Integer, indx_y::Integer, object )
    if (indx_x, indx_y) in positions
        xs = [ indx_x, indx_x-1, indx_x ]
        ys = [ indx_y+1, indx_y, indx_y-1 ]
  
        expand!( positions, concentrations, dtm, indx_x+1, indx_y, object )
        expand!.( Ref(positions), Ref(concentrations), Ref(dtm), xs, ys, deepcopy(object) )
        return nothing
    else
        result = compue_result!(
            Δx * cosdir - Δx * sindir,
            Δy * sindir + Δy * cosdir,
            object
        )
        if condition(result)
            push!( positions, (ind_x, ind_y) )
            push!( results, result )
            xs = [ indx_x, indx_x-1, indx_x ]
            ys = [ indx_y+1, indx_y, indx_y-1 ]
            expand!( positions, concentrations, dtm, indx_x+1, indx_y, object )
            expand!.( Ref(positions), Ref(results), Ref(dtm), xs, ys, deppcopy(object) )
        end
        return nothing
    end
end


# FORSE SI PUO' PASSARE DIRETTAMENTE LE COORDINATE E USARE QUELLE
function compute_result!( x0, y0, xi, yi, sediment::Sediment )
    Δx, Δy = Functions.toCoords(dtm, x0, y0) - Functions.toCoords(dtm, xi, yi)
    dir = deg2rad(object.dir)
    sindir = sin(dir)
    cosdir = cos(dir)
    sediment.x = x
    sediment.y = y
    return calcSediment!(sediment)
end

function compute_result!( x, y, daf::DAF )
    daf.x = x
    daf.y = y
    return calcDAF!(daf)
end

function compute_result!( x, y, lake::Lake )
    lake.distance_x = y
    lake.distance_y = x
    return calc_concentration!(lake)
end





#   PLUMES
function expand!( plume::Plume )
    if (indx_x, indx_y) in positions
    else

      plume.d = Δy * sindir + Δy * cosdir
      plume.y = Δx * cosdir - Δy * sindir
      plume.z = agd.getband(dtm, 1)[indx_x, indx_y]
      concentration = calcPlume!(plume)

    end
end


# SERVONO IL DTM E GLI INDICI DELLA CELLA CORRENTE
function compute_result!( x, y, plume::Plume )
    plume.d = y
    plume.y = x
    plume.z = agd.getband(dtm, 1)[]
    return calc_concentration!(lake)
end

















end # module