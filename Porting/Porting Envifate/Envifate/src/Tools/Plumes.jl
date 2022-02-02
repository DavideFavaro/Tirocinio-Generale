module Plumes
"""
Module for the modelling of pollutants' dispersion in the atmosphere
"""



import ArchGDAL as agd
using ArgParse
using Dates


include("..\\Library\\Functions.jl")



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

    concentration::Float64
    d::Float64
    y::Float64
    z::Float64
    stability::String
    outdoor::String
    stack_height::Float64
    stack_diameter::Float64
    wind_direction::Integer
    wind_speed::Float64
    gas_speed::Float64
    smoke_temperature::Float64
    temperature::Float64
    max_domain::Float64

    H
    σy
    σz
    g1
    g2
    
    Plume(concentration,d,y,z,stability,outdoor,stack_height,stack_diameter,wind_direction,wind_speed,gas_speed,smoke_temperature,temperature,max_domain)=new(concentration,d,y,z,stability,outdoor,stack_height,stack_diameter,wind_direction,wind_speed,gas_speed,smoke_temperature,temperature,max_domain)
end



function calc_h!( p::Plume )
    try
        fb = 9.81 * ( (p.stack_diameter * p.gas_speed) / 4 ) * ( ( p.smoke_temperature / p.temperature ) / p.smoke_temperature )
        Δh = 1.6 * fb^0.333333 * p.d^0.666667
        p.H = p.stack_height + Δh
    catch
        p.H = p.stack_height
    end
    return p.H
end
#= VERSIONE SENZA LO STRUCT
function calc_h( d::Real, stack_diameter::Real, gas_speed::Real, smoke_temperature::Real, temperature::Real )
    try
        fb = 9.81 * ( (stack_diameter * gas_speed) / 4 ) * ( ( smoke_temperature / temperature ) / smoke_temperature )
        Δh = 1.6fb^(1/3) * d^(2/3)
        return stack_height + Δh
    catch
        return stack_height
    end
end
=#


function calc_σ!(p::Plume)
    σ_values = Functions.air_extract( p.c_stability, p.outdoor )
    σy1 = σ_values[0]
    σy2 = σ_values[1]
    σyexp = σ_values[2]
    σz1 = σ_values[3]
    σz2 = σ_values[4]
    σzexp = σ_values[5]

    p.σy = ( σy1 * p.d ) / ( 1 + σy2 * p.d )^σyexp
    p.σz = ( σz1 * p.d ) / ( 1 + σz2 * p.d )^σzexp
    return p.σy, p.σz
end
#= VERSIONE SENZA LO STRUCT
function calc_σ( d::Real, stability_class::AbstractString, outdoor::AbstractString )
    σ_values = Functions.air_extract( stability_class, outdoor )
    σy = ( σ_values[0] * d ) / ( 1 + σ_values[1] * d )^σ_values[2]
    σz = ( σ_values[3] * d ) / ( 1 + σ_values[4] * d )^σ_values[5]
    return σy, σz
end
=#


function calc_g!(p::Plume)
    p.g1 = ℯ^( ( -0.5 * p.y^2 ) / p.σy^2 )
    p.g2 = ℯ^( ( -0.5 * (p.z - p.stack_height)^2 ) / p.σz^2 ) + ℯ^( ( -0.5 * (p.z + p.stack_height)^2 ) / p.σz^2 )
    return p.g1, p.g2
end
#= VERSIONE SENZA LO STRUCT
function calc_g( y::Real, z::Real, σy::Real, σz::Real, stack_height::Real )
    g1 = ℯ^( ( -0.5y^2 ) / σy^2 )
    g2 = ℯ^( ( -0.5(z - stack_height)^2 ) / σz^2 ) + ℯ^( ( -0.5(z + stack_height)^2 ) / σz^2 )
    return g1, g2
end
=#


function calc_C!( p::Plume )
    p.C = ( 100p.concentration / 3600p.wind_speed ) * ( (p.g1 * p.g2) / ( 2π * p.σy * p.σz ) )
    return p.C
end
#= VERSIONE SENZA LO STRUCT
function calc_C( concentration::Real, wind_speed::Real, σy::Real, σz::Real, g1::Real, g2::Real )
    return ( 100concentration / 3600wind_speed ) * ( (g1 * g2) / ( 2π * σy * σz ) )
end
=#


function calcPlume!(p::Plume)
    if p.d <= 0
        return 0.0
    else
        calc_σ!(element)
        calc_g!(element)
        calc_h!(element)
        cfinal = calc_C!(element)
        if cfinal > max_dominio
            p.max_domain = cfinal
        end
        return cfinal
    end
end



"""
Given the raster `dtm` and the indexes (`r0`, `c0`) of the source, modify the postion values of object `sediment` and return the concentration at indexes (`ri`, `ci`)
"""
function compute_result!( dtm::AbstractArray, r0::Integer, c0::Integer, ri::Integer, ci::Integer, plume::Plume )
    plume.d, plume.y = Functions.compute_position(dtm, r0, c0, ri, ci, plume.wind_direction)
    plume.z = agd.getband(dtm, 1)[ri, ci]
    return calc_concentration!(lake)
end



"""
    run_plume( dem, source, stability::AbstractString, outdoor::AbstractString, resolution::Integer, wind_direction, concentration::Real, wind_speed::Real, stack_height::Real, 
               gas_speed::Real=0.0, stack_diameter::Real=0.0, smoke_temperature::Real=0.0,  temperature::Real=0.0, output_path::AbstractString=".\\otput_model_plume.tiff" )

# Arguments
- `dem`: raster containing the height of the terrain in each cell.
- `source`: source point of the plume.
- `stability::AbstractString`: information on the weather.
- `outdoor::AbstractString`: 
- `resolution::Integer`: size of the cell in meters.
- `wind_direction``: angle of direction of the wind in degrees.
- `concentration::Real`: concentration of contaminants at the source.
- `wind_speed::Real`: average wind speed.
- `stack_height::Real `: height of the stack, or height of the source of the plume.
- `gas_speed::Real=0.0`: movement speed of the gas.
- `stack_diameter::Real=0.0`: diameter of the stack emiting the plume.
- `smoke_temperature::Real=0.0`: temperature of the fumes.
- `temperature::Real=0.0`: average temperature of the environment.
- `output_path::AbstractString=".\\otput_model_plume.tiff"`: output file path. 
"""
         #                                                                                                x_w             q / text_conc        u / wspeed        h_s / height
function run_plume( dem, source, stability::AbstractString, outdoor::AbstractString, resolution::Integer, wind_direction, concentration::Real, wind_speed::Real, stack_height::Real, 
                  # v_s / gspeed         d_s / diameter            t_s / temp                    t_a / etemp
                    gas_speed::Real=0.0, stack_diameter::Real=0.0, smoke_temperature::Real=0.0,  temperature::Real=0.0, output_path::AbstractString=".\\otput_model_plume.tiff" )

    geom = agd.getgeom( collect(agd.getlayer(source, 0))[1] )

    if agd.geomdim(geom) != 0
        throw(DomainError(source, "`source` must be a point"))
    end

    refsys = agd.getproj(dem)

    if agd.importWKT(refsys) != agd.getspatialref(geom)
        throw(DomainError("The reference systems are not uniform. Aborting analysis."))
    end
  
 #= PENSO SIA ROBA UTILIZZATA IN DEI PRINT CHE SONO STATI RIMOSSI
    lst_fields = [
      "sf_ing",
      "sf_inal",
      "iur",
      "rfd_ing",
      "rfd_inal",
      "rfc"
    ]
    lst_toxic_msg = [ 
        "Slope Factor per ingestione", 
        "Slope Factor per inalazione",
        "Inhalation Unit Risk",
        "Reference Dose per ingestione",
        "Reference Dose per inalazione",
        "Reference Concentration"
    ]
    toxic = Functions.substance_extract(contaminant, lst_fields, "..\\Library\\")
 =#

    x_source = agd.getx(geom, 0)
    y_source = agd.gety(geom, 0)
    r_source, c_source = toIndexes(geomtransform, x_source, y_source)

    # start_time = time.time()

    points = [ (r_source, c_source) ]
    values = [ concentration ]
 # NON CI INTERESSA DARE DEI VALORI COERENTI A d, y, E z PERCHE' AD OGNI CHIAMATA expand! LI RESETTA
    plume = Plume(concentration, 0.0, 0.0, 0.0, stability, outdoor, wind_speed, stack_height, stack_diameter, gas_speed, smoke_temperature, temperature, wind_direction,0.0) 
    expand!(points, values, dtm, r_source, c_source, plume)

    maxR = maximum( point -> point[1], points )
    minR = minimum( point -> point[1], points )
    maxC = maximum( point -> point[2], points )
    minC = minimum( point -> point[2], points )

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