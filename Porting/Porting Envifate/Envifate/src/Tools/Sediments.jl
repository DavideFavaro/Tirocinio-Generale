module Sediments
"""
Module for marine sedimentation analysis
"""


import ArchGDAL as agd
using ArgParse
using Parameters
using Dates


include("../Library/Functions.jl")



#   NON SO COME SI DOCUMENTANO GLI STRUCT
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
    time::Int64
    mean_depth::Float64
    x_dispersion_coeff::Float64
    y_dispersion_coeff::Float64
    x::Float64
    y::Float64
    mean_flow_speed::Float64
    flow_direction::Float64
    mean_sedimentation_velocity::Float64
    time_intreval::Int64
    current_oscillatory_amplitude::Float64 = 0.0
    tide::Int64 = 0
  
    ω = 0
    ew


    function Sediment(dredged_mass, time, mean_depth, x_dispersion_coeff, y_dispersion_coeff, x, y, mean_flow_speed, flow_direction, mean_sedimentation_velocity, time_intreval, current_oscillatory_amplitude, tide)
        if current_oscillatory_amplitude > 0 && tide > 0
            ω = 2π/tide
            return new(dredged_mass, time, mean_depth, x_dispersion_coeff, y_dispersion_coeff, x, y, mean_flow_speed, flow_direction, mean_sedimentation_velocity, time_intreval, current_oscillatory_amplitude, tide, ω)
        end
    end
end



function calc_q( s::Sediment )
    return s.dredged_mass / ( 4π *s.mean_depth * isqrt(s.x_dispersion_coeff * s.y_dispersion_coeff) )
end


function calc_e!( s::Sediment, i )
    s.ew = s.ω > 0 ? s.current_oscillatory_amplitude / ( s.ω * cos(deg2rad(s.ω)) - cos(deg2rad(s.ω * i *s.time_intreval)) ) : 0
    e1 = ℯ^(-(( s.x - s.mean_flow_speed * ( s.time - i * s.time_intreval) + s.ew ) / ( 4s.x_dispersion_coeff * (s.time - i * s.time_intreval) ) ))
    e2 = ℯ^(-( s.y^2 / ( 4s.y_dispersion_coeff * (s.time - i * s.time_intreval) ) ) - ( (s.mean_sedimentation_velocity * (s.time - i * s.time_intreval)) / s.mean_depth ) )
    return e1*e2
end


function calcSediment!( s::Sediment )
    if s.x <= 0
        return 0.0
    else
        q = calc_q(s)
        n = round( Int64, s.time / s.time_intreval )
        csum = 0
        for i in 1:n
            #   csum += calc_e!(s, i) * ( 1 / ( s.time - ( i * s.time_intreval ) ) )
            csum += calc_e!(s, i) / ( s.time - ( i * s.time_intreval ) )
        end
        return q * csum * s.time_intreval
    end
end



"""
    compute_result!( dtm::AbstractArray, r0::Integer, c0::Integer, ri::Integer, ci::Integer, sediment::Sediment )

Given the raster `dtm` and the indexes (`r0`, `c0`) of the source, modify the postion values of object `sediment` and return the concentration at indexes (`ri`, `ci`)
"""
function compute_result!( dtm::AbstractArray, r0::Integer, c0::Integer, ri::Integer, ci::Integer, sediment::Sediment )
    sediment.x, sediment.y = Functions.compute_position(dtm, r0, c0, ri, ci, sediment.flow_direction)
    return calcSediment!(sediment)
end



"""
    run_sediment( source, resolution::Integer, mean_flow_speed::Real, mean_depth::Real, x_dispersion_coeff::Real, y_dispersion_coeff::Real,
                  dredged_mass::Real, flow_direction::Real, mean_sedimentation_velocity::Real, time::Integer, time_intreval::Integer,
                  current_oscillatory_amplitude::Integer=0, tide::Integer=0, output_path::AbstractString=".\\output_model.tiff" )

Run a simulation of plumes of turbidity induced by dredging.

# Arguments
- `source`: dredging source point.
- `resolution::Integer`: size of a cell in meters.
- `mean_flow_speed::Real`: speed of the flowing water.
- `mean_depth::Real`: depth in meters.
- `x_dispersion_coeff::Real`: coefficient of dispersion along the x axis.
- `y_dispersion_coeff::Real,`: coefficient of dispersion along y axis.
- `dredged_mass::Real`: initial mass of the dredged substance.
- `flow_direction::Real`: direction of the as an angle, in degrees.
- `mean_sedimentation_velocity::Real`: velocity of sedimentation.
- `time::Integer`: start time for the model.
- `time_intreval::Integer`: length of an epoch.
- `current_oscillatory_amplitude::Integer=0`: water oscillatory amplitude.
- `tide::Integer=0`: value of tide.
- `output_path::AbstractString=".\\output_model.tiff"`: path of the resulting raster.
"""
                    #                                     v                      h                 dx                        dy                        q                   dir
function run_sediment( source, resolution::Integer, mean_flow_speed::Real, mean_depth::Real, x_dispersion_coeff::Real, y_dispersion_coeff::Real, dredged_mass::Real, flow_direction::Real,
                    #  w                                  t / time       dt                      u
                       mean_sedimentation_velocity::Real, time::Integer, time_intreval::Integer, current_oscillatory_amplitude::Integer=0, tide::Integer=0, output_path::AbstractString=".\\output_model.tiff" )

 # messaggio+='ALGORITMO UTILIZZATO: Shao (Shao, Dongdong, et al. "Modeling dredging-induced turbidity plumes in the far field under oscillatory tidal currents." Journal of Waterway, Port, Coastal, and Ocean Engineering 143.3 (2016))\n\n'

    feature = collect(agd.getlayer(source, 0))[1]
    geom = agd.getgeom(feature)

    if agd.geomdim(geom) != 0
        throw(DomainError(source, "`source` must be a point"))
    end
    
    x_source = agd.getx(geom, 0)
    y_source = agd.gety(geom, 0)
    r_source, c_source = toCoords(dtm, x_source, y_source)
    
    refsys = agd.getspatialref(geom)

    #   start_time = time.time()


    points = [ (r_source, c_source) ]
 # NON SO SE SIA IL VALORE CORRETTO DA INSERIRE
    values = [ dredged_mass ]
    element = Sediment(dredged_mass, time, mean_depth, x_dispersion_coeff, y_dispersion_coeff, 0.0, 0.0, mean_flow_speed, flow_direction, mean_sedimentation_velocity, time_intreval, current_oscillatory_amplitude, tide)
    expand!( points, values, dem, r_source, c_source, element )


    points = [ (r_source, c_source) ]
    values = [ concentration ]
 # NON CI INTERESSA DARE DEI VALORI COERENTI A d, y, E z PERCHE' AD OGNI CHIAMATA expand! LI RESETTA
    plume = Plume(concentration, 0.0, 0.0, 0.0, stability, outdoor, wind_speed, stack_height, stack_diameter, gas_speed, smoke_temperature, temperature, wind_direction,0.0) 
    expand!(points, values, dtm, r_source, c_source, plume)

    maxR = maximum( point -> point[1], points )
    minR = minimum( point -> point[1], points )
    maxC = maximum( point -> point[2], points )
    minC = minimum( point -> point[2], points )

    rows = maxR - minR
    cols = maxC - minC
    minX, maxY = toCoords(dtm, minX, maxY)

    gtiff_driver = agd.getdriver("GTiff")
    target_ds = agd.create( output_path, gtiff_driver, rows, cols, 1, agd.GDAL.GDT_Float32 )
    agd.setgeotransform!(target_ds, [ minX, resolution, 0.0, maxY, 0.0, -resolution ])
    agd.setproj!(target_ds, refsys)
 """ NON SO QUALE SIA IL COMANDO PER SETTARE I METADATI CON `ArchGDAL`
    target_ds.SetMetadata(
        Dict(
            "credits" => "Envifate - Francesco Geri, Oscar Cainelli, Paolo Zatelli, Gianluca Salogni, Marco Ciolli - DICAM Università degli Studi di Trento - Regione Veneto",
            "modulo" => "Analisi sedimentazione marina",
            "descrizione" => "Simulazione di sedimento disperso in ambiente marino",
            "srs" => refsys,
            "data" => today()
        )
    )
 """
    valNoData = -9999.0
    band1 = agd.getband(target_ds, 1)
    agd.setnodatavalue!( band1, Float64(valNoData) )
    agd.fillraster!(band, valNoData)
    band = agd.read(band1)

    for (point, value) in zip(points[i], values[i])
        r, c = point - (minR, minC)
        band[r, c] = value
    end
end



end # module
