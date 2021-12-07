module Lakes

# -*- coding: utf-8 -*-
#=
/***************************************************************************
 OpenRisk
                                 A QGIS plugin
 Open Risk: Open source tool for environmental risk analysis
                              -------------------
        begin                : 2016-07-15
        git sha              : $Format:%H$
        copyright            : (C) 2016 by Francesco Geri
        email                : fgeri@icloud.com
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
=#


import ArchGDAL as agd
using ArgParse
using Dates
using Sys

include("../Library/Functions.jl")

export run_lake



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
Recursively compute the concentration of each point and add the value and its indexes to positions
"""
function expand!( positions::AbstractVector, results::AbstractVector, dtm, indx_x::Integer, indx_y::Integer, lake::Lake )
    if (indx_x, indx_y) in positions
        xs = [ indx_x+1, indx_x, indx_x-1, indx_x ]
        ys = [ indx_y, indx_y+1, indx_y, indx_y-1 ]
        expand!.( Ref(positions), Ref(concentrations), Ref(dtm), xs, ys, lake )
        return nothing
    else
        Δx, Δy = Functions.toCoords(dtm, positions[1][1], positions[1][2]) - Functions.toCoords(dtm, indx_x, indx_y)
        dir = deg2rad(plume.wind_direction)
        sindir = sin(dir)
        cosdir = cos(dir)
        lake.distance_x = Δy * sindir + Δy * cosdir
        lake.distance_y = Δx * cosdir - Δy * sindir
        result = calc_concentration!(lake)

        if round(result, digits=5) > 0
            push!( positions, (ind_x, ind_y) )
            push!( results, result )
            xs = [ indx_x+1, indx_x, indx_x-1, indx_x ]
            ys = [ indx_y, indx_y+1, indx_y, indx_y-1 ]
            expand!.( Ref(positions), Ref(results), Ref(dtm), xs, ys, lake )
        end
        return nothing
    end
end



function run_lake( source, xw, pollutant_mass, flow_mean_speed, resolution::In64, hours::Int64,
                   fickian_x::Real=0.05, fickian_y::Real=0.05, λk::Real=0.0, output_path::AbstractString=".\\lake_otput_model.tiff" )

    hours *= 3600
    refsys = agd.getspatialref(source)
    velocity_x = √( round( flow_mean_speed * cos(deg2rad(xw)), digits=3 )^2 )
    velocity_y = √( round( flow_mean_speed * sin(deg2rad(xw)), digits=3 )^2 )

    feature = collect(agd.getfeature(source))
    geom = agd.getgeom(feature[1])
    x_source = agd.getx(geom, 0)
    y_source = agd.gety(geom, 0)
    r_source, c_source = toIndexes(dem, x_source, y_source)

    points = [ (r_source, c_source) ]
    values = [ pollutant_mass ]
    lake = Lake( pollutant_mass, hours, 0, 0, fickian_x, fickian_y, velocity_x, velocity_y, λk )
    expand!(points, values, dem,  r_source, c_source, lake)

    maxR = maximum( point -> point[1], points )
    minR = minimum( point -> point[1], points )
    maxC = maximum( point -> point[2], points )
    minC = minimum( point -> point[2], points )
  
    rows = maxR - minR
    cols = maxC - minC
    minX, maxY = toCoords(dtm, minR, maxC)

    valNoData = -9999

    gtiff_driver = agd.getdriver("GTiff")
    target_ds = agd.create( output_path, gtiff_driver, rows, cols, 1, agd.GDAL.GDT_Float32 )
    agd.setgeotransform!( target_ds, [ minX, resolution, 0.0, maxY, 0.0, -resolution ] )
    agd.setproj!(target_ds, refsys)
 """ NON SO QUALE SIA IL COMANDO PER SETTARE I METADATI CON `ArchGDAL`
    target_ds.SetMetadata(
        Dict( 
            "credits" => "Envifate - Francesco Geri, Oscar Cainelli, Paolo Zatelli, Gianluca Salogni, Marco Ciolli - DICAM Università degli Studi di Trento - Regione Veneto",
            "modulo" => "Dispersione in laghi e bacini",
            "descrizione" => "Simulazione di dispersione inquinante all\'interno di un corpo idrico superficiale fermo",
            "srs" => refsys,
            "data" => today()
        )
    )
 """
    band1 = agd.getband(target_ds, 1)
    agd.setnodatavalue!( band1, convert(Float64, valNoData) )
    band = agd.read(band1)
    agd.fillraster!(band, valNoData)

    for (point, value) in zip(points, values)
        r, c = point - (minR, minC)
        band[r, c] = value
    end

    agd.write!( target_ds, band )
end



end # module