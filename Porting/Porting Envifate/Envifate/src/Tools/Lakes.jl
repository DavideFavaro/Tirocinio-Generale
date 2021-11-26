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

#= CON LA MODIFICA ALLE FUNZIONI QUESTO NON SERVE
mutable struct Lake
    """docstring for element"""
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
=#

# QUESTI SI POSSONO SPOSTARE SU UN FILE A PARTE 
srid = [ 3003, 3004, 32632, 32633, 3857, 4326 ]
classiwind = Dict(
    "N" => 0,
    "NE" => 45,
    "E" => 90,
    "SE" => 135,
    "S" => 180,
    "SW" => 225,
    "W" => 270,
    "NW" => 315
)
classiwind2 = Dict(
    0 => 180,
    45 => 225,
    90 => 270,
    135 => 315,
    180 => 0,
    225 => 45,
    270 => 90,
    315 => 135
)


#= VERSIONE CON LO STRUCT
function calc_concentration!( l::Lake )
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
    return l.C
end
=#
function calc_concentration( concentration::Real, distance_x::Real, distance_y::Real, fickian_x::Real, fickian_y::Real, velocity_x::Real, velocity_y::Real, λk::Real, time )
    c1, c2 = ( (distance_x, distance_y) - ( (velocity_x, velocity_y) * time ) )^2 / ( 4 * (fickian_x, fickian_y) * time )
    c3 = ℯ^( -(c1 + c2) - (λk * time) )
    c4 = concentration / ( 4π * time * √(fickian_x * fickian_y) )
    return c4 * c3
end



function run_lake( source, area, xw, pollutant_mass, flow_mean_speed, resolution::In64, hours::Int64,
                   fickian_x::Real=0.05, fickian_y::Real=0.05, λk::Real=0.0, output_path::AbstractString=".\\otput_model.tiff" )

    if agd.getspatialref(area) != agd.getspatialref(source)
        throw(DomainError("The reference systems are not uniform. Aborting analysis."))
    end


    hours *= 3600
    refsys = agd.importEPSG(agd.fromWKT(agd.getspatialref(source)))
    velocity_x = √( round( flow_mean_speed * cos(deg2rad(xw)), digits=3 )^2 )
    velocity_y = √( round( flow_mean_speed * sin(deg2rad(xw)), digits=3 )^2 )


 """
    path_layer = dialog.areastudio.dataProvider().dataSourceUri()
    path = split( path_layer, "|" )
    source_ds = ogr.Open(path[0])
    area_layer = source_ds.GetLayer()

    x_min, y_min, x_max, y_max = round.( Int64, [ area_layer.GetExtent()[1], area_layer.GetExtent()[3], area_layer.GetExtent()[2], area_layer.GetExtent()[4] ] )
 """
    area_layer = agd.getlayer(area, 0)
 # NON FUNZIONANTE / DA ELIMINARE
    x_min, y_min, x_max, y_max = agd.envelope(area_layer)
    valNoData = -9999

    # Create the destination data source
    x_res = ( x_max - x_min ) / resolution
    y_res = ( y_max - y_min ) / resolution

    gtiff_driver = agd.getdriver("GTiff")
    target_ds = agd.create( path_output, gtiff_driver, round(Int64, x_res), round(Int64, y_res), 1, agd.GDAL.GDT_Float32 )
    agd.setgeotransform!( target_ds, [ x_min, resolution, 0.0, y_max, 0.0, -resolution ] )
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
    xsize = agd.width(band1)
    ysize = agd.height(band1)

    polygons = collect(agd.getfeature(area))
    feature = collect(agd.getfeature(source))
    geom = agd.getgeom(feature[1])
    x_source = agd.getx(geom, 0)
    y_source = agd.gety(geom, 0)

    rows = ysize - 1
    cols = xsize - 1

    for row in 1:rows
        for col in 1:cols
            x, y = (col, row) * resolution + (x_min, y_min) .+ (resolution/2)
            control_point = agd.createpoint(x,y)
            for polygon in polygons
                if agd.within( control_point, agd.getgeom(polygon) )
                    Δx = x - x_source
                    Δy = y - y_source
                    true_x = Δx * cos(deg2rad(xw)) - Δy * sin(deg2rad(xw))
                    true_y = Δx * sin(deg2rad(xw)) + Δy * cos(deg2rad(xw))

                 #= VERSIONE CON STRUCT `Lake`
                    element = Lakes.Lake( pollutant_mass, hours, true_x, true_y, fickian_x, fickian_y, velocity_x, velocity_y, λk )
                    cfinal = calc_concentration!(element)
                 =#
                    cfinal = calc_concentration( pollutant_mass, true_x,  true_y, fickian_x, fickian_y, velocity_x, velocity_y, λk, hours)

                    outData[row,col] = cfinal
                else
                    outData[row,col] = 0
                end
            end
        end
    end



    
    outData_raster = outData[::-1]
    band.WriteArray(outData_raster)
    astats = band.GetStatistics(0, 1)

    band = nothing
    target_ds = nothing

    base_raster_name = basename(path_output)
    raster_name = splitext(base_raster_name)[0]
    dialog.iface.addRasterLayer( dialog.path_output, raster_name )

    layer = nothing
    for lyr in collect( QgsProject.instance().mapLayers().values() )
        if lyr.name() == raster_name
            layer = lyr
        end
    end

    functions.applystyle( layer, "gr", 0.5 )
end



end # module