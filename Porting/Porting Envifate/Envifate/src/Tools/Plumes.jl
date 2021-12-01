module Plumes

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



function calc_h!( p::Plume )
  try
    fb = 9.81 * ( (p.stack_diameter * p.gas_speed) / 4 ) * ( ( p.smoke_temperature / p.temperature ) / p.smoke_temperature )

    delta_h = 1.6 * fb^0.333333 * p.d^0.666667

    p.H = p.stack_height + delta_h
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

  p.σy = ( σy1 * p.d ) / (1 + σy2 * p.d )^σyexp
  p.σz = ( σz1 * p.d ) / (1 + σz2 * p.d)^σzexp
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
  p.g1 = exp( ( -0.5 * p.y^2 ) / p.σy^2 )
  p.g2 = exp( ( -0.5 * (p.z - p.stack_height)^2 ) / p.σz^2 ) + exp( ( -0.5 * (p.z + p.stack_height)^2 ) / p.σz^2 )
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
Recursively compute the concentration of each point and add the value and its indexes to positions
"""
function expand!( positions::AbstractVector, results::AbstractVector, dtm, indx_x::Integer, indx_y::Integer, plume::Plume )
  if (indx_x, indx_y) in positions
    xs = [ indx_x+1, indx_x, indx_x-1, indx_x ]
    ys = [ indx_y, indx_y+1, indx_y, indx_y-1 ]
    expand!.( Ref(positions), Ref(concentrations), Ref(dtm), xs, ys, plume )
    return nothing
  else
    Δx, Δy = Functions.toCoords(dtm, positions[1][1], positions[1][2]) - Functions.toCoords(dtm, indx_x, indx_y)
    dir = deg2rad(plume.wind_direction)
    sindir = sin(dir)
    cosdir = cos(dir)
 # SETTANO d A true_y E y A true_x, NON SO SE SIA GIUSTO
    plume.d = Δy * sindir + Δy * cosdir
    plume.y = Δx * cosdir - Δy * sindir
    plume.z = agd.getband(dtm, 1)[indx_x, indx_y]
    concentration = calcPlume!(plume)
    if round(concentration, digits=5) > 0
        push!( positions, (ind_x, ind_y) )
        push!( results, concentration )
        xs = [ indx_x+1, indx_x, indx_x-1, indx_x ]
        ys = [ indx_y, indx_y+1, indx_y, indx_y-1 ]
        expand!.( Ref(positions), Ref(results), Ref(dtm), xs, ys, plume )
    end
    return nothing
  end
end







#=
    def popolacombo(self):
        self.progressBar.setValue(0)
        self.combo_source.clear()
        self.combo_outdoor.clear()
        self.combo_stability.clear()
        self.combo_maindirwind.clear()
        self.line_output.clear()
        self.combo_contaminant.clear()

        self.allLayers = self.canvas.layers()
        self.listalayers=dict()
        elementovuoto="No file"

        for i in self.allLayers:
            self.listalayers[i.name()]=i
            if i.type() == QgsMapLayer.VectorLayer:
                self.combo_source.addItem(str(i.name()))
                self.combo_bound.addItem(str(i.name()))
            if i.type()==QgsMapLayer.RasterLayer:
                self.combo_dem.addItem(str(i.name()))

        # pyqtRemoveInputHook()
        # pdb.set_trace()
        # conn = sqlite3.connect(os.path.dirname(__file__)+"/../library/substance.db")
        # cursor=conn.cursor()
        # query_substance="select id,nome from substance"
        # conn.close()
        # cursor.execute(query_substance)
        # sql_fetch=cursor.fetchall()

        # self.inquinanti=dict()
        # for row in sql_fetch:
        #     self.inquinanti[row[1]]=row[0]
        #     self.combo_contaminant.addItem(row[1])
        self.combo_stability.addItem("class a")
        self.combo_stability.addItem("class b")
        self.combo_stability.addItem("class c")
        self.combo_stability.addItem("class d")
        self.combo_stability.addItem("class e")
        self.combo_stability.addItem("class f")

        self.combo_outdoor.addItem("country")
        self.combo_outdoor.addItem("urban")


        self.combo_maindirwind.addItem("N")
        self.combo_maindirwind.addItem("NE")
        self.combo_maindirwind.addItem("E")
        self.combo_maindirwind.addItem("SE")
        self.combo_maindirwind.addItem("S")
        self.combo_maindirwind.addItem("SW")
        self.combo_maindirwind.addItem("W")
        self.combo_maindirwind.addItem("NW")


        conn = sqlite3.connect(os.path.dirname(__file__)+"/../library/substance.db")
        cursor=conn.cursor()
        query_substance="select id,nome from substance"
        cursor.execute(query_substance)
        self.sql_fetch_inq=cursor.fetchall()

        self.inquinanti=dict()
        for row in self.sql_fetch_inq:
            self.inquinanti[row[1]]=row[0]
            self.combo_contaminant.addItem(row[1])
        conn.close()
=#



         #                                                                                                x_w             q / text_conc        u / wspeed        h_s / height
function run_plume( dem, source, stability::AbstractString, outdoor::AbstractString, resolution::Integer, wind_direction, concentration::Real, wind_speed::Real, stack_height::Real, 
                  # v_s / gspeed         d_s / diameter            t_s / temp                    t_a / etemp
                    gas_speed::Real=0.0, stack_diameter::Real=0.0, smoke_temperature::Real=0.0,  temperature::Real=0.0, output_path::AbstractString=".\\otput_model.tiff" )

  if agd.geomdim(source) != 0
      throw(DomainError(source, "`source` must be a point"))
  end

  refsys = agd.getspatialref(source)

  if agd.importWKT(agd.getproj(dem)) != refsys
      throw(DomainError("The reference systems are not uniform. Aborting analysis."))
  end
  
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

  toxic = Functions.substance_extract(contaminant, lst_fields, ".\\..\\Library\\")

  feature = collect(agd.getfeature(source))
  geom = agd.getgeom(feature[1])
  x_source = agd.getx(geom, 0)
  y_source = agd.gety(geom, 0)
  r_source, c_source = toCoords(dtm, x_source, y_source)

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

  rows = maxR - minR
  cols = maxC - minC
  minX, maxY = toCoords(dtm, minX, maxY)

  gtiff_driver = agd.getdriver("GTiff")
  target_ds = agd.create( path, gtiff_driver, rows, cols, 1, agd.GDAL.GDT_Float32 )
 # NON SONO CERTO CHE IL GEOTRASFORM VADA BENE
  agd.setgeotransform!( target_ds, [ minX, resolution, 0.0, maxY, 0.0, -resolution ] )
  agd.setproj!( target_ds, refsys )
 """ NON SO QUALE SIA IL COMANDO PER SETTARE I METADATI CON `ArchGDAL`
  target_ds.SetMetadata(
    Dict(
      "credits" => "Envifate - Francesco Geri, Oscar Cainelli, Paolo Zatelli, Gianluca Salogni, Marco Ciolli - DICAM Università degli Studi di Trento - Regione Veneto",
      "modulo" => "Dispersione in falda",
      "descrizione" => "Simulazione di dispersione inquinante in falda",
      "srs" => refsys,
      "data" => today()
    )
  )
 """
  valNoData = -9999.0
  band1 = agd.getband( target_ds, 1 )
  agd.setnodatavalue!(band1, valNoData)
  agd.fillraster!(band1, valNoData)
  band = agd.read(band1)

  for (point, value) in zip(points[i], values[i])
      r, c = point - (minR, minC)
      band[r, c] = value
  end
end



end # module