module Sediments

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
using Parameters
using Dates

include("../Library/Functions.jl")



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
Recursively compute the concentration of each point and add the value and its indexes to positions
"""
function expand!( positions::AbstractVector, results::AbstractVector, dtm, indx_x::Integer, indx_y::Integer, sediment::Sediment )
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
    sediment.y = Δy * sindir + Δy * cosdir
    sediment.x = Δx * cosdir - Δy * sindir
    concentration = calcSediment!(sediment)
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
    def help(self):
        #self.credits = u"Università della Tuscia\n Viterbo - Italy\nRaffaele Pelorosso, Federica Gobattoni\nDeveloper: Francesco Geri"
        #QMessageBox.about(self.dlg,"Credits", self.credits )
        if platform.uname()[0]=="Windows":
            os.system("start "+os.path.dirname(__file__)+"/../tutorial/manuale_envifate_sedimentazione_marina.pdf")
        if platform.uname()[0]=="Linux":
            os.system("xdg-open "+os.path.dirname(__file__)+"/../tutorial/manuale_envifate_sedimentazione_marina.pdf")
        else:
            # pyqtRemoveInputHook()
            # pdb.set_trace()
            os.system("open "+os.path.join(os.path.dirname(__file__), "../tutorial/manuale_envifate_sedimentazione_marina.pdf"))


    def popolacombo(self):
        self.progressBar.setValue(0)
        self.combo_source.clear()
        self.combo_bound.clear()
        self.combo_maindirwind.clear()
        self.line_output.clear()

        self.allLayers = self.canvas.layers()
        self.listalayers=dict()
        elementovuoto="No file"

        for i in self.allLayers:
            self.listalayers[i.name()]=i
            if i.type() == QgsMapLayer.VectorLayer:
                self.combo_source.addItem(str(i.name()))
                self.combo_bound.addItem(str(i.name()))


        self.combo_maindirwind.addItem("N")
        self.combo_maindirwind.addItem("NE")
        self.combo_maindirwind.addItem("E")
        self.combo_maindirwind.addItem("SE")
        self.combo_maindirwind.addItem("S")
        self.combo_maindirwind.addItem("SW")
        self.combo_maindirwind.addItem("W")
        self.combo_maindirwind.addItem("NW")
=#


                    #                                     v                      h                 dx                        dy                        q                   dir
function run_sediment( source, resolution::Integer, mean_flow_speed::Real, mean_depth::Real, x_dispersion_coeff::Real, y_dispersion_coeff::Real, dredged_mass::Real, flow_direction::Real,
                    #  w                                  t / time       dt                      u
                       mean_sedimentation_velocity::Real, time::Integer, time_intreval::Integer, current_oscillatory_amplitude::Integer=0, tide::Integer=0, output_path::AbstractString=".\\output_model.tiff" )

    if agd.geomdim(source) != 0
        throw(DomainError(source, "`source` must be a point"))
    end
 
    refsys = agd.importEPSG(agd.fromWKT(agd.getspatialref(source)))


 # SERVE UN ELEMENTO DI classwind
 #   self.dir=self.classiwind[self.combo_maindirwind.currentText()]

 # messaggio+='ALGORITMO UTILIZZATO: Shao (Shao, Dongdong, et al. "Modeling dredging-induced turbidity plumes in the far field under oscillatory tidal currents." Journal of Waterway, Port, Coastal, and Ocean Engineering 143.3 (2016))\n\n'



    feature = collect(agd.getfeature(source))
    geom = agd.getgeom(feature[1])
    x_source = agd.getx(geom, 0)
    y_source = agd.gety(geom, 0)
    r_source, c_source = toCoords(dtm, x_source, y_source)

    #   start_time = time.time()


    points = [ (r_source, c_source) ]
 # NON CREDO SIA IL VALORE CORRETTO DA INSERIRE
    values = [ dredged_mass ]
    element = Sediment( dredged_mass, time, mean_depth, x_dispersion_coeff, y_dispersion_coeff, 0.0, 0.0, mean_flow_speed, mean_sedimentation_velocity, time_intreval, current_oscillatory_amplitude, tide)
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
