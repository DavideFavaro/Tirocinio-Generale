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

include("../Library/Functions.jl")


#=
mutable struct Plume
  """docstring for element"""
  # concentration: concentrazione inquinante m3/sec
  # d: distanza (coordinata x)
  # y,z coordinate in metri
  # stability: classe di stabilità
  # wind_speed: velocità del vento nella direzione x
  # stack_height: altezza camino
  # gas_speed: velocità gas
  # stack_diameter: diametro camino
  # smoke_temperature: temperatura gas all'uscita del camino
  # temperature: temperatura ambiente

  concentration::Real
  d::Integer
  y::Integer
  z::Integer
  stability::AbstractString
  wind_speed::Real
  stack_height::Real
  gas_speed::Real
  stack_diameter::Real
  smoke_temperature::Real
  temperature::Real
  x_w::Integer
  outdoor::AbstractString

  H
  σy
  σz
  g1
  g2

  Plume(concentration,d,y,z,stability,wind_speed,stack_height,gas_speed,stack_diameter,smoke_temperature,temperature,x_w,outdoor) = new(concentration,d,y,z,stability,outdoor,wind_speed,gas_speed,stack_diameter,smoke_temperature,temperature,stack_height,x_w)
end
=#

#= VERSIONE CON LO STRUCT
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
=#
function calc_h( d::Real, stack_diameter::Real, gas_speed::Real, smoke_temperature::Real, temperature::Real )
    try
      fb = 9.81 * ( (stack_diameter * gas_speed) / 4 ) * ( ( smoke_temperature / temperature ) / smoke_temperature )
      Δh = 1.6fb^(1/3) * d^(2/3)
      return stack_height + Δh
    catch
      return stack_height
    end
end


#= VERSIONE CON LO STRUCT
function calc_sigma!( p::Plume )
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
=#
function calc_σ( d::Real, stability_class::AbstractString, outdoor::AbstractString )
    σ_values = Functions.air_extract( stability_class, outdoor )
    σy = ( σ_values[0] * d ) / ( 1 + σ_values[1] * d )^σ_values[2]
    σz = ( σ_values[3] * d ) / ( 1 + σ_values[4] * d )^σ_values[5]
    return σy, σz
  end


#= VERSIONE CON LO STRUCT
function calc_g!( p::Plume )
  p.g1 = exp( ( -0.5 * p.y^2 ) / p.σy^2 )
  p.g2 = exp( ( -0.5 * (p.z - p.stack_height)^2 ) / p.σz^2 ) + exp( ( -0.5 * (p.z + p.stack_height)^2 ) / p.σz^2 )
  return p.g1, p.g2
end
=#
function calc_g( y::Real, z::Real, σy::Real, σz::Real, stack_height::Real )
    g1 = ℯ^( ( -0.5y^2 ) / σy^2 )
    g2 = ℯ^( ( -0.5(z - stack_height)^2 ) / σz^2 ) + ℯ^( ( -0.5(z + stack_height)^2 ) / σz^2 )
    return g1, g2
  end


#= VERSIONE CON LO STRUCT
function calc_C!( p::Plume )
  p.C = ( 100p.concentration / 3600p.wind_speed ) * ( (p.g1 * p.g2) / ( 2π * p.σy * p.σz ) )
  return p.C
end
=#
function calc_C( concentration::Real, wind_speed::Real, σy::Real, σz::Real, g1::Real, g2::Real )
    return ( 100concentration / 3600wind_speed ) * ( (g1 * g2) / ( 2π * σy * σz ) )
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



         #                                                                                                         q / text_conc        u / wspeed
function run_plume( dem, source, area, stability::AbstractString, outdoor::AbstractString, resolution::Int64, x_w, concentration::Real, wind_speed::Real, 
                  # h_s / height        v_s / gspeed         d_s / diameter            t_s / temp                    t_a / etemp
                    stack_height::Real, gas_speed::Real=0.0, stack_diameter::Real=0.0, smoke_temperature::Real=0.0,  temperature::Real=0.0, output_path::AbstractString=".\\otput_model.tiff" )

    if agd.geomdim(source) != 0
        throw(DomainError(source, "`source` must be a point"))
    end
    if agd.geomdim(area) != 2
        throw(DomainError(source, "`area` must be a polygon"))
    end

 """ CONTROLLO SULLA VALIDITA' DEL RASTER
    if not self.dem.isValid():
        QMessageBox.warning(self,"Warning", "Il file DEM non è valido" )
        return
 """
 
    if agd.getspatialref(area) != agd.getspatialref(source) || agd.getspatialref(dem) != agd.getspatialref(source)
        throw(DomainError("The reference systems are not uniform. Aborting analysis."))
    end


    # self.text_line_srid=str(self.line_srid.text())
    
    lst_fields = [ "sf_ing", "sf_inal", "iur", "rfd_ing", "rfd_inal", "rfc" ]

    lst_toxic_msg = [ 
        "Slope Factor per ingestione", 
        "Slope Factor per inalazione",
        "Inhalation Unit Risk",
        "Reference Dose per ingestione",
        "Reference Dose per inalazione",
        "Reference Concentration"
    ]

 # DA SISTEMARE IL PATH
    toxic = Functions.substance_extract( contaminant, lst_fields, os.path.dirname(__file__)+"/../library/" )


    # try:
    #     self.srid=int(self.text_line_srid)
    #     if self.srid not in self.list_srid :
    #         QMessageBox.warning(self,"Warning", u"Errore codice srid" )
    #         return
    # except Exception as e:
    #     QMessageBox.warning(self,"Warning", u"Errore codice srid" )
    #     return





    refsys = agd.importEPSG(agd.fromWKT(agd.getspatialref(source)))
 """ PRINT DI COSE
    messaggio="Inizio elaborazione plume atmosferico Envifate\n"
    messaggio+="---------------------------\n\n"
    messaggio+="FILE DI INPUT:\n"
    messaggio+="Vettoriale sorgente: "+str(self.text_vector)+"\n"
    messaggio+="Vettoriale confine: "+str(self.text_area)+"\n"
    messaggio+="DTM: "+str(self.text_dem)+"\n\n"
    messaggio+="VARIABILI:\n"
    messaggio+=u"Concentrazione inquinante: "+str(self.text_conc)+" Kg/h\n"
    messaggio+=u"Classe stabilità atmosferica: "+str(self.combo_stability.currentText())+"\n"
    messaggio+=u"Classe outdoor: "+str(self.class_outdoor)+"\n"
    messaggio+="Direzione vento: "+str(self.combo_maindirwind.currentText())+"\n"
    messaggio+=u"Velocità vento: "+str(self.text_wspeed)+" m/s\n"
    messaggio+="Temperatura: "+str(self.text_etemp)+"\n"
    messaggio+=u"Diametro stack: "+str(self.text_diameter)+" m\n"
    messaggio+=u"Altezza stack: "+str(self.text_height)+" m\n"
    messaggio+="Temperatura gas: "+str(self.text_temp)+" °C\n"
    messaggio+=u"Velocità gas: "+str(self.text_gspeed)+" m/s\n"
    messaggio+="Risoluzione: "+str(self.res)+"\n\n"
    messaggio+='ALGORITMO UTILIZZATO: Pasquill-Gifford (Pasquill, F., 1961: The estimation of the dispersion of windborne material. Meteor. Mag.,90, 33–49.; Gifford, F. A., Jr., 1961: Use of routine observations for estimating atmospheric dispersion. Nucl. Saf.,2, 47–57; Turner, D. B., 1967: Workbook of atmospheric dispersion estimates. PHS Publ. 999 AP-26, 84 pp.)\n\n'
    messaggio+="---------------------------\n\n"
 """

    area_layer = agd.getlayer(area, 0)
 # NON FUNZIONANTE / DA ELIMINARE
    x_min, y_min, x_max, y_max = agd.envelope(area_layer)
    valNoData = -9999
    # Create the destination data source
    x_res = ( x_max - x_min ) / resolution
    y_res = ( y_max - y_min ) / resolution

    
    gtiff_driver = agd.getdriver("GTiff")
    target_ds = agd.create( output_path, gtiff_driver, round(Int64, x_res), round(Int64, y_res), 1, agd.GDAL.GDT_Float32 )
    agd.setgeotransform!(target_ds, [ x_min, resolution, 0.0, y_max, 0.0, -resolution ])
    agd.setproj!(target_ds, refsys)
 """ NON SO QUALE SIA IL COMANDO PER SETTARE I METADATI CON `ArchGDAL`
    target_ds.SetMetadata(
        Dict(
            "credits" => "Envifate - Francesco Geri, Oscar Cainelli, Paolo Zatelli, Gianluca Salogni, Marco Ciolli - DICAM Università degli Studi di Trento - Regione Veneto",
            "modulo" => "Dispersione atmosferica",
            "descrizione" => "Simulazione di dispersione di un inquinante in atmosfera da una sorgente tipo "ciminiera", modello gaussiano di Pasquill-Gifford",
            "srs" => refsys,
            "data" => today()
        )
    )
 """
    band1 = agd.getband(target_ds, 1)
    agd.setnodatavalue!( band1, Float64(valNoData) )
    band = agd.read(band1)
    agd.fillraster!(band, valNoData)
    xsize = agd.width(band)
    ysize = agd.height(band)
 # NON SONO CERTO SIA IL METODO GIUSTO 
    # outData = deepcopy(band)
    outData = band

    polygons = collect(agd.getfeature(area))
    feature = collect(agd.getfeature(source))
    geom = agd.getgeom(feature[1])
    x_source = agd.getx(geom, 0)
    y_source = agd.gety(geom, 0)

    rows = ysize - 1
    cols = xsize - 1

    # start_time = time.time()

    max_dominio=0.0

    for row in 1:rows
        for col in 1:cols
            x, y = (col, row) * resolution + (x_min, y_min) .+ (resolution/2)
         # POTREBBE NON ESSERE COSI' SEMPLICE
            z = dem[x, y]
            Δx = x - x_source
            Δy = y - y_source
            true_x = Δx * cos(deg2rad(x_w)) -  Δy * sin(de2rad(x_w))
            true_y = Δx * sin(deg2rad(x_w)) +  Δy * cos(de2rad(x_w))

            if true_y > 0
             #= VERSIONE CON STRUCT `Plume`
                element = Plume(concentration, true_y, true_x, z[1], stability, outdoor, wind_speed, stack_height, gas_speed, stack_diameter, smoke_temperature, temperature, x_w)
                calc_sigma!(element)
                calc_g!(element)
                calc_h!(element)
                cfinal = calc_C!(element)
             =#
                σy, σz = calc_σ( true_y, stability, outdoor )
                g1, g2 = calc_g( true_y, z[1], σy, σz, stack_height )
             # IL CALCOLO QUI SOTTO SEMBRA ESSERE INUTILE
                s_height = calc_h( true_y, stack_diameter, gas_speed, smoke_temperature, temperature )
                cfinal = calc_C( concentration, wind_speed, σy, σz, g1, g2 )
                if cfinal > max_dominio
                    max_dominio = cfinal
                end
                outData[row, col] = cfinal
            else
                outData[row, col] = 0.0
            end
        end
    end







    outData_raster=outData[::-1]
    band.WriteArray(outData_raster)
    band= None
    target_ds = None

    base_raster_name=os.path.basename(self.path_output)
    raster_name=os.path.splitext(base_raster_name)[0]
    self.outputlayer=self.iface.addRasterLayer(self.path_output, raster_name)

    layer=None
    for lyr in list(QgsProject.instance().mapLayers().values()):
        if lyr.name() == raster_name:
            layer = lyr
        end
    end

    functions.applystyle(layer,'gr',0.5)


 """ PRINT DI COSE
    max_round_dominio = round(max_dominio*1000000,5)
  
    tempoanalisi=time.time() - start_time
    tempostimato=time.strftime("%H:%M:%S", time.gmtime(tempoanalisi))
    messaggio="---------------------------------\n"
    messaggio+="Fine modellazione\n"
    messaggio+="\nTempo di analisi: "+tempostimato+"\n"

    max_round_dominio=round(max_dominio*1000000,10)

    messaggio+="\nMassima concentrazione rilevata nel dominio di analisi: "+str(max_round_dominio)+" \u03BCg /  m\u00b3\n\n\n"

    check_alert_msg=0
    for idx, tox in enumerate(toxic):
        if tox:
            if float(tox)<=max_round_dominio:
                if check_alert_msg==0:
                    messaggio+="*** ALLERTA RILEVATA!!! ***\n\n"
                    check_alert_msg=1
                messaggio+="Concetrazione limite superata per: "+lst_toxic_msg[idx]+" (valore soglia: "+str(tox)+")\n\n"

    messaggio+="---------------------------------\n\n"
    self.console.appendPlainText(messaggio)

    self.label_status.setText("In attesa di dati")
    self.label_status.setStyleSheet('color : green; font-weight:bold')
"""



end # module