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
#=
    def run_sediment(self):


        self.text_v=str(self.tableWidget.item(3,0).text())
        self.text_h=str(self.tableWidget.item(4,0).text())
        self.text_dx=str(self.tableWidget.item(5,0).text())
        self.text_dy=str(self.tableWidget.item(6,0).text())
        self.text_q=str(self.tableWidget.item(7,0).text())
        self.text_w=str(self.tableWidget.item(8,0).text())

        self.dir=self.classiwind[self.combo_maindirwind.currentText()]

        self.text_t=str(self.tableWidget.item(9,0).text())
        self.text_dt=str(self.tableWidget.item(10,0).text())
        self.text_u=str(self.tableWidget.item(11,0).text())
        self.text_tide=str(self.tableWidget.item(12,0).text())


        self.text_vector = str(self.combo_source.currentText())
        self.text_area = str(self.combo_bound.currentText())


        self.res=int(self.spinRes.text())



        try:
            self.v=float(self.text_v)
        except Exception as e:
            QMessageBox.warning(self,"Warning", "La velocità media della corrente è obbligatoria" )
            return

        try:
            self.h=float(self.text_h)
        except Exception as e:
            QMessageBox.warning(self,"Warning", "La profondità media è obbligatoria" )
            return

        try:
            self.dx=float(self.text_dx)
        except Exception as e:
            QMessageBox.warning(self,"Warning", "Il coefficiente di dispersione sull'asse X è obbligatorio" )
            return

        try:
            self.dy=float(self.text_dy)
        except Exception as e:
            QMessageBox.warning(self,"Warning", "Il coefficiente di dispersione sull'asse Y è obbligatorio" )
            return

        try:
            self.q=float(self.text_q)
        except Exception as e:
            QMessageBox.warning(self,"Warning", "Il quantitativo di materia dragata è obbligatorio")
            return

        try:
            self.w=float(self.text_w)
        except Exception as e:
            QMessageBox.warning(self,"Warning", "La velocità media di sedimentazione marina è obbligatoria")
            return

        # try:
        #     self.dir=int(self.text_dir)
        # except Exception as e:
        #     QMessageBox.warning(self,"Warning", "La direzione media della corrente è obbligatoria")
        #     return


        try:
            self.time=int(self.text_t)
        except Exception as e:
            QMessageBox.warning(self,"Warning", "Il tempo di analisi è obbligatorio")
            return


        try:
            self.dt=int(self.text_dt)
        except Exception as e:
            QMessageBox.warning(self,"Warning", "L'intervallo di analisi è obbligatorio")
            return

        if self.text_u!='':
            try:
                self.u=int(self.text_u)
            except Exception as e:
                QMessageBox.warning(self,"Warning", "Errore nel dato di ampiezza dell'oscillazione della corrente" )
                return
        else:
            self.u=0

        if self.text_tide!='':
            try:
                self.tide=int(self.text_tide)
            except Exception as e:
                QMessageBox.warning(self,"Warning", "Errore nel dato di ciclo di marea" )
                return
        else:
            self.tide=0

        self.source=self.listalayers[self.text_vector]

        if self.source.wkbType()!=1:
            QMessageBox.warning(self,"Warning", "The source file must have point geometry" )
            return

        self.areastudio=self.listalayers[self.text_area]


        if self.areastudio.wkbType()!=6:
            QMessageBox.warning(self,"Warning", "The boundaries file must have polygon geometry" )
            return


        self.path_output=self.line_output.text()
        if self.path_output=="":
            self.path_output=os.path.dirname(__file__)+"/output_model.tif"


        if self.areastudio.crs().authid()!=self.source.crs().authid():
            QMessageBox.warning(self,"Warning", "Errore: i sistemi di riferimento non sono uniformi. Impossibile continuare con l'analisi." )
            return

        self.refsys=self.source.crs().authid().split(':')[1]

        messaggio="Inizio elaborazione plume atmosferico Envifate\n"
        messaggio+="---------------------------\n\n"
        messaggio+="FILE DI INPUT:\n"
        messaggio+="Vettoriale sorgente: "+str(self.text_vector)+"\n"
        messaggio+="Vettoriale confine: "+str(self.text_area)+"\n"
        messaggio+="VARIABILI:\n"
        messaggio+=u"Quantità di sedimento dragata: "+str(self.text_q)+"\n"
        messaggio+=u"Profondità media: "+str(self.text_h)+"\n"
        messaggio+=u"Velocità media della corrente: "+str(self.text_v)+"\n"
        messaggio+=u"Direzione media della corrente: "+str(self.combo_maindirwind.currentText())+"\n"
        messaggio+=u"Coefficienti di diffusione (x,y): "+str(self.text_dx)+" "+str(self.text_dy)+"\n"
        messaggio+=u"Velocità di sedimentazione marina: "+str(self.text_w)+"\n"
        messaggio+=u"Tempo di analisi: "+str(self.text_t)+"\n"
        messaggio+=u"Intervallo di analisi: "+str(self.text_dt)+"\n"
        if self.text_u!="":
            messaggio+=u"Ampiezza della marea: "+str(self.text_u)+"\n"
        if self.text_tide!="":
            messaggio+=u"Ciclo della marea: "+str(self.text_tide)+"\n"
        messaggio+="Risoluzione: "+str(self.res)+"\n\n"
        messaggio+='ALGORITMO UTILIZZATO: Shao (Shao, Dongdong, et al. "Modeling dredging-induced turbidity plumes in the far field under oscillatory tidal currents." Journal of Waterway, Port, Coastal, and Ocean Engineering 143.3 (2016))\n\n'
        messaggio+="---------------------------\n\n"

        self.console.appendPlainText(messaggio)


        self.label_status.setText("Preparazione dati")
        self.label_status.setStyleSheet('color : #e8b445;font-weight:bold')




        path_layer=self.areastudio.dataProvider().dataSourceUri()
        path=path_layer.split("|")
        source_ds = ogr.Open(path[0])
        area_layer = source_ds.GetLayer()
        #x_min, x_max, y_min, y_max = area_layer.GetExtent()
        x_min=int(area_layer.GetExtent()[0])
        y_min=int(area_layer.GetExtent()[2])
        x_max=int(area_layer.GetExtent()[1])
        y_max=int(area_layer.GetExtent()[3])

        drivermem = gdal.GetDriverByName('MEM')
        # Define pixel_size and NoData value of new raster
        pixel_size = self.res
        NoData_value = -9999



        #metadata
        #source_ds.SetMetadata({'descrizione':'Simulazione di sedimento disperso in ambiente marino','srs':self.source.crs().authid()})



        #fine metadata

        # pyqtRemoveInputHook()
        # pdb.set_trace()
        # Create the destination data source
        x_res = int((x_max - x_min) / pixel_size)
        y_res = int((y_max - y_min) / pixel_size)
        #target_ds = drivermem.Create('', x_res, y_res, 1, gdal.GDT_Byte)
        target_ds = gdal.GetDriverByName('GTiff').Create(self.path_output, x_res, y_res, 1, gdal.GDT_Float32)
        target_ds.SetGeoTransform((x_min, pixel_size, 0, y_max, 0, -pixel_size))

        projectionfrom = target_ds.GetProjection()
        #geotransform = target_ds.GetGeoTransform()
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(int(self.refsys))
        target_ds.SetProjection( srs.ExportToWkt() )
        target_ds.SetMetadata({'credits':'Envifate - Francesco Geri, Oscar Cainelli, Paolo Zatelli, Gianluca Salogni, Marco Ciolli - DICAM Università degli Studi di Trento - Regione Veneto',
                               'modulo':'Analisi sedimentazione marina',
                               'descrizione':'Simulazione di sedimento disperso in ambiente marino',
                               'srs':self.source.crs().authid(),
                               'data':datetime.datetime.now().strftime("%d-%m-%y")})


        band = target_ds.GetRasterBand(1)
        band.SetNoDataValue(float(NoData_value))
        band.Fill(NoData_value)
        xsize = band.XSize
        ysize = band.YSize

        outData = np.array(band.ReadAsArray(0, 0, xsize,ysize).astype(np.float))

        feature = next(self.source.getFeatures())
        geom = feature.geometry().asPoint()
        x_source=geom[0]
        y_source=geom[1]

        rows=ysize-1
        cols=xsize-1


        max_progress=rows*cols
        self.progressBar.setMaximum(max_progress)
        start_time = time.time()

        self.label_status.setText("Processing data")
        self.label_status.setStyleSheet('color : #e8b445;font-weight:bold')



        index_progress=0
        controllo=1
        if controllo==1:
            #values = []
            for row in range(rows):
                #row_data = []
                for col in range(cols):
                    index_progress+=1
                    self.progressBar.setValue(index_progress)
                    x = col*pixel_size+x_min+(pixel_size/2)
                    y = row*pixel_size+y_min+(pixel_size/2)

                    deltax=x-x_source
                    deltay=y-y_source
                    xvero=deltax*math.cos(math.radians(self.dir))-deltay*math.sin(math.radians(self.dir))
                    yvero=deltax*math.sin(math.radians(self.dir))+deltay*math.cos(math.radians(self.dir))


                    # pyqtRemoveInputHook()
                    # pdb.set_trace()


                    #controllo=1
                    #if controllo==1:
                    if yvero>0:

                        element=sediment.sediment(self.time,self.q,self.h,self.dx,self.dy,yvero,xvero,self.v,self.w,self.dt,self.u,self.tide)

                        # pyqtRemoveInputHook()
                        # pdb.set_trace()

                        q=element.calc_q()
                        n=int(element.t/element.dt)
                        csum=0
                        for i in range(n):
                            e=element.calcolo_e(i)
                            csum+=e*(1 / (element.t-i*element.dt) )
                        outData[row,col]=q*csum*element.dt
                    else:
                        outData[row,col]=0




            self.label_status.setText("Preparazione output")
            self.label_status.setStyleSheet('color : #e8b445;font-weight:bold')

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



        #functions.applystyle(layer,'gr',0.5)


        tempoanalisi=time.time() - start_time
        tempostimato=time.strftime("%H:%M:%S", time.gmtime(tempoanalisi))
        messaggio="---------------------------------\n"
        messaggio+="Fine modellazione\n"
        messaggio+="\nTempo di analisi: "+tempostimato+"\n"
        messaggio+="---------------------------------\n\n"
        self.console.appendPlainText(messaggio)

        self.label_status.setText("In attesa di dati")
        self.label_status.setStyleSheet('color : green; font-weight:bold')

        ax1f1 = self.figure.gca(projection='3d')
        (x, y) = np.meshgrid(outData_raster.shape[0], outData_raster.shape[1])
        #fig = plt.figure()
        #ax1f1 = fig.add_subplot(111, projection='3d')
        surf = ax1f1.plot_wireframe(x, y, outData_raster)
        #ax1f1.plot(self.list_result)

        self.canvas_mat.draw()
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
