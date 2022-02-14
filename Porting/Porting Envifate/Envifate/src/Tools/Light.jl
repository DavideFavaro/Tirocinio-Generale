module Lights

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
#=
    def help(self):
        #self.credits = u"Università della Tuscia\n Viterbo - Italy\nRaffaele Pelorosso, Federica Gobattoni\nDeveloper: Francesco Geri"
        #QMessageBox.about(self.dlg,"Credits", self.credits )
        if platform.uname()[0]=="Windows":
            os.system("start "+os.path.dirname(__file__)+"/../tutorial/manuale_envifate_inquinamento_luminoso.pdf")
        if platform.uname()[0]=="Linux":
            os.system("xdg-open "+os.path.dirname(__file__)+"/../tutorial/manuale_envifate_inquinamento_luminoso.pdf")
        else:
            os.system("open "+os.path.dirname(__file__)+"/../tutorial/manuale_envifate_inquinamento_luminoso.pdf")
=#

import ArchGDAL as agd

include("..\\Library\\Functions.jl")


#=
def run_light(self):

        self.text_line_intensity=str(self.tableWidget.item(3,0).text())
        self.text_line_hsource=str(self.tableWidget.item(3,0).text())
        self.text_line_htarget=str(self.tableWidget.item(4,0).text())
        self.text_line_rarefraction=str(self.tableWidget.item(5,0).text())
        self.text_line_memory=str(self.tableWidget.item(6,0).text())


        self.text_vector = str(self.combo_source.currentText())
        self.text_area = str(self.combo_bound.currentText())
        self.text_dem = str(self.combo_dem.currentText())


        self.res=int(self.spinRes.text())



        if self.text_line_hsource!='':
            try:
                self.hsource=float(self.text_line_hsource)
            except Exception as e:
                QMessageBox.warning(self,"Warning", "Errore nell'altezza della sorgente" )
                return
        else:
            self.hsource=0.00

        if self.text_line_htarget!='':
            try:
                self.htarget=float(self.text_line_htarget)
            except Exception as e:
                QMessageBox.warning(self,"Warning", "Errore nell'altezza dell'osservatore" )
                return
        else:
            self.htarget=1.75


        if self.text_line_rarefraction!='':
            try:
                self.rarefraction=float(self.text_line_rarefraction)
            except Exception as e:
                QMessageBox.warning(self,"Warning", "Errore nel coefficiente di rarefrazione" )
                return
        else:
            self.rarefraction=0.14286


        if self.text_line_memory!='':
            try:
                self.memory=int(self.text_line_memory)
            except Exception as e:
                QMessageBox.warning(self,"Warning", "Errore nell'indicazione della memoria di processing" )
                return
        else:
            self.memory=500


        # wkbType: 1:point, 6:multipolygon, 2: Linestring

        self.dem=self.listalayers[self.text_dem]

        if not self.dem.isValid():
            QMessageBox.warning(self,"Warning", "The dem file is not valid" )
            return

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
            self.path_output=os.path.dirname(__file__)+"/light_intensity.tif"


        if self.areastudio.crs().authid()!=self.source.crs().authid() or self.dem.crs().authid()!=self.source.crs().authid():
            QMessageBox.warning(self,"Warning", "Errore: i sistemi di riferimento non sono uniformi. Impossibile continuare con l'analisi." )
            return

        self.refsys=self.source.crs().authid().split(':')[1]


        #recupero dati database


        messaggio="Inizio elaborazione analisi inquinamento luminoso\n"
        messaggio+="---------------------------\n\n"
        messaggio+="FILE DI INPUT:\n"
        messaggio+="Vettoriale sorgente: "+str(self.text_vector)+"\n"
        messaggio+="Vettoriale confine: "+str(self.text_area)+"\n"
        messaggio+="DTM: "+str(self.text_dem)+"\n\n"

        messaggio+="VARIABILI:\n"
        messaggio+="Altezza sorgente: "+str(self.text_line_hsource)+" m\n"
        messaggio+="Altezza osservatore: "+str(self.text_line_htarget)+" m\n"
        messaggio+="Coefficiente rarefrazione: "+str(self.rarefraction)+"\n"
        messaggio+="Risoluzione: "+str(self.res)+"\n\n"
        messaggio+='ALGORITMO UTILIZZATO: decadimento dell\'onda luminosa in funzione della distanza; analisi di intervisibilità r.viewshed\n\n'
        messaggio+="---------------------------\n\n"
        self.console.appendPlainText(messaggio)


        self.label_status.setText("Preparazione dati")
        self.label_status.setStyleSheet('color : #e8b445;font-weight:bold')



        path_layer=self.areastudio.dataProvider().dataSourceUri()
        path=path_layer.split("|")
        source_ds = ogr.Open(path[0])
        area_layer = source_ds.GetLayer()
        x_min=int(area_layer.GetExtent()[0])
        y_min=int(area_layer.GetExtent()[2])
        x_max=int(area_layer.GetExtent()[1])
        y_max=int(area_layer.GetExtent()[3])

        drivermem = gdal.GetDriverByName('MEM')
        pixel_size = self.res
        NoData_value = -9999

        # Create the destination data source
        x_res = (x_max - x_min) / pixel_size
        y_res = (y_max - y_min) / pixel_size

        target_ds = gdal.GetDriverByName('GTiff').Create(self.path_output, int(x_res), int(y_res), 1, gdal.GDT_Float32)
        target_ds.SetGeoTransform((x_min, pixel_size, 0, y_max, 0, -pixel_size))
        projectionfrom = target_ds.GetProjection()

        srs = osr.SpatialReference()
        srs.ImportFromEPSG(int(self.refsys))
        target_ds.SetProjection( srs.ExportToWkt() )

        target_ds.SetMetadata({'credits':'Envifate - Francesco Geri, Oscar Cainelli, Paolo Zatelli, Gianluca Salogni, Marco Ciolli - DICAM Università degli Studi di Trento - Regione Veneto',
                               'modulo':'Analisi inquinamento luminoso',
                               'descrizione':'Simulazione di inquinamento luminoso da sorgente puntuale singola o multipla',
                               'srs':self.source.crs().authid(),
                               'data':datetime.datetime.now().strftime("%d-%m-%y")})
        # geotransform = target_ds.GetGeoTransform()


        band = target_ds.GetRasterBand(1)
        band.SetNoDataValue(float(NoData_value))
        band.Fill(NoData_value)
        xsize = band.XSize
        ysize = band.YSize


        outData = np.array(band.ReadAsArray(0, 0, xsize,ysize).astype(np.float))


        nfeature=0

        features=self.source.getFeatures()

        for feature in features:

            geom = feature.geometry().asPoint()
            x_source=geom[0]
            y_source=geom[1]


            idxlevel = self.source.fields().indexFromName('level')
            intensity=feature.attributes()[idxlevel]

            nfeature+=1

            polygons = [feature for feature in self.areastudio.getFeatures()]

            rows=ysize-1
            cols=xsize-1

            max_progress=rows*cols
            self.progressBar.setMaximum(max_progress)
            start_time = time.time()


            self.label_status.setText("Processing data")
            self.label_status.setStyleSheet('color : #e8b445;font-weight:bold')

            grass_area=str(x_min)+','+str(x_max)+','+str(y_min)+','+str(y_max)+' ['+str(self.areastudio.crs().authid())+']'
            grass_coord=str(x_source)+','+str(y_source)+' ['+str(self.source.crs().authid())+']'

            nameviewshed='viewshedanalysis'+str(nfeature)

            params = { '-b' : True, '-c' : False, '-e' : False, '-r' : False, 'GRASS_RASTER_FORMAT_META' : '', 'GRASS_RASTER_FORMAT_OPT' : '',
                       'GRASS_REGION_CELLSIZE_PARAMETER' : 0, 'GRASS_REGION_PARAMETER' :grass_area, 'coordinates' : grass_coord,
                       'input' : self.dem.dataProvider().dataSourceUri(), 'max_distance' : -1, 'memory' : self.memory, 'observer_elevation' : self.hsource, 'output' : nameviewshed,
                       'refraction_coeff' : self.rarefraction, 'target_elevation' : self.htarget }
            viewshed_proc = processing.run('grass7:r.viewshed', params)
            #viewshed=viewshed_proc['output']


            #QgsProject.instance().mapLayersByName("memory:viewshed")[0]




            #aggiungo per controllo la viewshed alla toc
            iface.addRasterLayer(viewshed_proc['output'])

            viewshed=QgsProject.instance().mapLayersByName(nameviewshed)[0]


            index_progress=0
            controllo=1
            if controllo==1:

                for row in range(rows):

                    for col in range(cols):
                        index_progress+=1
                        self.progressBar.setValue(index_progress)
                        x = col*pixel_size+x_min+(pixel_size/2)
                        y = row*pixel_size+y_min+(pixel_size/2)


                        punto_controllo = QgsPointXY(x,y)



                        for pol in polygons:
                            poly = pol.geometry()
                            if poly.contains(punto_controllo):

                                cfr_viewshed=viewshed.dataProvider().identify(QgsPointXY(x, y),QgsRaster.IdentifyFormatValue)


                                if cfr_viewshed.results()[1]==1:
                                    deltax=x-x_source
                                    deltay=y-y_source
                                    dist=math.sqrt(math.pow(deltay,2)+math.pow(deltax,2))


                                    #new_intensity=(1/math.pow(dist,2))*intensity
                                    new_intensity=(1/dist)*intensity


                                    if nfeature==1:
                                        if new_intensity>0:
                                            outData[row,col]=new_intensity
                                        else:
                                            outData[row,col]=0
                                    else:

                                        if new_intensity>0:
                                            outData[row,col]=outData[row,col]+new_intensity
                                        else:
                                            outData[row,col]=outData[row,col]

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
        self.iface.addRasterLayer(self.path_output, raster_name)


        layer=None
        for lyr in list(QgsProject.instance().mapLayers().values()):
            if lyr.name() == raster_name:
                layer = lyr


        functions.applystyle(layer,'gr',0.5)


        tempoanalisi=time.time() - start_time
        tempostimato=time.strftime("%H:%M:%S", time.gmtime(tempoanalisi))
        messaggio="---------------------------------\n"
        messaggio+="Fine modellazione\n"
        messaggio+="\nTempo di analisi: "+tempostimato+"\n"
        messaggio+="---------------------------------\n\n"
        self.console.appendPlainText(messaggio)

        self.label_status.setText("In attesa di dati")
        self.label_status.setStyleSheet('color : green; font-weight:bold')

=#



# === VIEWSHED TEST ===========================================================================================================================

profile1 = [
    (0,12),
    (25,12),
    (50,8),
    (75,10),
    (100,13),
    (125,11),
    (150,15),
    (175,14),
    (200,23)
]
profile2 = [
    (0,500),
    (25,12),
    (50,8),
    (75,10),
    (100,9),
    (125,11),
    (150,15),
    (175,14),
    (200,23)
]
profile3 = [
    (0,12),
    (25,12),
    (50,8),
    (75,10),
    (100,16),
    (125,11),
    (150,500),
    (175,14),
    (200,23)
]
profile = profile3

visible = [ profile[1], profile[2] ]
vali = abs( ( profile[2][2] - profile[1][2] ) / ( profile[2][1] - profile[1][1] ) )
for i in 3:length(profile)
    println("vali: $vali")
    val = ( profile[i][2] - profile[i-1][2] ) / ( profile[i][1] - profile[i-1][1] )
    print("$(profile[i]): $val")
    if val >= abs(vali)
        print("    PUSH")
        push!( visible, profile[i] )
    end
    println("\n")
    vali = vali < 0 ? vali+val : vali-val
end
visible

# ============================================================================================================================================




function veiwshed( profile::AbstractVector, result::Symbol=:visible )::AbstractVector
    if result ∉ [:visible, :invisible]
        throw(DomainError(result, "result must either be :visible or :invisible"))
    end
    if result == :visible
        visible = [ profile[1], profile[2] ]
        slope = ( profile[2][2] - profile[1][2] ) / ( profile[2][1] - profile[1][1] )
        for i in 3:length(profile)
            new_slope = ( profile[i][2] - profile[1][2] ) / ( profile[i][1] - profile[1][1] )
            if new_slope >= slope
                push!( visible, profile[i] )
                slope = new_slope
            end
        end
        return visible
    else
        invisible = []
        for i in 3:length(profile)
            new_slope = ( profile[i][2] - profile[1][2] ) / ( profile[i][1] - profile[1][1] )
            if new_slope < slope
                push!( nonvisible, profile[i] )
                slope = new_slope
            end
        end
        return invisible
    end
end



function run_light( dem, source, resolution::Integer, intenisty::Real, source_height::Real=0.0, observer_height::Real=1.75, rarefraction::Real=0.14286,
                    output_path::AbstractString=".\\light_intensity.tiff"  )

 """ CONTROLLO SUL RASTER dem
    if not self.dem.isValid():
        QMessageBox.warning(self,"Warning", "The dem file is not valid" )
        return
 """
    if intensity < 0
        throw(DomainError(intenisty, "`intenisty` must be positive"))
    end

    src_layer = agd.getlayer(source, 0)
    src_geoms = agd.getgeom.(collect(src_layer))

    if any( geom -> agd.geomdim(geom) != 0, src_geoms )
        throw(DomainError(source, "`source` must be a point"))
    end

    refsys = agd.getspatialref(src_layer)

    if agd.importWKT(agd.getproj(dem)) != refsys
        throw(DomainError("The reference systems are not uniform. Aborting analysis."))
    end

 """ PRINT DI COSE
    #recupero dati database
    messaggio="Inizio elaborazione analisi inquinamento luminoso\n"
    messaggio+="---------------------------\n\n"
    messaggio+="FILE DI INPUT:\n"
    messaggio+="Vettoriale sorgente: "+str(self.text_vector)+"\n"
    messaggio+="Vettoriale confine: "+str(self.text_area)+"\n"
    messaggio+="DTM: "+str(self.text_dem)+"\n\n"

    messaggio+="VARIABILI:\n"
    messaggio+="Altezza sorgente: "+str(self.text_line_hsource)+" m\n"
    messaggio+="Altezza osservatore: "+str(self.text_line_htarget)+" m\n"
    messaggio+="Coefficiente rarefrazione: "+str(self.rarefraction)+"\n"
    messaggio+="Risoluzione: "+str(self.res)+"\n\n"
    messaggio+='ALGORITMO UTILIZZATO: decadimento dell\'onda luminosa in funzione della distanza; analisi di intervisibilità r.viewshed\n\n'
    messaggio+="---------------------------\n\n"
    self.console.appendPlainText(messaggio)
 """

    nfeature = 0
    features = agd.getgeom.( agd.getlayer(source, 0) )
    for feature in features
        x_source = agd.getx(feature, 0)
        y_source = agd.gety(feature, 0)


     # LA/E INTENSITA' LA/E INSERIAMO COME PARAMETRO (SINGOLO VALORE/ARRAY)
        #   idxlevel = self.source.fields().indexFromName('level')
        #   intensity=feature.attributes()[idxlevel]

        nfeature += 1

        #   start_time = time.time()


        # BISOGNA TROVARE IL MODO DI CALCOLARE IL RANGE MASSIMO POSSIBILE DELLA LUCE
        max_points = ?
        profiles = Viewshed.generate_profiles( max_points, x_source, y_source )
        viewshed = Viewshed.viewshed(profiles)

     """ SERVE LA VIEWSHED COME RASTER CREDO """
        viewshed=QgsProject.instance().mapLayersByName(nameviewshed)[0]
     """                                     """

     # QUESTO O QUALCHE ALTRO MODO PER OTTENERE I LIMITI DEL RASTER VIEWSHED
        # Define the area of the `dtm` that coincides with the area of the `viewshed`
        minX, maxY, maxX, maxY = getSidesDistances()
        row_begin, col_begin = toIndexes(dtm, minX, minY)
        row_end, col_end = toIndexes(dtm, maxX, maxY)



        valNoData = -9999
        gtiff_driver = agd.getdriver("GTiff")
        target_ds = agd.create( output_path, gtiff_driver, row_end-row_begin, col_end-col_begin, 1, agd.GDAL.GDT_Float32 )
        agd.setgeotransform!(target_ds, [ minX, resolution, 0.0, maxY, 0.0, -resolution ])
        agd.setproj!( target_ds, refsys )
     """ NON SO QUALE SIA IL COMANDO PER SETTARE I METADATI CON `ArchGDAL`
        target_ds.SetMetadata(
            Dict(
                "credits" => "Envifate - Francesco Geri, Oscar Cainelli, Paolo Zatelli, Gianluca Salogni, Marco Ciolli - DICAM Università degli Studi di Trento - Regione Veneto",
                "modulo" => "Analisi inquinamento luminoso",
                "descrizione" => "Simulazione di inquinamento luminoso da sorgente puntuale singola o multipla",
                "srs" => refsys,
                "data" => today()
            )
        )
     """
        band1 = agd.getband(target_ds, 1)
        agd.setnodatavalue!( band1, Float64(valNoData) )
        agd.fillraster!(band1, valNoData)
        band = agd.read(band1) 

        
        # Compute the light intensity in said area of the dtm
        for row in row_begin:row_end
            for col in col_begin:col_end
                if viewshed[row, col] == 1
                    x, y = toCoords(dtm, row, col)
                    dist = √( (y - y_source)^2 + (x - x_source)^2 )
                    new_intensity = intensity / dist

                    res = new_intensity > 0 ? new_intensity : 0
                    if nfeature == 1
                        band[row, col] = res
                    else
                        band[row, col] += res
                    end
                else
                    band[row, col] = 0
                end
            end
        end
    end


    # MANCA SCRIVERE IL RASTER




end



end # module





















function run_light( dem, source, resolution::Integer, srource_height::Real=0.0, observer_height::Real=1.75, rarefraction::Real=0.14286, processing_memory::Integer=500, output_path::AbstractString=".\\light_intensity.tiff"  )

    if agd.geomdim(source) != 0
        throw(DomainError(source, "`source` must be a point"))
    end

    refsys = agd.getspatialref(source)

    if agd.importWKT(agd.getproj(dem)) != refsys
        throw(DomainError("The reference systems are not uniform. Aborting analysis."))
    end


    nfeature = 0
    features = agd.getgeom.( agd.getfeature(source) )
    for feature in features
        x_source = agd.getx(geom, 0)
        y_source = agd.gety(geom, 0)


        idxlevel = self.source.fields().indexFromName('level')
        intensity=feature.attributes()[idxlevel]

        nfeature += 1


        
        


        grass_area=str(x_min)+','+str(x_max)+','+str(y_min)+','+str(y_max)+' ['+str(self.areastudio.crs().authid())+']'
        grass_coord=str(x_source)+','+str(y_source)+' ['+str(self.source.crs().authid())+']'

        nameviewshed='viewshedanalysis'+str(nfeature)

        params = { '-b' : True, '-c' : False, '-e' : False, '-r' : False, 'GRASS_RASTER_FORMAT_META' : '', 'GRASS_RASTER_FORMAT_OPT' : '',
                   'GRASS_REGION_CELLSIZE_PARAMETER' : 0, 'GRASS_REGION_PARAMETER' :grass_area, 'coordinates' : grass_coord,
                   'input' : self.dem.dataProvider().dataSourceUri(), 'max_distance' : -1, 'memory' : self.memory, 'observer_elevation' : self.hsource, 'output' : nameviewshed,
                   'refraction_coeff' : self.rarefraction, 'target_elevation' : self.htarget }
        viewshed_proc = processing.run('grass7:r.viewshed', params)




        #aggiungo per controllo la viewshed alla toc
        iface.addRasterLayer(viewshed_proc['output'])

        viewshed=QgsProject.instance().mapLayersByName(nameviewshed)[0]


        for row in 1:rows
            for col in 1:cols
                x = col*pixel_size+x_min+(pixel_size/2)
                y = row*pixel_size+y_min+(pixel_size/2)

                punto_controllo = QgsPointXY(x,y)

                for pol in polygons
                    poly = pol.geometry()
                    if poly.contains(punto_controllo)

                        cfr_viewshed=viewshed.dataProvider().identify(QgsPointXY(x, y),QgsRaster.IdentifyFormatValue)

                        if cfr_viewshed.results()[1] == 1
                            deltax=x-x_source
                            deltay=y-y_source
                            dist=math.sqrt(math.pow(deltay,2)+math.pow(deltax,2))

                            #new_intensity=(1/math.pow(dist,2))*intensity
                            new_intensity=(1/dist)*intensity

                            if nfeature == 1
                                if new_intensity > 0
                                    outData[row,col]=new_intensity
                                else
                                    outData[row,col]=0
                                end
                            else
                                if new_intensity > 0
                                    outData[row,col]=outData[row,col]+new_intensity
                                else
                                    outData[row,col]=outData[row,col]
                                end
                            end
                        else
                            outData[row,col]=0
                        end
                    end
                end
            end
        end
    end



    valNoData = -9999

    gtiff_driver = agd.getdriver("GTiff")
    target_ds = agd.create( output_path, gtiff_driver, round(Int64, x_res), round(Int64, y_res), 1, agd.GDAL.GDT_Float32 )
    agd.setgeotransform!(target_ds, [ x_min, resolution, 0.0, y_max, 0.0, -resolution ])
    agd.setproj!( target_ds, refsys )
 """ NON SO QUALE SIA IL COMANDO PER SETTARE I METADATI CON `ArchGDAL`
    target_ds.SetMetadata(
        Dict(
            "credits" => "Envifate - Francesco Geri, Oscar Cainelli, Paolo Zatelli, Gianluca Salogni, Marco Ciolli - DICAM Università degli Studi di Trento - Regione Veneto",
            "modulo" => "Analisi inquinamento luminoso",
            "descrizione" => "Simulazione di inquinamento luminoso da sorgente puntuale singola o multipla",
            "srs" => refsys,
            "data" => today()
        )
    )
 """
    band1 = agd.getband(target_ds, 1)
    agd.setnodatavalue!( band1, Float64(valNoData) )
    agd.fillraster!(band1, valNoData)
    xsize = agd.width(band1)
    ysize = agd.height(band1)
    band = agd.read(band1)
    # outData = deepcopy(band)
    outData = band
    rows=ysize-1
    cols=xsize-1


        self.label_status.setText("Preparazione output")
        self.label_status.setStyleSheet('color : #e8b445;font-weight:bold')

        outData_raster=outData[::-1]
        band.WriteArray(outData_raster)


    band= None
    target_ds = None



    base_raster_name=os.path.basename(self.path_output)
    raster_name=os.path.splitext(base_raster_name)[0]
    self.iface.addRasterLayer(self.path_output, raster_name)


    layer=None
    for lyr in list(QgsProject.instance().mapLayers().values()):
        if lyr.name() == raster_name
            layer = lyr
        end
    end


    functions.applystyle(layer,'gr',0.5)


    tempoanalisi=time.time() - start_time
    tempostimato=time.strftime("%H:%M:%S", time.gmtime(tempoanalisi))
    messaggio="---------------------------------\n"
    messaggio+="Fine modellazione\n"
    messaggio+="\nTempo di analisi: "+tempostimato+"\n"
    messaggio+="---------------------------------\n\n"
    self.console.appendPlainText(messaggio)

    self.label_status.setText("In attesa di dati")
    self.label_status.setStyleSheet('color : green; font-weight:bold')
end