module WaterNoise

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
            os.system("start "+os.path.dirname(__file__)+"/../tutorial/manuale_envifate_rumore_in_acqua.pdf")
        if platform.uname()[0]=="Linux":
            os.system("xdg-open "+os.path.dirname(__file__)+"/../tutorial/manuale_envifate_rumore_in_acqua.pdf")
        else:
            os.system("open "+os.path.dirname(__file__)+"/../tutorial/manuale_envifate_rumore_in_acqua.pdf")
=#


function run_f1(S,T)
    return √(0.78(S / 35)) * ℯ^(T/26)
end



function run_f2(T)
    return 42ℯ^(T/17)
end



function run_spread( source, area, depth::Real, salinity::Real, pH::Real, temperature::Real, frequency::Real, resolution::Integer, output_path::AbstractString=".\\" )

    if pH < 0 || pH > 14
        throw(DomainError(pH, "`pH` must be a number between 0 and 14"))
    end

    if agd.geomdim(source) != 0
        throw(DomainError(source, "`source` must be a point"))
    end

    if agd.geomdim(area) != 2
        throw(DomainError(source, "`area` must be a polygon"))
    end

    refsys = agd.getspatialref(source)
    layers = agd.getlayer(area, 0)

    if agd.getspatialref(layers) != refsys
        throw(DomainError("The reference systems are not uniform. Aborting analysis." ))
    end

    measure_dist = 15.0

    path_temp_lc=self.path_working+"/temp_lc.tif"


    listafrequenze=[]



    
    if self.check_freq_list==True:
        with open(self.line_freqList.text(), newline='') as csvfile:
            spamreader = csv.reader(csvfile, delimiter=',')
            for row in spamreader:
                listafrequenze.append(row[0])
    else:
        listafrequenze.append(str(self.freq))

    for frequenza in listafrequenze:


        self.path_output=self.path_working+"/"+self.line_output+frequenza+".tif"
        if self.line_output=="":
            #self.path_output=os.path.dirname(__file__)+"/sound_level"+frequenza+".tif"
            self.path_output=self.path_working+"/sound_level"+frequenza+".tif"

        messaggio="Inizio elaborazione Analisi del rumore in acqua\n"
        messaggio+="---------------------------\n\n"
        messaggio+="FILE DI INPUT:\n"
        messaggio+="Vettoriale sorgente: "+str(self.text_vector)+"\n"
        messaggio+="Vettoriale confine: "+str(self.text_area)+"\n"
        messaggio+="VARIABILI:\n"
        messaggio+="Salinità: "+str(self.line_salinity)+"\n"
        messaggio+="Profondità: "+str(self.line_depth)+"\n"
        messaggio+=u"Ph: "+str(self.line_ph)+"\n"
        messaggio+="Temperatura: "+str(self.line_temperature)+"\n"
        messaggio+="Risoluzione: "+str(self.res)+"\n\n"
        messaggio+='ALGORITMO UTILIZZATO: Ainslie, M. A., & McColm, J. G. (1998). A simplified formula for viscous and chemical absorption in sea water. The Journal of the Acoustical Society of America, 103(3), 1671-1672.)\n\n'
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


        pixel_size = self.res
        NoData_value = -9999


        # Create the destination data source
        x_res = int((x_max - x_min) / pixel_size)
        y_res = int((y_max - y_min) / pixel_size)
        #target_ds = drivermem.Create('', x_res, y_res, 1, gdal.GDT_Byte)
        target_ds = gdal.GetDriverByName('GTiff').Create(self.path_output, int(x_res), int(y_res), 1, gdal.GDT_Float32)
        target_ds.SetGeoTransform((x_min, pixel_size, 0, y_max, 0, -pixel_size))
        # if self.srid!="":
        #     srs = osr.SpatialReference()
        #     srs.ImportFromEPSG(self.srid)
        #     target_ds.SetProjection( srs.ExportToWkt() )
        # else:
        #     target_ds.SetProjection(projectionfrom)
        # self.refsys
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(int(self.refsys))
        target_ds.SetProjection( srs.ExportToWkt() )

        target_ds.SetMetadata({'credits':'Envifate - Francesco Geri, Oscar Cainelli, Paolo Zatelli, Gianluca Salogni, Marco Ciolli - DICAM Università degli Studi di Trento - Regione Veneto',
                               'modulo':'Dispersione rumore in acqua',
                               'descrizione':'Analisi della dispersione acustica in acqua',
                               'srs':self.source.crs().authid(),
                               'data':datetime.datetime.now().strftime("%d-%m-%y")})

        band = target_ds.GetRasterBand(1)
        band.SetNoDataValue(float(NoData_value))
        band.Fill(NoData_value)
        xsize = band.XSize
        ysize = band.YSize
        outData = np.array(band.ReadAsArray(0, 0, xsize,ysize).astype(np.float))


        nfeature=0

        features=self.source.getFeatures()

        for feature in features:
            #feature = next(self.source.getFeatures())
            geom = feature.geometry().asPoint()
            x_source=geom[0]
            y_source=geom[1]


            idxlevel = self.source.fields().indexFromName('level')
            soundlevel=feature.attributes()[idxlevel]

            nfeature+=1


            rows=ysize-1
            cols=xsize-1


            npts = 100


            freq_f = float(frequenza)

            max_progress=rows*cols
            self.progressBar.setMaximum(max_progress)
            start_time = time.time()

            self.label_status.setText("Processing data for point n°"+str(nfeature))
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

                        # calcolo distanza
                        deltax=x-x_source
                        deltay=y-y_source
                        dist=math.sqrt(math.pow(deltay,2)+math.pow(deltax,2))


                        f1=self.run_f1(self.salinity,self.temperature)
                        f2=self.run_f2(self.temperature)

                        alfa1=0.106*((f1*(freq_f**2))/((freq_f**2)+(f1**2)))*math.exp((self.ph-8)/0.56)
                        alfa2=0.52*(1+self.temperature/43)*(self.salinity/35)*((f2*(freq_f**2))/((freq_f**2)+(f1**2)))*math.exp((-self.depth)/6)
                        alfa3=0.00049*(freq_f**2)*math.exp(-((self.temperature/27)+(self.depth/17)))

                        alfa=alfa1+alfa2+alfa3

                        tl=(20*math.log(dist))+alfa

                        total_loss=soundlevel-tl

                        # pyqtRemoveInputHook()
                        # pdb.set_trace()
                        if nfeature==1:
                            if total_loss>0:
                                outData[row,col]=total_loss
                            else:
                                outData[row,col]=0
                        else:

                            if total_loss>0:
                                outData[row,col]=outData[row,col]+total_loss
                            else:
                                outData[row,col]=outData[row,col]

                        ##### inizio scrittura file temporanei
                        # outData_eucdist[row,col]=self.soundlevel
                        # outData_aal[row,col]=ssl_loss
                        # outData_mvl[row,col]=sslaal
                        # outData_bar[row,col]=sslaal1
                        # outData_wind[row,col]=wind_loss
                        ##### fine scrittura file temporanei

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


        functions.applystyle(layer,'viridis',0.5)

        # pyqtRemoveInputHook()
        # pdb.set_trace()


        layer.triggerRepaint()

        tempoanalisi=time.time() - start_time
        tempostimato=time.strftime("%H:%M:%S", time.gmtime(tempoanalisi))
        messaggio="---------------------------------\n"
        messaggio+="Fine modellazione\n"
        messaggio+="\nTempo di analisi: "+tempostimato+"\n"
        messaggio+="---------------------------------\n\n"
        self.console.appendPlainText(messaggio)

        self.label_status.setText("In attesa di dati")
        self.label_status.setStyleSheet('color : green; font-weight:bold')
