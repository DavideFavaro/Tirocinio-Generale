module NoiseSingleSource

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


import functions, noise



function transmission_loss( r::Real )
    return 20log10(r)
end



function atmospheric_absorpion_loss( r::Real, height_m::Real, relative_humidity::Real, temperature_k::Real, frequency::Real )
    # Calculate atmospheric absorption coefficient using ANSI S1.26-1995 standard
    # Convert elevation to atmospheric pressure
    atmospheric_pressure = 101.325( 1 - ( 2.25577 * 10^(-5) * height_m ) )^5.25588
    # Calculate derived values for subsequent equations
    p_atm_pressure = atmospheric_pressure / 101.325
    t_tr = temperature_k / 293.15
    # Convert relative humidity to molar concentration of water vapor
    C = ( -6.8346( 273.16 / temperature_k )^1.261 ) + 4.6151
    p_saturation_pressure = 10^C
    humidity = relative_humidity * p_saturation_pressure * p_atm_pressure^(-1)
    # Calculate relaxation frequency of O (equation 3)
    #   frO₂ = ( p_atm_pressure * ( (24 + 4.04e04) * humidity ) * (0.02 + humidity) ) / (0.391 + humidity)
    frO₂ = p_atm_pressure * ( 24 + ( 4.04e04humidity * ( (0.02 + humidity) / (0.391 + humidity) ) ) )
    # Calculate relaxation frequency of N (equation 4)
    frN₂ = p_atm_pressure * √t_tr * ( 9 + ( 280humidity * ℯ^( -4.170(t_tr^(-1/3) - 1) ) ) )
    # Calculate alpha (equation 5)
    term1 = 1.84 * 10^(-11) * p_atm_pressure^(-1) * √t_tr
    #   term2 = t_tr^(-2.5) * ( 0.01275 * ℯ^(-2239.1 / temp_k) * ( frO₂ / (frO₂^2 + freq^2) ) )
    term2 = t_tr^(-2.5)( 0.01275ℯ^(-2239.1 / temp_k) / ( frO₂ + (freq^2 / frO₂) ) )
    #   term3 = 0.1068 * ℯ^(-3352 / temp_k) * ( frN₂ / (frN₂^2 + freq^2) )
    term3 = t_tr^(-2.5) * ( 0.1068ℯ^(-3352 / temp_k) / ( frN₂ + (freq^2 / frN₂) ) ) 
    #   α = 8.686 * (frequency^2) * ( term1 + term2 + term3 )
    α = frequency^2 * (term1 + term2 + term3)

    return α * r / 100
end


function convert_seasonal_conditions(season_conditions)
    conditions = season_conditions[:2]
    meteophi = Dict(
        "01" => 180,
        "02" => 180,
        "03" => 0,
        "04" => 180,
        "05" => 144,
        "06" => 144,
        "07" => 62,
        "08" => 70,
        "09" => 90,
        "10" => 90
    )
    return meteophi[conditions]
end


#=
    def help(self):
        #self.credits = u"Università della Tuscia\n Viterbo - Italy\nRaffaele Pelorosso, Federica Gobattoni\nDeveloper: Francesco Geri"
        #QMessageBox.about(self.dlg,"Credits", self.credits )
        if platform.uname()[0]=="Windows":
            os.system("start "+os.path.dirname(__file__)+"/../tutorial/manuale_envifate_spreadgis.pdf")
        if platform.uname()[0]=="Linux":
            os.system("xdg-open "+os.path.dirname(__file__)+"/../tutorial/manuale_envifate_spreadgis.pdf")
        else:
            os.system("open "+os.path.dirname(__file__)+"/../tutorial/manuale_envifate_spreadgis.pdf")

    def popolacombo(self):
        self.progressBar.setValue(0)
        self.combo_bound.clear()
        self.combo_source.clear()
        self.combo_lc.clear()
        self.combofield_lc.clear()
        self.combo_dem.clear()
        self.combo_stability.clear()
        self.combo_maindirwind.clear()


        self.combo_stability.addItem("01 Estate, sereno, vento debole o assente, giorno")
        self.combo_stability.addItem("02 Inverno, sereno, vento debole o assente, giorno")
        self.combo_stability.addItem("03 Estate, sereno, vento debole o assente, notte")
        self.combo_stability.addItem("04 Inverno, sereno, vento debole o assente, notte")
        self.combo_stability.addItem("05 Estate, sereno, vento presente, giorno")
        self.combo_stability.addItem("06 Inverno, sereno, vento presente, giorno")
        self.combo_stability.addItem("07 Estate, sereno, vento presente, notte")
        self.combo_stability.addItem("08 Inverno, sereno, vento presente, notte")
        self.combo_stability.addItem("09 Cielo nuvoloso, vento debole o assente")
        self.combo_stability.addItem("10 Cielo nuvoloso, vento presente")



        for (key, nome) in list(self.classiwind.items()):
            self.combo_maindirwind.addItem(key)


        self.allLayers = self.canvas.layers()
        self.listalayers=dict()
        #elementovuoto="No required"
        for i in self.allLayers:
            if i.type() == QgsMapLayer.VectorLayer:
                self.listalayers[i.name()]=i
                self.combo_source.addItem(str(i.name()))
                self.combo_bound.addItem(str(i.name()))
                self.combo_lc.addItem(str(i.name()))
            if i.type()==QgsMapLayer.RasterLayer:
                self.listalayers[i.name()]=i
                self.combo_dem.addItem(str(i.name()))

        self.popolafields(self.combo_lc,self.combofield_lc)

        conn = sqlite3.connect(os.path.dirname(__file__)+"/../library/substance.db")
        cursor=conn.cursor()
        query_substance="select id,nome from substance"
        cursor.execute(query_substance)
        sql_fetch=cursor.fetchall()

        self.inquinanti=dict()
        # for row in sql_fetch:
        #     self.inquinanti[row[1]]=row[0]
        #     self.combo_contaminant.addItem(row[1])

        query_texture="select id,nome from texture"
        cursor.execute(query_texture)
        sql_texture=cursor.fetchall()

        self.texture=dict()
        # for rowt in sql_texture:
        #     self.texture[rowt[1]]=rowt[0]
        #     self.combo_texture.addItem(rowt[1])


        conn.close()
    
    def extract_values(self, raster,x,y):
        z=raster.dataProvider().identify(QgsPointXY(x, y),QgsRaster.IdentifyFormatValue)
        zresult=z.results()
        zvalue=zresult[1]
        return(zvalue)
=#
                   #                                  dirwind               windspeed        misdist                      umidita         temperatura        freq
function run_spread(dem, source, landcover, weather, wind_direction::Real, wind_speed::Real, measure_distance::Real=15.0, humidity::Real, temperature::Real, frequency::Real,
                    soundlevel::Real, resolution::Integer, folder::AbstractString=".\\"  )

 """ FORSE MANCENO QUESTI PARAMETRI
    self.text_lcfield = str(self.combofield_lc.currentText())
    self.text_freq=str(self.tableWidget.item(7,0).text())
 """

    if agd.geomdim(source) != 0
        throw(DomainError(source, "`source` must be a point"))
    end

    landcover_layer = agd.getlayer(landcover, 0) 

    if agd.geomdim(landcover_layer) != 2
        throw(DomainError(source, "Not a valid `landcover` geometry"))
    end

    # wind_direction = 180 - wind_direction

    output_pat = folder * "output_model.tiff"
    temp_landcover_path = folder * "\\temp_landcover.tiff"

    # messaggio+='ALGORITMO UTILIZZATO: Spread (Harrison, Robin T., Roger N. Clark, and George H. Stankey. "Predicting impact of noise on recreationists." Predicting impact of noise on recreationists. (1980).)\n\n'

    ϕ = convert_seasonal_conditions(weather)





 """ """
    lc_clip_proc = processing.run('qgis:clip', {'INPUT':self.lc, 'OVERLAY':self.areastudio, 'OUTPUT':self.path_working+'/clip.gpkg'})


    # pyqtRemoveInputHook()
    # pdb.set_trace()

    lc_clip=QgsVectorLayer(lc_clip_proc['OUTPUT'], 'lc_clip', 'ogr')

    path_layer=self.areastudio.dataProvider().dataSourceUri()
    path=path_layer.split("|")
    source_ds = ogr.Open(path[0])
    area_layer = source_ds.GetLayer()
    #x_min, x_max, y_min, y_max = area_layer.GetExtent()
    x_min=int(area_layer.GetExtent()[0])
    y_min=int(area_layer.GetExtent()[2])
    x_max=int(area_layer.GetExtent()[1])
    y_max=int(area_layer.GetExtent()[3])
""" """



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






    path__layer_lc=lc_clip.dataProvider().dataSourceUri()
    path_lc=path__layer_lc.split("|")
    source_ds_lc = ogr.Open(path_lc[0])
    lc_layer = source_ds_lc.GetLayer()
    lc_ds = gdal.GetDriverByName('GTiff').Create(self.path_temp_lc, x_res, y_res, 1, gdal.GDT_Float32)
    lc_ds.SetGeoTransform((x_min, 25, 0, y_max, 0, -25))
    lc_ds.SetProjection(projectionfrom)
    band_lc = lc_ds.GetRasterBand(1)
    band_lc.SetNoDataValue(float(-9999))
    band_lc.Fill(-9999)
    xsize = band_lc.XSize
    ysize = band_lc.YSize

    #pdb.set_trace()

    gdal.RasterizeLayer(lc_ds, [1], lc_layer,options=["ATTRIBUTE="+self.text_lcfield])
    lc_ds=None

    lc_layer=QgsRasterLayer(self.path_temp_lc,"lc_layer")



    ###### inizio file di controllo

    # eucdist = ssl

    self.path_eucdist=self.path_working+"/step1.tif"

    target_ds_eucdist = gdal.GetDriverByName('GTiff').Create(self.path_eucdist, x_res, y_res, 1, gdal.GDT_Float32)
    target_ds_eucdist.SetGeoTransform((x_min, pixel_size, 0, y_max, 0, -pixel_size))
    projectionfrom = target_ds_eucdist.GetProjection()
    band_eucdist = target_ds_eucdist.GetRasterBand(1)
    band_eucdist.SetNoDataValue(float(NoData_value))
    band_eucdist.Fill(NoData_value)
    xsize = band_eucdist.XSize
    ysize = band_eucdist.YSize
    outData_eucdist = np.array(band_eucdist.ReadAsArray(0, 0, xsize,ysize).astype(np.float))



    #aal
    self.path_aal=self.path_working+"/step2.tif"
    target_ds_aal = gdal.GetDriverByName('GTiff').Create(self.path_aal, x_res, y_res, 1, gdal.GDT_Float32)
    target_ds_aal.SetGeoTransform((x_min, pixel_size, 0, y_max, 0, -pixel_size))
    projectionfrom = target_ds_aal.GetProjection()
    band_aal = target_ds_aal.GetRasterBand(1)
    band_aal.SetNoDataValue(float(NoData_value))
    band_aal.Fill(NoData_value)
    xsize = band_aal.XSize
    ysize = band_aal.YSize
    outData_aal = np.array(band_aal.ReadAsArray(0, 0, xsize,ysize).astype(np.float))


    #max_veg_loss
    self.path_mvl=self.path_working+"/step3.tif"
    target_ds_mvl = gdal.GetDriverByName('GTiff').Create(self.path_mvl, x_res, y_res, 1, gdal.GDT_Float32)
    target_ds_mvl.SetGeoTransform((x_min, pixel_size, 0, y_max, 0, -pixel_size))
    projectionfrom = target_ds_mvl.GetProjection()
    band_mvl = target_ds_mvl.GetRasterBand(1)
    band_mvl.SetNoDataValue(float(NoData_value))
    band_mvl.Fill(NoData_value)
    xsize = band_aal.XSize
    ysize = band_aal.YSize
    outData_mvl = np.array(band_mvl.ReadAsArray(0, 0, xsize,ysize).astype(np.float))


    #bar
    self.path_bar=self.path_working+"/step4.tif"
    target_ds_bar = gdal.GetDriverByName('GTiff').Create(self.path_bar, x_res, y_res, 1, gdal.GDT_Float32)
    target_ds_bar.SetGeoTransform((x_min, pixel_size, 0, y_max, 0, -pixel_size))
    projectionfrom = target_ds_bar.GetProjection()
    band_bar = target_ds_bar.GetRasterBand(1)
    band_bar.SetNoDataValue(float(NoData_value))
    band_bar.Fill(NoData_value)
    xsize = band_aal.XSize
    ysize = band_aal.YSize
    outData_bar = np.array(band_bar.ReadAsArray(0, 0, xsize,ysize).astype(np.float))


    #windloss
    self.path_wind=self.path_working+"/temp_wind.tif"
    target_ds_wind = gdal.GetDriverByName('GTiff').Create(self.path_wind, x_res, y_res, 1, gdal.GDT_Float32)
    target_ds_wind.SetGeoTransform((x_min, pixel_size, 0, y_max, 0, -pixel_size))
    projectionfrom = target_ds_wind.GetProjection()
    band_wind = target_ds_wind.GetRasterBand(1)
    band_wind.SetNoDataValue(float(NoData_value))
    band_wind.Fill(NoData_value)
    xsize = band_aal.XSize
    ysize = band_aal.YSize
    outData_wind = np.array(band_wind.ReadAsArray(0, 0, xsize,ysize).astype(np.float))




    ###### fine file temporanei

    #array_area = np.array(band.ReadAsArray())
    feature = next(self.source.getFeatures())
    geom = feature.geometry().asPoint()
    x_source=geom[0]
    y_source=geom[1]


    rows=ysize-1
    cols=xsize-1


    npts = 100

    m2ft = 3.28084

    misdist_ft=self.misdist*m2ft

    temp_k = float(self.temperatura) + 273.15

    rh = float(self.umidita)
    freq_f = float(self.freq)

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

                # calcolo distanza (metri e piedi)
                deltax=x-x_source
                deltay=y-y_source
                dist=math.sqrt(math.pow(deltay,2)+math.pow(deltax,2))
                dist_ft=dist*m2ft #distanza espressa in piedi (attenzione equivale a distance in spreadgishlpr.py)

                #inizio analisi rumore

                dist_vec = noise.find_distances(dist, npts) # vettore distanza

                elev_source=self.extract_values(self.dem,x_source,y_source)
                elev_target=self.extract_values(self.dem,x,y)

                t_slope=(elev_target-elev_source)/dist #calcolo della linea di pendenza tra sorgente e punto

                x_coords = noise.find_coords(x_source, x, npts)
                y_coords = noise.find_coords(y_source, y, npts)
                hgt_vec=[]
                veg_cut=[]
                max_height = 0


                # pdb.set_trace()

                bar_dist=dist_ft

                #spherical spreading loss
                ssl=self.spherical_spreading_loss(dist_ft,misdist_ft)

                #atmospheric adsorption loss

                aal_loss=self.atmospheric_absorption_loss(elev_target,rh,temp_k,freq_f)

                for i in range(len(dist_vec)):

                    hgt_vec.append(self.extract_values(self.dem,x_coords[i],y_coords[i])) # vettore elevazione
                    veg_cut.append(self.extract_values(lc_layer,x_coords[i],y_coords[i])) # vettore corine

                    #prendi il massimo ostacolo nella traiettoria sorgente punto e memorizza distanza ostacolo e altezza ostacolo
                    slope_elevation=t_slope*dist_vec[i]*m2ft+elev_source
                    height_above_slope = dist_vec[i] - slope_elevation
                    if height_above_slope > max_height:
                        max_height = height_above_slope
                        bar_dist = dist_vec[i] * m2ft



                term1 = (max_height**2 + bar_dist**2)**(0.5)

                term2 = (max_height**2 + (dist_ft - bar_dist)**2)**(0.5)

                BPD_1 = term1 + term2 - dist_ft


                BPD = BPD_1 if BPD_1>0 else 0

                # pyqtRemoveInputHook()
                # pdb.set_trace()

                con = list([x for x in veg_cut if x==312])
                hwd = list([x for x in veg_cut if x==311 or x==313])
                heb = list([x for x in veg_cut if x==321 or x==322 or x==323 or x==324 ])


                if len(con) == 0:
                    max_con_loss = 0
                else:
                    distance_con = (float(len(con)) / npts) * dist
                    max_con_loss = 5.2504 * math.log(distance_con) - 9.8094 # R2 = 0.99
                    if max_con_loss < 0:
                        max_con_loss = 0

                if len(hwd) == 0:
                    max_hwd_loss = 0
                else:
                    dist_hwd = (float(len(hwd)) / npts) * dist
                    max_hwd_loss = 6.6224 * math.log(dist_hwd) - 16.762 # R2 = 0.99
                    if max_hwd_loss < 0:
                        max_hwd_loss = 0

                max_heb_loss = 0
                if len(heb) > 0:
                    max_heb_loss = 4


                # Add sources of vegetation loss
                max_veg_loss = max_con_loss + max_hwd_loss + max_heb_loss

                # Cap total vegetation loss at 14 dB
                if max_veg_loss > 14:
                    max_veg_loss = 14


                L = ((0.0000000000005*self.freq**4) - (0.000000001*self.freq**3) - (0.0000004*self.freq**2) + (0.0028*self.freq) - (0.3051))

                bar_factor=L*BPD

                bar=13.573* (bar_factor**0.2299)




                #inizio analisi gradi
                m_gradi=0


                try:
                    m=(y-y_source)/(x-x_source)

                    if x>x_source:
                        m_gradi=(math.atan(m)*180)/math.pi+270
                    if x<x_source:
                        m_gradi=(math.atan(m)*180)/math.pi+90
                except:
                    m_gradi=0
                m_gradi2=m_gradi+self.dirwind
                if m_gradi2>360:
                    m_gradi3=m_gradi2-360
                else:
                    if m_gradi2<0:
                        m_gradi3=m_gradi2+360
                    else:
                        m_gradi3=m_gradi2

                trueValue = 360 - m_gradi3
                if m_gradi3>180:
                    gradireali=trueValue
                else:
                    gradireali=m_gradi3

                updownwind=phi-gradireali

                if updownwind>0:
                    if updownwind>=50:
                        wind_loss=25
                    else:
                        wind_loss=5.7642 * math.log(updownwind) + 2.5664
                elif updownwind<=0:
                    freq_dist = dist_ft*self.freq
                    if freq_dist <= 406237:
                        wind_loss=0
                    else:
                        wind_loss=4.2598 * math.log(freq_dist) - 55.014
                else:
                    wind_loss=0





                bar_wind_loss=0
                if (bar+wind_loss)> 25:
                    bar_wind_loss=25
                else:
                    bar_wind_loss=bar+wind_loss

                #fine analisi gradi

                #bar_veg_wind=max_veg_loss+bar_wind_loss #effetto comulativo vegetazione e terreno



                if (self.soundlevel-ssl)>0:
                    ssl_loss=self.soundlevel-ssl
                else:
                    ssl_loss=0

                if (ssl_loss -  aal_loss)>0:
                    sslaal = ssl_loss -  aal_loss
                else:
                    sslaal=0

                if (sslaal-bar_wind_loss)>0:
                    sslaal1=sslaal-bar_wind_loss
                else:
                    sslaal1=0

                total_loss=sslaal1-max_veg_loss


                if total_loss>0:
                    outData[row,col]=total_loss
                else:
                    outData[row,col]=0


                ##### inizio scrittura file temporanei
                outData_eucdist[row,col]=self.soundlevel
                outData_aal[row,col]=ssl_loss
                outData_mvl[row,col]=sslaal
                outData_bar[row,col]=sslaal1
                outData_wind[row,col]=wind_loss
                ##### fine scrittura file temporanei

        self.label_status.setText("Preparazione output")
        self.label_status.setStyleSheet('color : #e8b445;font-weight:bold')

        outData_raster=outData[::-1]
        band.WriteArray(outData_raster)

        ##### inizio rotazione file temporanei

        outData_eucdist_raster=outData_eucdist[::-1]
        band_eucdist.WriteArray(outData_eucdist_raster)

        outData_aal_raster=outData_aal[::-1]
        band_aal.WriteArray(outData_aal_raster)

        outData_mvl_raster=outData_mvl[::-1]
        band_mvl.WriteArray(outData_mvl_raster)

        outData_bar_raster=outData_bar[::-1]
        band_bar.WriteArray(outData_bar_raster)

        outData_wind_raster=outData_wind[::-1]
        band_wind.WriteArray(outData_wind_raster)

        ##### fine rotazione file temporanei


    band= None
    target_ds = None


    ####### inizio reset file temporanei
    band_eucdist= None
    target_ds_eucdist= None


    band_aal= None
    target_ds_aal= None

    band_mvl= None
    target_ds_mvl= None

    band_bar= None
    target_ds_bar= None

    band_wind= None
    target_ds_wind= None


    ####### fine reset file temporanei



    base_raster_name=os.path.basename(self.path_output)
    raster_name=os.path.splitext(base_raster_name)[0]
    self.outputlayer=self.iface.addRasterLayer(self.path_output, raster_name)

    layer=None
    for lyr in list(QgsProject.instance().mapLayers().values()):
        if lyr.name() == raster_name:
            layer = lyr


    renderer = layer.renderer()
    transparency = renderer.rasterTransparency()
    ltr = QgsRasterTransparency.TransparentSingleValuePixel()


    tr_list = []
    ltr.min = 0
    ltr.max = 0
    ltr.percentTransparent = 100

    tr_list.append(ltr)

    transparency.setTransparentSingleValuePixelList(tr_list)
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


end # module