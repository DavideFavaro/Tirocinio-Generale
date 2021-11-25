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

class Dialog(EnviDialog):

    def __init__(self, iface):
        QDialog.__init__(self, iface.mainWindow())
        self.iface = iface
        self.canvas=self.iface.mapCanvas()
        #self.registry = QgsMapLayerRegistry.instance()
        self.msgBar = self.iface.messageBar()
        # Set up the user interface from Designer.
        self.setupUi(self)

        self.tabWidget.setCurrentIndex(0)

        self.tabWidget.removeTab(2)
        #pyqtRemoveInputHook()
        #pdb.set_trace()
        #self.reset_field_button.clicked.connect(self.reset_fields)
        #self.buttonBox.accepted.connect(self.leach)
        #self.button_run.clicked.connect(self.run1)
        self.label_title.setText("Analisi rumore")
        self.label_title.setStyleSheet('background-color : qlineargradient(spread:pad, x1:0, y1:0, x2:1, y2:0, stop:0 #4D4E53, stop:1 rgba(0, 0, 0, 0)); color : white')
        self.tableWidget.setRowCount(16)
        self.tableWidget.horizontalHeader().setStretchLastSection(True)

        self.combo_bound = QComboBox()
        self.combo_source = QComboBox()
        self.combo_lc = QComboBox()
        self.combofield_lc = QComboBox()
        self.combo_dem = QComboBox()
        self.combo_stability = QComboBox()
        self.combo_maindirwind = QComboBox()

        self.tableWidget.setCellWidget(0,0, self.combo_source)
        self.tableWidget.setCellWidget(1,0, self.combo_bound)
        self.tableWidget.setCellWidget(2,0, self.combo_lc)
        self.tableWidget.setCellWidget(3,0, self.combofield_lc)
        self.tableWidget.setCellWidget(4,0, self.combo_dem)
        self.tableWidget.setCellWidget(5,0, self.combo_stability)
        self.tableWidget.setCellWidget(6,0, self.combo_maindirwind)
        self.tableWidget.setItem(7 , 0, QTableWidgetItem(""))
        self.tableWidget.setItem(8 , 0, QTableWidgetItem(""))
        self.tableWidget.setItem(9 , 0, QTableWidgetItem(""))
        self.tableWidget.setItem(10 , 0, QTableWidgetItem(""))
        self.tableWidget.setItem(11 , 0, QTableWidgetItem(""))
        self.tableWidget.setItem(12 , 0, QTableWidgetItem(""))

        hbox = QHBoxLayout()
        hbox.addStretch(1)
        self.line_output = QLineEdit()
        self.line_output.setFixedHeight(25)
        self.saveButton = QPushButton("Scegli")
        self.saveButton.setFixedHeight(25)
        hbox.addWidget(self.line_output)
        hbox.addWidget(self.saveButton)
        cellWidget = QWidget()
        cellWidget.setLayout(hbox)

        self.tableWidget.setCellWidget(13,0, cellWidget)

        hbox1 = QHBoxLayout()
        hbox1.addStretch(1)
        self.line_folder = QLineEdit()
        self.line_folder.setFixedHeight(25)
        self.saveFolder = QPushButton("Scegli")
        self.saveFolder.setFixedHeight(25)
        hbox1.addWidget(self.line_folder)
        hbox1.addWidget(self.saveFolder)
        cellWidget1 = QWidget()
        cellWidget1.setLayout(hbox1)

        self.tableWidget.setCellWidget(14,0, cellWidget1)


        self.spinRes=QSpinBox()
        self.spinRes.setValue(25)

        self.tableWidget.setCellWidget(15,0, self.spinRes)


        self.tableWidget.setVerticalHeaderLabels((u'Vettoriale sorgente*', u'Vettoriale confine*',u'Vettoriale landcover*',
                                                  u'Campo landcover*',u'Mappa DTM*', u'Condizioni meteo',u'Direzione vento',u'Frequenza* (Khz)',
                                                  u'Velocità vento* (m/s)',u'Distanza misurazione (deafult 15 m)',u'Umidità* (%)',
                                                  u'Temperatura* (C°)',u'Livello sonoro sorgente* (dB)',u'Output file',u'Working folder',u'Risoluzione'))



        self.tableWidget.resizeRowsToContents();

        self.tabWidget.removeTab(1)

        self.menufile=self.menuFile_6

        self.spreadmenu=QAction("Spread",self)

        self.menufile.addAction(self.spreadmenu)

        self.saveButton.clicked.connect(lambda: self.scegli_file("salvaraster"))
        self.saveFolder.clicked.connect(lambda: self.scegli_file("folder"))
        self.reset_field_button.clicked.connect(self.reset_fields)
        self.actionCredits.triggered.connect(self.about)
        self.actionManuale.triggered.connect(self.help)
        self.spreadmenu.triggered.connect(self.spread_or)



        self.classiwind={}

        self.classiwind['N']=180
        self.classiwind['NW']=225
        self.classiwind['W']=270
        self.classiwind['SW']=315
        self.classiwind['S']=0
        self.classiwind['SE']=45
        self.classiwind['E']=90
        self.classiwind['NE']=135

        self.buttonBox.accepted.connect(self.run_spread)

        self.clear_out_button.clicked.connect(self.reset_output)
        self.save_out_button.clicked.connect(self.esporta_output)

        self.label_status.setText("In attesa di dati")
        self.label_status.setStyleSheet('color : green; font-weight:bold')

        self.popolacombo()

        self.combo_lc.currentIndexChanged[str].connect(self.checkfields)

    def run1(self):
        # fix_print_with_import
        print("run effettuato")


    def help(self):
        #self.credits = u"Università della Tuscia\n Viterbo - Italy\nRaffaele Pelorosso, Federica Gobattoni\nDeveloper: Francesco Geri"
        #QMessageBox.about(self.dlg,"Credits", self.credits )
        if platform.uname()[0]=="Windows":
            os.system("start "+os.path.dirname(__file__)+"/../tutorial/manuale_envifate_spreadgis.pdf")
        if platform.uname()[0]=="Linux":
            os.system("xdg-open "+os.path.dirname(__file__)+"/../tutorial/manuale_envifate_spreadgis.pdf")
        else:
            os.system("open "+os.path.dirname(__file__)+"/../tutorial/manuale_envifate_spreadgis.pdf")

    def spread_or(self):
        #self.credits = u"Università della Tuscia\n Viterbo - Italy\nRaffaele Pelorosso, Federica Gobattoni\nDeveloper: Francesco Geri"
        #QMessageBox.about(self.dlg,"Credits", self.credits )
        if platform.uname()[0]=="Windows":
            os.system("start "+os.path.dirname(__file__)+"/../tutorial/Harrisonetal1980.pdf")
        if platform.uname()[0]=="Linux":
            os.system("xdg-open "+os.path.dirname(__file__)+"/../tutorial/Harrisonetal1980.pdf")
        else:
            os.system("open "+os.path.dirname(__file__)+"/../tutorial/Harrisonetal1980.pdf")


    def checkfields(self):
        self.combofield_lc.clear()
        self.popolafields(self.combo_lc,self.combofield_lc)

    def popolafields(self,combo_in,combo_out):
        vect_source_text=combo_in.currentText()
        if vect_source_text!="":
            #vfields = self.allLayers[mainvect].pendingFields()
            mainvect = QgsProject.instance().mapLayersByName( vect_source_text )[0]
            vfields = mainvect.fields()
            # combo_out.addItem("No field")
            for field in vfields:
                combo_out.addItem(field.name())

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

    def reset_fields(self):


        self.console.clear()
        self.tableWidget.setItem(7 , 0, QTableWidgetItem(""))
        self.tableWidget.setItem(8 , 0, QTableWidgetItem(""))
        self.tableWidget.setItem(9 , 0, QTableWidgetItem(""))
        self.tableWidget.setItem(10 , 0, QTableWidgetItem(""))
        self.tableWidget.setItem(11 , 0, QTableWidgetItem(""))
        self.tableWidget.setItem(12 , 0, QTableWidgetItem(""))


        # for i in range(5,11):
        #     self.tableWidget.setItem(i , 0, QTableWidgetItem(""))

        # self.popolacombo()


    def esporta_output(self):
        try:
            resultmodel=self.console.toPlainText()
            self.f_res_export=open("risultatospread.txt","w")
            self.f_res_export.write(resultmodel.encode('utf8'))
            self.f_res_export.close()
            QMessageBox.information(self,"Info", "File risultatospread.txt exported in the working folder" )
        except:
            QMessageBox.information(self,"Info", "Nothing to export" )



    def reset_output(self):
        ret = QMessageBox.warning(self,"Attenzione", "Vuoi davvero eliminare i risultati del modello?",QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        if ret== QMessageBox.Yes:
            self.console.clear()
        else:
            return False

    def scegli_file(self,tipofile):
        if tipofile=="sqlite":
            self.fname = QFileDialog.getOpenFileName(None, 'Open file', '/home','sqlite3 files (*.sqlite);;all files (*.*)')
            self.dlg_conf.pathtodb.setText(self.fname)
        if tipofile=="csv":
            self.fname = QFileDialog.getOpenFileName(None, 'Open file', '/home','csv files (*.csv);;all files (*.*)')
            self.dlg_conf.path_to_kmean.setText(self.fname)
        if tipofile=="tif":
            self.fname = QFileDialog.getOpenFileName(None, 'Open file', '/home','GeoTiff files (*.tif);;all files (*.*)')
            self.dlg_reclass.output_raster_class.setText(self.fname)
        if tipofile=="tutti":
            self.fname = QFileDialog.getOpenFileName(None, 'Open file', '/home','all files (*.*)')
            self.dlg_reclass.input_reclass.setText(self.fname)
        if tipofile=="salvaraster":
            self.fname = QFileDialog.getSaveFileName(None, 'Save file', '/home','GeoTiff files (*.tif);;all files (*.*)')
            self.line_output.setText(self.fname)
        if tipofile=="folder":
            self.folder = QFileDialog.getExistingDirectory(self, "Select Directory")
            self.line_folder.setText(self.folder)
        # if self.tipofile=="salvacsv":
        #     self.fname = QFileDialog.getSaveFileName(None, 'Save file', '/home','csv files (*.csv);;all files (*.*)')
        #     self.dlg.lineEdit_csv.setText(self.fname)

    def about(self):
        QMessageBox.about(self, "Credits EnviFate",u"""<p>EnviFate: Open source tool for environmental risk analysis<br />Release 1.0<br />13-1-2017<br />License: GPL v. 3<br /><a href='https://bitbucket.org/fragit/envifate'>Home page plugin</a></p><hr><p>Lavoro svolto nell’ambito del  Progetto  di   ricerca   scientifica  “Definizione  di   metodi   standard     e  di strumenti applicativi   informatici per   il calcolo degli effetti dei fattori di perturbazione   ai sensi della decisione  2011/484/Ue,  da impiegarsi  nell’ambito  della valutazione di incidenza” finanziato dalla Regione Veneto. Partner principale è il DICAM, Dipartimento di Ingegneria Civile Ambientale e Meccanica dell’Università di Trento (Italia).</p><hr><p>Autori: Francesco Geri, Marco Ciolli</p><p>Universita' di Trento, Trento - Dipartimento di Ingegneria Civile Ambientale e Meccanica (DICAM) <a href="http://www.dicam.unitn.it/">www.dicam.unitn.it/</a></p><hr><p>Consulenti: Paolo Zatelli, Oscar Cainelli</p>""")


    def convert_seasonal_conditions(self,seas_cond):

        condiz=seas_cond[:2]

        meteophi={}

        meteophi['01']=180
        meteophi['02']=180
        meteophi['03']=0
        meteophi['04']=180
        meteophi['05']=144
        meteophi['06']=144
        meteophi['07']=62
        meteophi['08']=70
        meteophi['09']=90
        meteophi['10']=90

        return meteophi[condiz]


    def extract_values(self, raster,x,y):
        z=raster.dataProvider().identify(QgsPointXY(x, y),QgsRaster.IdentifyFormatValue)
        zresult=z.results()
        zvalue=zresult[1]
        return(zvalue)


    def spherical_spreading_loss(self, dist_ft,mis_dist):

        ssl_patch=dist_ft/mis_dist
        ssl=20*math.log(ssl_patch)

        return ssl


    def atmospheric_absorption_loss(self,elev_m, rh, temp_k, freq):
           # Calculate atmospheric absorption coefficient using ANSI S1.26-1995 standard

        # Convert elevation to atmospheric pressure
        p_a = 101.325 * (1 - (2.25577 * (10 ** (-5)) * elev_m)) ** 5.25588

        # Convert relative humidity to molar concentration of water vapor
        C = (-6.8346 * ((273.16 / temp_k) ** 1.261)) + 4.6151
        psat_pr = 10 ** C
        h = (rh) * (psat_pr) * ((p_a / 101.325) ** (-1))
        # Calculate derived values for subsequent equations
        pa_pr = p_a / 101.325
        T_Tr = temp_k / 293.15
        e = 2.7182818284

        # Calculate frO (equation 3)
        frO = ((pa_pr) * ((24 + (4.04 * 10000)) * h) * (0.02 + h)) / (0.391 + h)

        # Calculate frN (equation 4)
        frN = pa_pr * (T_Tr ** (-0.5)) * (9 + (280 * h * (e ** (-4.170 * ((T_Tr ** (-0.33333)) - 1)))))

        # Calculate alpha (equation 5)
        term1 = 1.84 * (10 ** (-11)) * (pa_pr ** (-1)) * (T_Tr ** 0.5)
        term2 = (T_Tr ** (-2.5)) * (0.01275 * (e ** (-2239.1 / temp_k)) * (frO / ((frO ** 2) + (freq ** 2))))
        term3 = 0.1068 * (e ** (-3352 / temp_k)) * (frN / ((frN ** 2) + (freq ** 2)))
        alpha = 8.686 * (freq ** 2)*(term1 + term2 + term3)

        return alpha


    def run_spread(self):
        #pyqtRemoveInputHook()


        # self.text_conc=str(self.tableWidget.item(6,0).text())
        # self.text_etemp=str(self.tableWidget.item(7,0).text())
        # self.text_height=str(self.tableWidget.item(8,0).text())
        # self.text_diameter=str(self.tableWidget.item(9,0).text())
        # self.text_temp=str(self.tableWidget.item(10,0).text())
        # self.text_gspeed=str(self.tableWidget.item(11,0).text())
        # self.text_wspeed=str(self.tableWidget.item(12,0).text())


        self.text_vector = str(self.combo_source.currentText())
        self.text_area = str(self.combo_bound.currentText())
        self.text_dem = str(self.combo_dem.currentText())
        self.text_lc = str(self.combo_lc.currentText())
        self.text_lcfield = str(self.combofield_lc.currentText())
        self.text_stability = str(self.combo_stability.currentText())
        self.text_line_dirwind = str(self.combo_maindirwind.currentText())


        self.text_line_dirwind=self.classiwind[self.combo_maindirwind.currentText()]

        self.text_freq=str(self.tableWidget.item(7,0).text())
        self.text_line_windspeed=str(self.tableWidget.item(8,0).text())
        self.text_line_misdist=str(self.tableWidget.item(9,0).text())
        self.text_line_umidita=str(self.tableWidget.item(10,0).text())
        self.text_line_temperatura=str(self.tableWidget.item(11,0).text())
        self.text_line_soundlevel=str(self.tableWidget.item(12,0).text())


        # self.class_stability=str(self.combo_stability.currentText()).split(" ")
        # self.c_stability=self.class_stability[1]

        # self.class_outdoor=str(self.combo_outdoor.currentText())
        # self.outdoor=self.classioutdoor[self.class_outdoor]




        self.source=self.listalayers[self.text_vector]

        if self.source.wkbType()!=1:
            QMessageBox.warning(self,"Warning", "Il vettoriale sorgente deve essere di geometria puntuale" )
            return

        self.areastudio=self.listalayers[self.text_area]


        if self.areastudio.wkbType()!=6:
            QMessageBox.warning(self,"Warning", "Il vettoriale are deve essere di geometria poligonale" )
            return



        self.lc=self.listalayers[self.text_lc]


        if self.lc.wkbType()!=6:
            QMessageBox.warning(self,"Warning", "Not a valid landcover geometry" )
            return


        self.dem=self.listalayers[self.text_dem]

        if not self.dem.isValid():
            QMessageBox.warning(self,"Warning", "The dem file is not valid" )
            return


        try:
            self.dirwind=180-int(self.text_line_dirwind)
            #self.lf=float(self.text_line_acquifer_depth)
            #self.dz=float(self.text_line_sourcethick)
        except Exception as e:
            QMessageBox.warning(self,"Warning", "Errore nella variabile direzione del vento" )
            return

        try:
            self.windspeed=int(self.text_line_windspeed)
            #self.lf=float(self.text_line_acquifer_depth)
            #self.dz=float(self.text_line_sourcethick)
        except Exception as e:
            QMessageBox.warning(self,"Warning", u"Errore nella variabile velocità del vento" )
            return


        try:
            self.umidita=int(self.text_line_umidita)
            #self.lf=float(self.text_line_acquifer_depth)
            #self.dz=float(self.text_line_sourcethick)
        except Exception as e:
            QMessageBox.warning(self,"Warning", u"Errore nel dato di umidità (%)" )
            return


        try:
            self.temperatura=int(self.text_line_temperatura)
            #self.lf=float(self.text_line_acquifer_depth)
            #self.dz=float(self.text_line_sourcethick)
        except Exception as e:
            QMessageBox.warning(self,"Warning", u"Errore nel dato di temperatura (intero C°)" )
            return


        try:
            self.soundlevel=int(self.text_line_soundlevel)
            #self.lf=float(self.text_line_acquifer_depth)
            #self.dz=float(self.text_line_sourcethick)
        except Exception as e:
            QMessageBox.warning(self,"Warning", u"Errore nel dato del livello sonoro in sorgente (intero dB)" )
            return

        try:
            self.freq=int(self.text_freq)
        except Exception as e:
            QMessageBox.warning(self,"Warning", "Errore nella frequenza" )
            return


        if self.text_line_misdist!='':
            try:
                self.misdist=int(self.text_line_misdist)
            except Exception as e:
                QMessageBox.warning(self,"Warning", u"Errore nella distanza di misurazione (intero in metri)" )
                return
        else:
            self.misdist=15


        self.path_output=self.line_output.text()
        if self.path_output=="":
            self.path_output=os.path.dirname(__file__)+"/output_model.tif"


        self.path_working=self.line_folder.text()
        if self.path_working=="":
            self.path_working=os.path.dirname(__file__)

        self.path_temp_lc=self.path_working+"/temp_lc.tif"



        self.res=int(self.spinRes.text())

        # self.x_w=self.classiwind[self.combo_maindirwind.currentText()]

        messaggio="Inizio elaborazione SpreadGIS Envifate\n"
        messaggio+="---------------------------\n\n"
        messaggio+="FILE DI INPUT:\n"
        messaggio+="Vettoriale sorgente: "+str(self.text_vector)+"\n"
        messaggio+="Vettoriale confine: "+str(self.text_area)+"\n"
        messaggio+="Vettoriale uso del suolo: "+str(self.text_lc)+"\n"
        messaggio+="DTM: "+str(self.text_dem)+"\n\n"
        messaggio+="VARIABILI:\n"
        messaggio+="Condizioni meteo: "+str(self.text_stability)+"\n"
        messaggio+="Direzione vento: "+str(self.text_line_dirwind)+"\n"
        messaggio+=u"Velocità vento: "+str(self.text_line_windspeed)+"\n"
        messaggio+="Distanza di misurazione: "+str(self.misdist)+"\n"
        messaggio+=u"Umidità: "+str(self.text_line_umidita)+"\n"
        messaggio+="Temperatura: "+str(self.text_line_temperatura)+"\n"
        messaggio+="Livello sonoro sorgente: "+str(self.text_line_soundlevel)+"\n"
        messaggio+="Risoluzione: "+str(self.res)+"\n\n"
        messaggio+='ALGORITMO UTILIZZATO: Spread (Harrison, Robin T., Roger N. Clark, and George H. Stankey. "Predicting impact of noise on recreationists." Predicting impact of noise on recreationists. (1980).)\n\n'
        messaggio+="---------------------------\n\n"
        self.console.appendPlainText(messaggio)


        self.label_status.setText("Preparazione dati")
        self.label_status.setStyleSheet('color : #e8b445;font-weight:bold')


        phi=self.convert_seasonal_conditions(self.text_stability)

 
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


        pixel_size = self.res
        NoData_value = -9999


        # Create the destination data source
        x_res = int((x_max - x_min) / pixel_size)
        y_res = int((y_max - y_min) / pixel_size)
        #target_ds = drivermem.Create('', x_res, y_res, 1, gdal.GDT_Byte)
        target_ds = gdal.GetDriverByName('GTiff').Create(self.path_output, int(x_res), int(y_res), 1, gdal.GDT_Float32)
        target_ds.SetGeoTransform((x_min, pixel_size, 0, y_max, 0, -pixel_size))
        projectionfrom = target_ds.GetProjection()
        band = target_ds.GetRasterBand(1)
        band.SetNoDataValue(float(NoData_value))
        band.Fill(NoData_value)
        xsize = band.XSize
        ysize = band.YSize
        outData = np.array(band.ReadAsArray(0, 0, xsize,ysize).astype(np.float))





        #path__layer_lc=lc_clip['OUTPUT'].dataProvider().dataSourceUri()
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