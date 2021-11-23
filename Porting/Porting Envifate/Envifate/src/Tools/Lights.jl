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
    from __future__ import print_function

    from builtins import str
    from builtins import range
    from qgis.PyQt import QtCore, QtGui
    from PyQt5.QtCore import QSettings, QTranslator, QCoreApplication, Qt, QObject, pyqtSignal, pyqtRemoveInputHook
    from PyQt5.QtGui import QIcon
    from PyQt5.QtWidgets import QAction, QDialog, QFormLayout, QMenu, QComboBox, QTableWidgetItem, QHBoxLayout, QLineEdit, QPushButton, QWidget, QSpinBox, QTableWidgetItem, QMessageBox, QFileDialog


    import sys

    # Initialize Qt resources from file resources.py
    #import resources


    # Import the code for the dialog

    #from open_risk_dialog import OpenRiskDialog
    import os.path
    try:
      import sqlite3
    except:
      # fix_print_with_import
      print("librerie per la connessione al database sqlite non trovate")

    import qgis
    from qgis.core import *
    from qgis.gui import *
    from qgis.utils import iface

    import numpy as np
    import math
    import time
    from osgeo import gdal,ogr,osr
    import platform

    import processing
    from processing.core.Processing import Processing
    Processing.initialize()

    import datetime

    import pdb

    from envifate_dialog import EnviDialog

    from configuration_dialog import ConfigurationDialog

    sys.path.append( os.path.dirname(__file__)+"/../library" )

    import functions, lake, do_setting
=#
#= 
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

        self.tab_2.setEnabled(False)

        self.label_title.setText("Analisi inquinamento luminoso")
        self.label_title.setStyleSheet('background-color : qlineargradient(spread:pad, x1:0, y1:0, x2:1, y2:0, stop:0 #e5e94c, stop:1 rgba(0, 0, 0, 0)); color : black')
        self.tableWidget.setRowCount(9)
        self.tableWidget.horizontalHeader().setStretchLastSection(True)
        #self.tableWidget.horizontalHeaderItem(0).setText("newHeader")
        self.combo_bound = QComboBox()
        self.combo_source = QComboBox()
        self.combo_dem = QComboBox()
        self.tableWidget.setCellWidget(0,0, self.combo_source)
        self.tableWidget.setCellWidget(1,0, self.combo_bound)
        self.tableWidget.setCellWidget(2,0, self.combo_dem)

        self.tableWidget.setItem(3 , 0, QTableWidgetItem(""))
        self.tableWidget.setItem(4 , 0, QTableWidgetItem(""))
        self.tableWidget.setItem(5 , 0, QTableWidgetItem(""))
        self.tableWidget.setItem(6 , 0, QTableWidgetItem(""))



        hbox = QHBoxLayout()
        hbox.setContentsMargins(0, 0, 0, 0)
        hbox.setSpacing(0)
        self.line_output = QLineEdit()
        self.line_output.setFixedHeight(25)
        self.saveButton = QPushButton("Scegli")
        self.saveButton.setFixedHeight(25)
        hbox.addWidget(self.line_output)
        hbox.addWidget(self.saveButton)
        cellWidget = QWidget()
        cellWidget.setLayout(hbox)


        self.tableWidget.setCellWidget(7,0, cellWidget)

        self.spintime=QSpinBox()
        self.spinRes=QSpinBox()

        self.tableWidget.setCellWidget(8,0, self.spinRes)

        self.spinRes.setValue(25)
        self.spintime.setValue(10)


        self.tableWidget.resizeRowsToContents();

        # rowPosition = self.tableWidget.rowCount()
        # self.tableWidget.insertRow(rowPosition)
        # self.tableWidget.setItem(rowPosition , 0, QtGui.QTableWidgetItem("text1"))
        self.tableWidget.setVerticalHeaderLabels((u'Vettoriale sorgente*', u'Vettoriale confine*', u'DEM*',
                                                  u'Altezza sorgente (m)',u'Altezza osservatore (m)',
                                                  u'Coeff. di rarefrazione',u'Memoria',u'Output file',u'Risoluzione'))



        self.label_status.setText("In attesa di dati")
        self.label_status.setStyleSheet('color : green; font-weight:bold')

        self.clear_out_button.clicked.connect(self.reset_output)
        self.save_out_button.clicked.connect(self.esporta_output)

        # self.web = QWebView()
        # self.web.load(QUrl("https://grass.osgeo.org/grass70/manuals/addons/r.green.biomassfor.theoretical.html"))
        # self.web_layout.addWidget(self.web)

        self.popolacombo()

        self.saveButton.clicked.connect(lambda: self.scegli_file("salvaraster"))
        self.reset_field_button.clicked.connect(self.reset_fields)
        self.buttonBox.accepted.connect(self.run_light)
        self.actionManuale.triggered.connect(self.help)
        self.actionCredits.triggered.connect(self.about)
        self.actionSetting.triggered.connect(self.configuration)


        #pyqtRemoveInputHook()
        #pdb.set_trace()
        self.tabWidget.removeTab(1)
=#
#=
    def esporta_output(self):
        resultmodel=self.console.toPlainText()
        name = QFileDialog.getSaveFileName(self, 'Save File')
        file = open(name[0],'w')
        file.write(resultmodel)
        file.close()



    def reset_output(self):
        ret = QMessageBox.warning(self,"Attenzione", "Vuoi davvero eliminare i risultati del modello?",QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        if ret== QMessageBox.Yes:
            self.console.clear()
        else:
            return False


    def help(self):
        #self.credits = u"Università della Tuscia\n Viterbo - Italy\nRaffaele Pelorosso, Federica Gobattoni\nDeveloper: Francesco Geri"
        #QMessageBox.about(self.dlg,"Credits", self.credits )
        if platform.uname()[0]=="Windows":
            os.system("start "+os.path.dirname(__file__)+"/../tutorial/manuale_envifate_inquinamento_luminoso.pdf")
        if platform.uname()[0]=="Linux":
            os.system("xdg-open "+os.path.dirname(__file__)+"/../tutorial/manuale_envifate_inquinamento_luminoso.pdf")
        else:
            os.system("open "+os.path.dirname(__file__)+"/../tutorial/manuale_envifate_inquinamento_luminoso.pdf")


    def run1(self):
        # fix_print_with_import
        print("run effettuato")

    def popolafields(self,combo_in,combo_out):
        vect_source_text=combo_in.currentText()
        if vect_source_text!="":
            #vfields = self.allLayers[mainvect].pendingFields()
            mainvect = self.registry.mapLayersByName( vect_source_text )[0]
            vfields = mainvect.pendingFields()
            combo_out.addItem("No field")
            for field in vfields:
                combo_out.addItem(field.name())

    def configuration(self):
        d = do_setting.Dialog(self.iface)
        d.show()
        d.exec_()

    def popolacombo(self):
        self.combo_source.clear()
        self.combo_dem.clear()
        self.combo_bound.clear()
        self.line_output.clear()
        self.progressBar.setValue(0)


        self.allLayers = self.canvas.layers()
        self.listalayers=dict()
        #elementovuoto="No required"
        for i in self.allLayers:
            if i.type() == QgsMapLayer.VectorLayer:
                self.listalayers[i.name()]=i
                self.combo_source.addItem(str(i.name()))
                self.combo_bound.addItem(str(i.name()))
            if i.type()==QgsMapLayer.RasterLayer:
                self.listalayers[i.name()]=i
                self.combo_dem.addItem(str(i.name()))





    def reset_fields(self):
        self.console.clear()


        for i in range(3,6):
            self.tableWidget.setItem(i , 0, QTableWidgetItem(""))

        self.popolacombo()

        ##### da eliminare #####
        # self.line_speed.setText("4")
        # self.line_flickianx.setText("1000")
        # self.line_flickiany.setText("1000")
        # self.line_conc.setText("2000")
        # self.line_lambda.setText("0")
        # self.spintime.setValue(20)


        ##### fine da eliminare #####


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
            self.line_output.setText(self.fname[0])
        # if self.tipofile=="salvacsv":
        #     self.fname = QFileDialog.getSaveFileName(None, 'Save file', '/home','csv files (*.csv);;all files (*.*)')
        #     self.dlg.lineEdit_csv.setText(self.fname)


    def about(self):
        QMessageBox.about(self, "Credits EnviFate",u"""<p>EnviFate: Open source tool for environmental risk analysis<br />Release 1.0<br />13-1-2017<br />License: GPL v. 3<br /><a href='https://bitbucket.org/fragit/envifate'>Home page plugin</a></p><hr><p>Lavoro svolto nell’ambito del  Progetto  di   ricerca   scientifica  “Definizione  di   metodi   standard     e  di strumenti applicativi   informatici per   il calcolo degli effetti dei fattori di perturbazione   ai sensi della decisione  2011/484/Ue,  da impiegarsi  nell’ambito  della valutazione di incidenza” finanziato dalla Regione Veneto. Partner principale è il DICAM, Dipartimento di Ingegneria Civile Ambientale e Meccanica dell’Università di Trento (Italia).</p><hr><p>Autori: Francesco Geri, Marco Ciolli</p><p>Universita' di Trento, Trento - Dipartimento di Ingegneria Civile Ambientale e Meccanica (DICAM) <a href="http://www.dicam.unitn.it/">www.dicam.unitn.it/</a></p><hr><p>Consulenti: Paolo Zatelli, Oscar Cainelli</p>""")
=#
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

function run_light( source, boundary, resolution::Int64, src_h::Real=0.0, obsrv_h::Real=1.75, rarefraction::Real=0.14286, processing_memory::Int64=500, output_path::AbstractString=".\\light_intensity.tiff"  )
#=
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
=#



    path_layer=self.areastudio.dataProvider().dataSourceUri()
    path=path_layer.split("|")
    source_ds = ogr.Open(path[0])
    area_layer = source_ds.GetLayer()
    x_min, y_min, x_max, y_max = Int64.([ area_layer.GetExtent()[1], area_layer.GetExtent()[3], area_layer.GetExtent()[2], area_layer.GetExtent()[4] ])
    
    
    

    mem_driver = gdal.GetDriverByName('MEM')
    valNoData = -9999

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
