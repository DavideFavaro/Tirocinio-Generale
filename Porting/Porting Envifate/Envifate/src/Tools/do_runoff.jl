# -*- coding: utf-8 -*-
"""
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
"""
from __future__ import print_function

from builtins import str
from builtins import range
from qgis.PyQt import QtCore, QtGui
from PyQt5.QtCore import QSettings, QTranslator, QCoreApplication, Qt, QObject, pyqtSignal, pyqtRemoveInputHook,QVariant
from PyQt5.QtGui import QIcon
from PyQt5.QtWidgets import QAction, QDialog, QFormLayout, QMenu, QComboBox, QTableWidgetItem, QHBoxLayout, QLineEdit, QPushButton, QWidget, QSpinBox, QTableWidgetItem, QMessageBox, QFileDialog

import datetime
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

import pdb

from envifate_dialog import EnviDialog

from configuration_dialog import ConfigurationDialog

sys.path.append( os.path.dirname(__file__)+"/../library" )

import functions, do_setting

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

        self.label_title.setText("Analisi ruscellamento")
        self.label_title.setStyleSheet('background-color : qlineargradient(spread:pad, x1:0, y1:0, x2:1, y2:0, stop:0 #6b0200, stop:1 rgba(0, 0, 0, 0)); color : white')
        self.tableWidget.setRowCount(9)
        self.tableWidget.horizontalHeader().setStretchLastSection(True)
        #self.tableWidget.horizontalHeaderItem(0).setText("newHeader")
        self.combo_bound = QComboBox()
        self.combo_source = QComboBox()
        self.combo_dem = QComboBox()
        self.combo_lc = QComboBox()
        self.combofield_lc = QComboBox()
        self.combofield_soil = QComboBox()
        self.combo_target = QComboBox()
        self.combofield_target = QComboBox()
        self.combo_fieldp = QComboBox()
        self.tableWidget.setCellWidget(0,0, self.combo_source)
        self.tableWidget.setCellWidget(1,0, self.combo_fieldp)
        self.tableWidget.setCellWidget(2,0, self.combo_bound)


        hbox_lc = QHBoxLayout()
        hbox_lc.setContentsMargins(0, 0, 0, 0)
        hbox_lc.setSpacing(0)
        self.line_soil = QLineEdit()
        self.combo_soil = QComboBox()
        #self.line_freqList.setFixedHeight(25)
        hbox_lc.addWidget(self.combo_lc)
        hbox_lc.addWidget(self.combofield_lc)
        cellWidget_lc = QWidget()
        cellWidget_lc.setLayout(hbox_lc)

        self.tableWidget.setCellWidget(3,0, cellWidget_lc)

        hbox_subs = QHBoxLayout()
        hbox_subs.setContentsMargins(0, 0, 0, 0)
        hbox_subs.setSpacing(0)
        self.line_soil = QLineEdit()
        self.combo_soil = QComboBox()
        #self.line_freqList.setFixedHeight(25)
        hbox_subs.addWidget(self.combo_soil)
        hbox_subs.addWidget(self.combofield_soil)
        cellWidget_subs = QWidget()
        cellWidget_subs.setLayout(hbox_subs)

        self.tableWidget.setCellWidget(4,0, cellWidget_subs)


        self.tableWidget.setCellWidget(5,0, self.combo_dem)


        hbox_target = QHBoxLayout()
        hbox_target.setContentsMargins(0, 0, 0, 0)
        hbox_target.setSpacing(0)
        self.line_target = QLineEdit()
        hbox_target.addWidget(self.combo_target)
        hbox_target.addWidget(self.combofield_target)
        cellWidget_target = QWidget()
        cellWidget_target.setLayout(hbox_target)

        self.tableWidget.setCellWidget(6,0, cellWidget_target)


        #self.tableWidget.setCellWidget(6,0, self.combo_target)



        hbox = QHBoxLayout()
        hbox.setContentsMargins(0, 0, 0, 0)
        hbox.setSpacing(0)
        self.line_folder = QLineEdit()
        self.line_folder.setFixedHeight(25)
        self.saveButton = QPushButton("Scegli")
        self.saveButton.setFixedHeight(25)
        hbox.addWidget(self.line_folder)
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
        self.tableWidget.setVerticalHeaderLabels((u'Vettoriale sorgente*', u'Input quantità*', u'Vettoriale confine*',u'Vettoriale landcover*', u'Cat. suolo',
                                                  u'DEM*',u'Vettoriale target (campo nome)',u'Working folder',u'Risoluzione'))



        self.label_status.setText("In attesa di dati")
        self.label_status.setStyleSheet('color : green; font-weight:bold')

        self.clear_out_button.clicked.connect(self.reset_output)
        self.save_out_button.clicked.connect(self.esporta_output)

        # self.web = QWebView()
        # self.web.load(QUrl("https://grass.osgeo.org/grass70/manuals/addons/r.green.biomassfor.theoretical.html"))
        # self.web_layout.addWidget(self.web)

        self.popolacombo()

        self.combo_lc.currentIndexChanged[str].connect(self.checkfields)
        self.combo_source.currentIndexChanged[str].connect(self.checkfield_source)
        self.combo_target.currentIndexChanged[str].connect(self.checkfield_target)

        self.saveButton.clicked.connect(lambda: self.scegli_file("folder"))
        self.reset_field_button.clicked.connect(self.reset_fields)
        self.buttonBox.accepted.connect(self.run_runoff)
        self.actionManuale.triggered.connect(self.help)
        self.actionCredits.triggered.connect(self.about)
        self.actionSetting.triggered.connect(self.configuration)


        self.classisoil={}

        self.classisoil['A']=0
        self.classisoil['B']=1
        self.classisoil['C']=2
        self.classisoil['D']=3


        #pyqtRemoveInputHook()
        #pdb.set_trace()
        self.tabWidget.removeTab(1)




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
            os.system("start "+os.path.dirname(__file__)+"/../tutorial/manuale_envifate_ruscellamento.pdf")
        if platform.uname()[0]=="Linux":
            os.system("xdg-open "+os.path.dirname(__file__)+"/../tutorial/manuale_envifate_ruscellamento.pdf")
        else:
            os.system("open "+os.path.dirname(__file__)+"/../tutorial/manuale_envifate_ruscellamento.pdf")


    def run1(self):
        # fix_print_with_import
        print("run effettuato")


    def checkfields(self):
        self.combofield_lc.clear()
        self.combofield_soil.clear()
        self.popolafields(self.combo_lc,self.combofield_lc)
        self.popolafields(self.combo_lc,self.combofield_soil)


    def checkfield_source(self):
        self.combo_fieldp.clear()
        self.popolafields(self.combo_source,self.combo_fieldp)


    def checkfield_target(self):
        self.combofield_target.clear()
        self.popolafields(self.combo_target,self.combofield_target)

    def popolafields(self,combo_in,combo_out):
        vect_source_text=combo_in.currentText()
        if vect_source_text!="":
            #vfields = self.allLayers[mainvect].pendingFields()
            mainvect = QgsProject.instance().mapLayersByName( vect_source_text )[0]
            vfields = mainvect.fields()
            #combo_out.addItem("No field")
            for field in vfields:
                combo_out.addItem(field.name())

    def configuration(self):
        d = do_setting.Dialog(self.iface)
        d.show()
        d.exec_()

    def popolacombo(self):
        self.combo_source.clear()
        self.combo_dem.clear()
        self.combo_soil.clear()
        self.combo_bound.clear()
        self.combofield_lc.clear()
        self.combo_target.clear()
        self.combofield_soil.clear()
        self.combofield_target.clear()
        self.combo_fieldp.clear()
        self.combo_lc.clear()
        self.line_folder.clear()
        self.progressBar.setValue(0)


        self.allLayers = self.canvas.layers()
        self.listalayers=dict()
        #elementovuoto="No required"
        for i in self.allLayers:
            if i.type() == QgsMapLayer.VectorLayer:
                self.listalayers[i.name()]=i
                self.combo_source.addItem(str(i.name()))
                self.combo_bound.addItem(str(i.name()))
                self.combo_lc.addItem(str(i.name()))
                self.combo_target.addItem(str(i.name()))
            if i.type()==QgsMapLayer.RasterLayer:
                self.listalayers[i.name()]=i
                self.combo_dem.addItem(str(i.name()))

        self.combo_soil.addItem("Valore campo")
        self.combo_soil.addItem("A")
        self.combo_soil.addItem("B")
        self.combo_soil.addItem("C")
        self.combo_soil.addItem("D")

        self.popolafields(self.combo_lc,self.combofield_lc)
        self.popolafields(self.combo_lc,self.combofield_soil)
        self.popolafields(self.combo_source,self.combo_fieldp)
        self.popolafields(self.combo_target,self.combofield_target)



    def reset_fields(self):
        self.console.clear()



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
        if tipofile=="folder":
            self.folder = QFileDialog.getExistingDirectory(self, "Select Directory")
            self.line_folder.setText(self.folder)
        # if self.tipofile=="salvacsv":
        #     self.fname = QFileDialog.getSaveFileName(None, 'Save file', '/home','csv files (*.csv);;all files (*.*)')
        #     self.dlg.lineEdit_csv.setText(self.fname)


    def about(self):
        QMessageBox.about(self, "Credits EnviFate",u"""<p>EnviFate: Open source tool for environmental risk analysis<br />Release 1.0<br />13-1-2017<br />License: GPL v. 3<br /><a href='https://bitbucket.org/fragit/envifate'>Home page plugin</a></p><hr><p>Lavoro svolto nell’ambito del  Progetto  di   ricerca   scientifica  “Definizione  di   metodi   standard     e  di strumenti applicativi   informatici per   il calcolo degli effetti dei fattori di perturbazione   ai sensi della decisione  2011/484/Ue,  da impiegarsi  nell’ambito  della valutazione di incidenza” finanziato dalla Regione Veneto. Partner principale è il DICAM, Dipartimento di Ingegneria Civile Ambientale e Meccanica dell’Università di Trento (Italia).</p><hr><p>Autori: Francesco Geri, Marco Ciolli</p><p>Universita' di Trento, Trento - Dipartimento di Ingegneria Civile Ambientale e Meccanica (DICAM) <a href="http://www.dicam.unitn.it/">www.dicam.unitn.it/</a></p><hr><p>Consulenti: Paolo Zatelli, Oscar Cainelli</p>""")



    def calc_s(self,cn):
        return(254.0*((100/cn)-1))



    def extract_values(self, raster,x,y):
        z=raster.dataProvider().identify(QgsPointXY(x, y),QgsRaster.IdentifyFormatValue)
        zresult=z.results()
        zvalue=zresult[1]
        return(zvalue)

    def run_runoff(self):
        self.text_vector = str(self.combo_source.currentText())
        self.text_area = str(self.combo_bound.currentText())
        self.text_dem = str(self.combo_dem.currentText())
        self.text_lc = str(self.combo_lc.currentText())
        self.text_lcfield = str(self.combofield_lc.currentText())
        self.text_p = str(self.combo_fieldp.currentText())
        self.text_target = str(self.combo_target.currentText())
        self.text_targetfield = str(self.combofield_target.currentText())
        self.text_lcfield = str(self.combofield_lc.currentText())
        self.text_soil = str(self.combo_soil.currentText())
        self.text_soilfield = str(self.combofield_soil.currentText())

        self.res=int(self.spinRes.text())






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

        self.target=self.listalayers[self.text_target]

        if self.target.wkbType()!=6:
            QMessageBox.warning(self,"Warning", "L'area target deve avere geometria poligonale" )
            return

        self.lc=self.listalayers[self.text_lc]


        if self.lc.wkbType()!=6:
            QMessageBox.warning(self,"Warning", "Not a valid landcover geometry" )
            return

        # self.path_output=self.line_output.text()
        # if self.path_output=="":
        #     self.path_output=os.path.dirname(__file__)+"/runoff.tif"

        if self.areastudio.crs().authid()!=self.source.crs().authid() or self.lc.crs().authid()!=self.source.crs().authid() or self.dem.crs().authid()!=self.source.crs().authid() or self.target.crs().authid()!=self.source.crs().authid():
            QMessageBox.warning(self,"Warning", "Errore: i sistemi di riferimento non sono uniformi. Impossibile continuare con l'analisi." )
            return


        self.refsys=self.source.crs().authid().split(':')[1]


        self.path_working=self.line_folder.text()
        if self.path_working=="":
            self.path_working=os.path.dirname(__file__)

        self.path_temp_lc=self.path_working+"/temp_lc.tif"
        self.path_temp_soil=self.path_working+"/temp_soil.tif"


        #recupero dati database

        listaclc=functions.cn_list_extract()


        controllo_soil=0

        messaggio="Inizio elaborazione analisi dispersione per ruscellamento\n"
        messaggio+="---------------------------\n\n"
        messaggio+="FILE DI INPUT:\n"
        messaggio+="Vettoriale sorgente: "+str(self.text_vector)+"\n"
        messaggio+="Vettoriale confine: "+str(self.text_area)+"\n"
        messaggio+="Vettoriale target: "+str(self.text_target)+"\n"
        messaggio+="DTM: "+str(self.text_dem)+"\n\n"

        messaggio+="VARIABILI:\n"
        messaggio+="Risoluzione: "+str(self.res)+"\n\n"
        messaggio+='ALGORITMO UTILIZZATO: calcolo della separazione delle componenti infiltrazione e ruscellamento tramite metodo SCS-CN; US Department of Agriculture Soil Conservation Service, 1972. National Engineering Handbook, Section 4, Hydrology. US Government Printing Office, Washington, DC, 544pp.\n\n'
        messaggio+="---------------------------\n\n"
        self.console.appendPlainText(messaggio)


        self.label_status.setText("Preparazione dati")
        self.label_status.setStyleSheet('color : #e8b445;font-weight:bold')

        self.path_working=self.line_folder.text()
        if self.path_working=="":
            self.path_working=os.path.dirname(__file__)

        self.path_output=self.path_working+"/runoff.tif"



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
                               'modulo':'Analisi ruscellamento',
                               'descrizione':'Analisi di ruscellamento di un inquinante attraverso il metodo della separazione delle componenti',
                               'srs':self.source.crs().authid(),
                               'data':datetime.datetime.now().strftime("%d-%m-%y")})

        # geotransform = target_ds.GetGeoTransform()


        band = target_ds.GetRasterBand(1)
        band.SetNoDataValue(float(NoData_value))
        band.Fill(NoData_value)
        xsize = band.XSize
        ysize = band.YSize


        intervallo=int(self.dem.rasterUnitsPerPixelX())

        outData = np.array(band.ReadAsArray(0, 0, xsize,ysize).astype(np.float))

        lc_clip_proc = processing.run('qgis:clip', {'INPUT':self.lc, 'OVERLAY':self.areastudio, 'OUTPUT':self.path_working+'/clip.gpkg'})
        lc_clip=QgsVectorLayer(lc_clip_proc['OUTPUT'], 'lc_clip', 'ogr')
        lc_clip.setCrs(self.source.crs())


        #path__layer_lc=lc_clip['OUTPUT'].dataProvider().dataSourceUri()
        path__layer_lc=lc_clip.dataProvider().dataSourceUri()
        path_lc=path__layer_lc.split("|")
        source_ds_lc = ogr.Open(path_lc[0])
        lc_layer = source_ds_lc.GetLayer()
        lc_ds = gdal.GetDriverByName('GTiff').Create(self.path_temp_lc, int(x_res), int(y_res), 1, gdal.GDT_Float32)
        lc_ds.SetGeoTransform((x_min, 25, 0, y_max, 0, -25))
        lc_ds.SetProjection( srs.ExportToWkt() )
        band_lc = lc_ds.GetRasterBand(1)
        band_lc.SetNoDataValue(float(-9999))
        band_lc.Fill(-9999)
        xsize = band_lc.XSize
        ysize = band_lc.YSize
        gdal.RasterizeLayer(lc_ds, [1], lc_layer,options=["ATTRIBUTE="+self.text_lcfield])
        lc_ds=None

        lc_layer=QgsRasterLayer(self.path_temp_lc,"lc_layer")



        if self.text_soil=="Valore campo":
            controllo_soil=1

            source_ds_soil = ogr.Open(path_lc[0])
            soil_layer = source_ds_soil.GetLayer()
            soil_ds = gdal.GetDriverByName('GTiff').Create(self.path_temp_soil, int(x_res), int(y_res), 1, gdal.GDT_Float32)
            soil_ds.SetGeoTransform((x_min, 25, 0, y_max, 0, -25))
            soil_ds.SetProjection( srs.ExportToWkt() )
            band_soil = soil_ds.GetRasterBand(1)
            band_soil.SetNoDataValue(float(-9999))
            band_soil.Fill(-9999)
            gdal.RasterizeLayer(lc_ds, [1], soil_layer,options=["ATTRIBUTE="+self.text_soilfield])
            soil_ds=None

            soil_layer=QgsRasterLayer(self.path_temp_soil,"soil_layer")


        grass_area=str(x_min)+','+str(x_max)+','+str(y_min)+','+str(y_max)+' ['+str(self.areastudio.crs().authid())+']'
        # grass_coord=str(x_source)+','+str(y_source)+' ['+str(self.source.crs().authid())+']'

        namewatershed=self.path_working+'/watershed'
        namedrain=self.path_working+'/wshed.shp'

        params = { 'GRASS_RASTER_FORMAT_OPT' : '','GRASS_REGION_CELLSIZE_PARAMETER' : 0, 'GRASS_REGION_PARAMETER' :grass_area,
                   'GRASS_VECTOR_EXPORT_NOCAT' : False, '-a' : False, 'start_coordinates' : None,  '-n' : False,
                   'input' : self.dem.dataProvider().dataSourceUri(),'-c' : True, 'drain' : namedrain, 'GRASS_MIN_AREA_PARAMETER' : 0.0001,
                   'start_points' : self.source.dataProvider().dataSourceUri(), 'output' : namewatershed }



        waterwshed_proc = processing.run('grass7:r.drain', params)


        #aggiungo per controllo la viewshed alla toc
        #iface.addVectorLayer(namedrain,'watershed','ogr')
        #watershed=QgsProject.instance().mapLayersByName('watershed.shp')


        vdrain = QgsVectorLayer(namedrain, 'vdrain', 'ogr')

        idxlevel = self.source.fields().indexFromName(self.text_p)
        idxcat = vdrain.fields().indexFromName('cat')

        idxtargetname = self.target.fields().indexFromName(self.text_targetfield)



        # pyqtRemoveInputHook()
        # pdb.set_trace()

        features = vdrain.getFeatures()

        nfeat=0

        polygons_t = [feature_t for feature_t in self.target.getFeatures()]


        start_time = time.time()

        for f in features:
            geom = f.geometry()
            length = geom.length()
            currentdistance=intervallo
            nfeat+=1

            featlines=[]


            firstpoint=geom.interpolate(0)

            old_x=firstpoint.asPoint()[0]
            old_y=firstpoint.asPoint()[1]

            fileoutput=self.path_working+'/drain'+str(nfeat)+'.shp'


            max_progress=length/intervallo
            self.progressBar.setMaximum(max_progress)



            vline = QgsVectorLayer("LineString?crs=EPSG:"+self.refsys, "drain"+str(nfeat), "memory")

            prline = vline.dataProvider()
            prlfield=prline.addAttributes( [ QgsField("concentrazione", QVariant.Double) ] )

            self.label_status.setText("Processing data")
            self.label_status.setStyleSheet('color : #e8b445;font-weight:bold')

            idf=f.attributes()[idxcat]
            feat_drain = next(self.source.getFeatures(QgsFeatureRequest().setFilterFid(idf-1)))
            p00=feat_drain.attributes()[idxlevel]

            index_progress=0

            while currentdistance < length:
                if index_progress==0:
                    p0=p00
                else:
                    p0=pe


                point = geom.interpolate(currentdistance)
                x=point.asPoint()[0]
                y=point.asPoint()[1]
                clc=self.extract_values(lc_layer,x,y)
                if controllo_soil==1:
                    soil=self.extract_values(soil_layer,x,y)
                else:
                    soil=self.text_soil

                try:
                    cn=listaclc[clc][self.classisoil[soil]]
                    S=self.calc_s(int(cn))
                except Exception as e:
                    S=0

                pcheck=(math.pow((p0-(0.2*S)),2))/(p0-(0.2*S)+S)
                if pcheck>(0.2*S):
                    pe=pcheck
                    fetline = QgsFeature()
                    fetline.setGeometry( QgsGeometry.fromPolyline( [QgsPoint(old_x,old_y),QgsPoint(x,y)] ))
                    fetline.initAttributes(1)
                    fetline.setAttribute(0,pe)
                    vline.updateFeature(fetline)

                    featlines.append(fetline)
                    index_progress+=1

                    for pol_t in polygons_t:
                        poly_t = pol_t.geometry()
                        if poly_t.contains(point):
                            nometarget=pol_t.attributes()[idxtargetname]
                            messaggio="\nIl vettore drain"+str(nfeat)+" ha raggiunto l'area bersaglio denominata '"+str(nometarget)+"' con un volume pari a: "+str(round(pe,3))+" mm\n"
                            self.console.appendPlainText(messaggio)
                            currentdistance=length+1

                else:
                    pe=0
                    currentdistance=length+1
                self.progressBar.setValue(max_progress)
                # pyqtRemoveInputHook()
                # pdb.set_trace()

                old_x=x
                old_y=y

                currentdistance = currentdistance + intervallo


            prline.addFeatures(featlines)
            vline.updateFields()
            QgsProject.instance().addMapLayer(vline)




        tempoanalisi=time.time() - start_time
        tempostimato=time.strftime("%H:%M:%S", time.gmtime(tempoanalisi))
        messaggio="---------------------------------\n"
        messaggio+="Fine modellazione\n"
        messaggio+="\nTempo di analisi: "+tempostimato+"\n"
        messaggio+="---------------------------------\n\n"
        self.console.appendPlainText(messaggio)

        self.label_status.setText("In attesa di dati")
        self.label_status.setStyleSheet('color : green; font-weight:bold')
        self.progressBar.setValue(max_progress)
