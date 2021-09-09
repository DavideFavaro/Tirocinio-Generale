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
from PyQt5.QtCore import QSettings, QTranslator, QCoreApplication, Qt, QObject, pyqtSignal, pyqtRemoveInputHook
from PyQt5.QtGui import QIcon
from PyQt5.QtWidgets import QAction, QDialog, QFormLayout, QMenu, QComboBox, QTableWidgetItem, QHBoxLayout, QLineEdit, QPushButton, QWidget, QSpinBox, QTableWidgetItem, QMessageBox,QFileDialog

import os,sys
import time, datetime

# Initialize Qt resources from file resources.py
#import resources
import qgis
from qgis.core import *
from qgis.gui import *
from qgis.utils import iface

import numpy as np
import math
import struct
import platform


from osgeo import gdal,ogr,osr


# Import the code for the dialog

#from open_risk_dialog import OpenRiskDialog
import os.path
try:
  import sqlite3
except:
  # fix_print_with_import
  print("librerie per la connessione al database sqlite non trovate")

#from qgis.core import QgsMapLayerRegistry

import pdb

from envifate_dialog import EnviDialog

from configuration_dialog import ConfigurationDialog

sys.path.append( os.path.dirname(__file__)+"/../library" )

import functions,plume, do_setting

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

        self.label_title.setText("Dispersione atmosferica")
        self.label_title.setStyleSheet('background-color : qlineargradient(spread:pad, x1:0, y1:0, x2:1, y2:0, stop:0 #34CDEB, stop:1 rgba(0, 0, 0, 0)); color : black')
        self.tableWidget.setRowCount(16)
        self.tableWidget.horizontalHeader().setStretchLastSection(True)
        #self.tableWidget.horizontalHeaderItem(0).setText("newHeader")
        self.combo_bound = QComboBox()
        self.combo_source = QComboBox()
        self.combo_dem = QComboBox()
        self.combo_maindirwind = QComboBox()
        self.combo_stability = QComboBox()
        self.combo_outdoor = QComboBox()
        self.combo_contaminant = QComboBox()
        self.tableWidget.setCellWidget(0,0, self.combo_source)
        self.tableWidget.setCellWidget(1,0, self.combo_bound)
        self.tableWidget.setCellWidget(2,0, self.combo_dem)
        self.tableWidget.setCellWidget(3,0, self.combo_maindirwind)
        self.tableWidget.setCellWidget(4,0, self.combo_stability)
        self.tableWidget.setCellWidget(5,0, self.combo_outdoor)
        self.tableWidget.setCellWidget(6,0, self.combo_contaminant)

        self.tableWidget.setItem(7 , 0, QTableWidgetItem(""))
        self.tableWidget.setItem(8 , 0, QTableWidgetItem(""))
        self.tableWidget.setItem(9 , 0, QTableWidgetItem(""))
        self.tableWidget.setItem(10 , 0, QTableWidgetItem(""))
        self.tableWidget.setItem(11 , 0, QTableWidgetItem(""))
        self.tableWidget.setItem(12 , 0, QTableWidgetItem(""))
        self.tableWidget.setItem(13 , 0, QTableWidgetItem(""))



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
        self.tableWidget.setCellWidget(14,0, cellWidget)

        self.spinRes=QSpinBox()
        self.spinRes.setValue(25)

        self.tableWidget.setCellWidget(15,0, self.spinRes)

        hboxsrid = QHBoxLayout()
        hboxsrid.setContentsMargins(0, 0, 0, 0)
        hboxsrid.setSpacing(0)
        self.line_srid = QLineEdit()
        self.line_srid.setFixedHeight(25)
        hboxsrid.addWidget(self.line_srid)
        cellWidget1 = QWidget()
        cellWidget1.setToolTip('SRID permessi: 3003, 3004, 32632, 32633, 3857, 4326')
        cellWidget1.setLayout(hboxsrid)
        self.tableWidget.setCellWidget(16,0, cellWidget1)

        #self.tableWidget.setItem(15 , 0, QTableWidgetItem(""))


        self.tableWidget.resizeRowsToContents();

        # rowPosition = self.tableWidget.rowCount()
        # self.tableWidget.insertRow(rowPosition)
        # self.tableWidget.setItem(rowPosition , 0, QtGui.QTableWidgetItem("text1"))
        self.tableWidget.setVerticalHeaderLabels((u'Vettoriale sorgente*', u'Vettoriale confine*', u'Mappa DTM',
                                                  u'Direzione vento ',u'Classe stabilità',u'Classe outdoor',u'Contaminante',
                                                  u'Concentrazione* (Kg/h)',u'Temperatura aria (C°)',u'Altezza punto emissione (m)',
                                                  u'Diametro* (m)',u'Temperatura fumi (C°)',u'Velocità fumi (m/s)',
                                                  u'Velocità vento (m/s)',u'Output file',u'Risoluzione (m)'))



        self.label_status.setText("In attesa di dati")
        self.label_status.setStyleSheet('color : green; font-weight:bold')

        self.clear_out_button.clicked.connect(self.reset_output)
        self.save_out_button.clicked.connect(self.esporta_output)

        self.popolacombo()


        self.saveButton.clicked.connect(lambda: self.scegli_file("salvaraster"))
        self.reset_field_button.clicked.connect(self.reset_fields)
        self.buttonBox.accepted.connect(self.run_plume)
        self.actionManuale.triggered.connect(self.help)
        self.actionCredits.triggered.connect(self.about)
        self.actionSetting.triggered.connect(self.configuration)
        #self.button_run.clicked.connect(self.run1)

        self.classioutdoor={}
        self.classioutdoor['country']='c'
        self.classioutdoor['urban']='u'

        self.classiwind={}
        self.classiwind['N']=0
        self.classiwind['NE']=45
        self.classiwind['E']=90
        self.classiwind['SE']=135
        self.classiwind['S']=180
        self.classiwind['SW']=225
        self.classiwind['W']=270
        self.classiwind['NW']=315


        self.classiwind2={}
        self.classiwind2[0]=180
        self.classiwind2[45]=225
        self.classiwind2[90]=270
        self.classiwind2[135]=315
        self.classiwind2[180]=0
        self.classiwind2[225]=45
        self.classiwind2[270]=90
        self.classiwind2[315]=135

        self.tabWidget.removeTab(1)

        self.list_srid=[3003,3004,32632,32633,3857,4326]



    def run1(self):
        # fix_print_with_import
        print("run effettuato")


    def help(self):
        #self.credits = u"Università della Tuscia\n Viterbo - Italy\nRaffaele Pelorosso, Federica Gobattoni\nDeveloper: Francesco Geri"
        #QMessageBox.about(self.dlg,"Credits", self.credits )
        if platform.uname()[0]=="Windows":
            os.system("start "+os.path.dirname(__file__)+"/../tutorial/manuale_envifate_plume.pdf")
        if platform.uname()[0]=="Linux":
            os.system("xdg-open "+os.path.dirname(__file__)+"/../tutorial/manuale_envifate_plume.pdf")
        else:
            os.system("open "+os.path.dirname(__file__)+"/../tutorial/manuale_envifate_plume.pdf")

    def configuration(self):
        d = do_setting.Dialog(self.iface)
        d.show()
        d.exec_()


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


    def reset_fields(self):


        self.console.clear()


        for i in range(6,12):
            self.tableWidget.setItem(i , 0, QTableWidgetItem(""))


        #self.tableWidget.setItem(15 , 0, QTableWidgetItem(""))
        self.popolacombo()


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



    def run_plume(self):


        self.text_conc=str(self.tableWidget.item(7,0).text())
        self.text_etemp=str(self.tableWidget.item(8,0).text())
        self.text_height=str(self.tableWidget.item(9,0).text())
        self.text_diameter=str(self.tableWidget.item(10,0).text())
        self.text_temp=str(self.tableWidget.item(11,0).text())
        self.text_gspeed=str(self.tableWidget.item(12,0).text())
        self.text_wspeed=str(self.tableWidget.item(13,0).text())


        self.text_vector = str(self.combo_source.currentText())
        self.text_area = str(self.combo_bound.currentText())
        self.text_dem = str(self.combo_dem.currentText())

        self.class_stability=str(self.combo_stability.currentText()).split(" ")
        self.c_stability=self.class_stability[1]

        self.class_outdoor=str(self.combo_outdoor.currentText())
        self.outdoor=self.classioutdoor[self.class_outdoor]

        self.res=int(self.spinRes.text())

        self.x_w=self.classiwind[self.combo_maindirwind.currentText()]

        self.text_contaminant = str(self.combo_contaminant.currentText())

        sostanza=self.inquinanti[self.text_contaminant]

        lst_fields=['sf_ing','sf_inal','iur','rfd_ing','rfd_inal','rfc']

        lst_toxic_msg=['Slope Factor per ingestione','Slope Factor per inalazione','Inhalation Unit Risk','Reference Dose per ingestione','Reference Dose per inalazione','Reference Concentration']

        toxic=functions.substance_extract(sostanza,lst_fields,os.path.dirname(__file__)+"/../library/")


        # self.text_line_srid=str(self.line_srid.text())

        try:
            self.q=float(self.text_conc)
        except Exception as e:
            QMessageBox.warning(self,"Warning", "La concentrazione è obbligatoria" )
            return

        try:
            self.u=float(self.text_wspeed)
        except Exception as e:
            QMessageBox.warning(self,"Warning", "Le velocità del vento è obbligatoria" )
            return

        try:
            self.h_s=float(self.text_height)
        except Exception as e:
            QMessageBox.warning(self,"Warning", "L'altezza del camino è obbligatoria" )
            return


        if self.text_gspeed!='':
            try:
                self.v_s=float(self.text_gspeed)
            except Exception as e:
                QMessageBox.warning(self,"Warning", "Errore nella velocità del gas" )
                return
        else:
            self.v_s=0.0

        if self.text_diameter!='':
            try:
                self.d_s=float(self.text_diameter)
            except Exception as e:
                QMessageBox.warning(self,"Warning", "Errore nel diametro dello stack" )
                return
        else:
            self.d_s=0.0

        if self.text_temp!='':
            try:
                self.t_s=float(self.text_temp)
            except Exception as e:
                QMessageBox.warning(self,"Warning", "Errore nella temperatura dei fumi" )
                return
        else:
            self.t_s=0.0

        if self.text_etemp!='':
            try:
                self.t_a=float(self.text_etemp)
            except Exception as e:
                QMessageBox.warning(self,"Warning", "Errore nella temperatura" )
                return
        else:
            self.t_a=0.0



        # try:
        #     self.srid=int(self.text_line_srid)
        #     if self.srid not in self.list_srid :
        #         QMessageBox.warning(self,"Warning", u"Errore codice srid" )
        #         return
        # except Exception as e:
        #     QMessageBox.warning(self,"Warning", u"Errore codice srid" )
        #     return

        self.source=self.listalayers[self.text_vector]

        if self.source.wkbType()!=1:
            QMessageBox.warning(self,"Warning", "TIl file sorgentedeve avere geometria puntuale" )
            return

        self.areastudio=self.listalayers[self.text_area]


        if self.areastudio.wkbType()!=6:
            QMessageBox.warning(self,"Warning", "Il file area deve avere geometria poligonale" )
            return


        self.dem=self.listalayers[self.text_dem]

        if not self.dem.isValid():
            QMessageBox.warning(self,"Warning", "Il file DEM non è valido" )
            return

        self.path_output=self.line_output.text()
        if self.path_output=="":
            self.path_output=os.path.dirname(__file__)+"/output_model.tif"


        if self.areastudio.crs().authid()!=self.source.crs().authid() or self.dem.crs().authid()!=self.source.crs().authid():
            QMessageBox.warning(self,"Warning", "Errore: i sistemi di riferimento non sono uniformi. Impossibile continuare con l'analisi." )
            return

        self.refsys=self.source.crs().authid().split(':')[1]

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

        self.console.appendPlainText(messaggio)


        self.label_status.setText("Preparazione dati")
        self.label_status.setStyleSheet('color : #e8b445;font-weight:bold')
        #rasterizzazione layer boundaries
        #raster_fn = 'test.tif'

        # pyqtRemoveInputHook()
        # pdb.set_trace()


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
                               'modulo':'Dispersione atmosferica',
                               'descrizione':'Simulazione di dispersione di un inquinante in atmosfera da una sorgente tipo "ciminiera", modello gaussiano di Pasquill-Gifford',
                               'srs':self.source.crs().authid(),
                               'data':datetime.datetime.now().strftime("%d-%m-%y")})

        # if self.srid!="":
        #     srs = osr.SpatialReference()
        #     srs.ImportFromEPSG(self.srid)
        #     target_ds.SetProjection( srs.ExportToWkt() )
        # else:
        #    target_ds.SetProjection(projectionfrom)

        band = target_ds.GetRasterBand(1)
        band.SetNoDataValue(float(NoData_value))
        band.Fill(NoData_value)
        xsize = band.XSize
        ysize = band.YSize

        #outData1 = np.zeros((x_res,y_res), np.float16)
        #outData = np.array(target_ds.GetRasterBand(1).ReadAsArray())
        outData = np.array(band.ReadAsArray(0, 0, xsize,ysize).astype(np.float))

        #array_area = np.array(band.ReadAsArray())
        feature = next(self.source.getFeatures())
        geom = feature.geometry().asPoint()
        x_source=geom[0]
        y_source=geom[1]

        #a = band.ReadAsArray().astype(np.float)

        # a= np.array(band.ReadAsArray())
        # (y_index, x_index) = np.where(a==0)




        rows=ysize-1
        cols=xsize-1


        max_progress=rows*cols
        self.progressBar.setMaximum(max_progress)
        start_time = time.time()

        self.label_status.setText("Processing data")
        self.label_status.setStyleSheet('color : #e8b445;font-weight:bold')

        max_dominio=0.0

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

                    z= self.dem.dataProvider().identify(QgsPointXY(x, y),QgsRaster.IdentifyFormatValue)
                    # catx=x_source-x
                    # caty=y_source-y
                    deltax=x-x_source
                    deltay=y-y_source
                    #dist=math.sqrt(math.pow(deltay,2)+math.pow(deltax,2))
                    #tetha=self.classiwind2[self.x_w]-(math.atan(float(caty)/float(catx))*(180/math.pi))

                    # xvero=int(dist*math.cos(math.radians(tetha)))
                    # yvero=int(dist*math.sin(math.radians(tetha)))
                    xvero=deltax*math.cos(math.radians(self.x_w))-deltay*math.sin(math.radians(self.x_w))
                    yvero=deltax*math.sin(math.radians(self.x_w))+deltay*math.cos(math.radians(self.x_w))


                    # pyqtRemoveInputHook()
                    # pdb.set_trace()


                    #controllo=1
                    #if controllo==1:
                    if yvero>0:
                        zresult=z.results()
                        zvalue=zresult[1]
                        element=plume.plume(self.q,yvero,xvero,zvalue,self.c_stability,self.outdoor,
                                            self.u,self.h_s,self.v_s,self.d_s,self.t_s,self.t_a,self.x_w)
                        sigmay,sigmaz=element.calc_sigma()
                        g1,g2=element.calc_g()
                        hvero=element.calc_h()
                        cfinal=element.calc_C()
                        if cfinal>max_dominio:
                            max_dominio=cfinal
                        outData[row,col]=cfinal
                    else:
                        outData[row,col]=0


                    # row_data.append(cfinal)
                # values.append(row_data)
            # array=np.asarray(values, dtype=np.float16)
            #array = np.array(values, dtype=np.float16)


            self.label_status.setText("Preparazione output")
            self.label_status.setStyleSheet('color : #e8b445;font-weight:bold')

            outData_raster=outData[::-1]
            band.WriteArray(outData_raster)



        # fmt = "<" + ("d" * ysize)

        # for row in range(ysize):
        #     scanline = struct.pack(fmt, *array[row])
        #     band.WriteRaster(0, row, xsize, 1, scanline,buf_xsize=xsize,buf_ysize=1, buf_type=gdal.GDT_Float32)
        band= None
        target_ds = None

        base_raster_name=os.path.basename(self.path_output)
        raster_name=os.path.splitext(base_raster_name)[0]
        self.outputlayer=self.iface.addRasterLayer(self.path_output, raster_name)

        layer=None
        for lyr in list(QgsProject.instance().mapLayers().values()):
            if lyr.name() == raster_name:
                layer = lyr



        functions.applystyle(layer,'gr',0.5)

        max_round_dominio=round(max_dominio*1000000,5)

        
        # renderer = layer.renderer()
        # transparency = renderer.rasterTransparency()
        # ltr = QgsRasterTransparency.TransparentSingleValuePixel()


        # tr_list = []
        # ltr.min = 0  # Or another value
        # ltr.max = 0  # Or another value
        # ltr.percentTransparent = 100  # Or another value

        # tr_list.append(ltr)

        # transparency.setTransparentSingleValuePixelList(tr_list)
        # layer.triggerRepaint()


        # x = QgsRasterTransparency.TransparentSingleValuePixel()
        # x.pixelValue = 0
        # x.transparencyPercent = 100
        # transparency.setTransparentSingleValuePixelList([x])
        # renderer.setRasterTransparency(100)

        # pyqtRemoveInputHook()
        # pdb.set_trace()
        #passo alla funzione writeraster una serie di parametri
        #functions.writeraster(output_file,x_min,y_min,pixel_size,pixel_size,xsize,ysize,array_to_raster)
        tempoanalisi=time.time() - start_time
        tempostimato=time.strftime("%H:%M:%S", time.gmtime(tempoanalisi))
        messaggio="---------------------------------\n"
        messaggio+="Fine modellazione\n"
        messaggio+="\nTempo di analisi: "+tempostimato+"\n"

        max_round_dominio=round(max_dominio*1000000,10)

        messaggio+="\nMassima concentrazione rilevata nel dominio di analisi: "+str(max_round_dominio)+" \u03BCg /  m\u00b3\n\n\n"



        # pyqtRemoveInputHook()
        # pdb.set_trace()

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
