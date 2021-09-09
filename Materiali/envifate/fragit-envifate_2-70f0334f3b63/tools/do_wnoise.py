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
from PyQt5.QtGui import QIcon,QColor
from PyQt5.QtWidgets import QAction, QDialog, QFormLayout, QMenu, QComboBox, QTableWidgetItem, QHBoxLayout, QLineEdit, QPushButton, QWidget, QSpinBox, QTableWidgetItem, QMessageBox, QFileDialog

import os,sys
import time
import platform
import csv

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
import processing

import numpy as np
import math
import struct
import datetime

from osgeo import gdal,ogr,osr

import pdb


from envifate_dialog import EnviDialog

from configuration_dialog import ConfigurationDialog

sys.path.append( os.path.dirname(__file__)+"/../library" )

import functions, noise, do_setting

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
        self.label_title.setText("Analisi rumore in acqua")
        self.label_title.setStyleSheet('background-color : qlineargradient(spread:pad, x1:0, y1:0, x2:1, y2:0, stop:0 #4D4E53, stop:1 rgba(0, 0, 0, 0)); color : white')
        self.tableWidget.setRowCount(10)
        self.tableWidget.horizontalHeader().setStretchLastSection(True)

        self.combo_bound = QComboBox()
        self.combo_source = QComboBox()


        self.tableWidget.setCellWidget(0,0, self.combo_source)
        self.tableWidget.setCellWidget(1,0, self.combo_bound)

        #self.tableWidget.setItem(7 , 0, QTableWidgetItem(""))


        hbox7 = QHBoxLayout()
        hbox7.setContentsMargins(0, 0, 0, 0)
        hbox7.setSpacing(0)
        self.line_freqList = QLineEdit()
        self.line_freqList.setFixedHeight(25)
        #self.line_freqList.setFixedWidth(100)
        self.saveFreqlist = QPushButton("CSV")
        self.saveFreqlist.setFixedHeight(25)
        hbox7.addWidget(self.line_freqList)
        hbox7.addWidget(self.saveFreqlist)
        cellWidget7 = QWidget()
        cellWidget7.setLayout(hbox7)

        # pyqtRemoveInputHook()
        # pdb.set_trace()


        self.tableWidget.setCellWidget(2,0, cellWidget7)



        self.tableWidget.setItem(3 , 0, QTableWidgetItem(""))
        self.tableWidget.setItem(4 , 0, QTableWidgetItem(""))
        self.tableWidget.setItem(5 , 0, QTableWidgetItem(""))
        self.tableWidget.setItem(6 , 0, QTableWidgetItem(""))
        self.tableWidget.setItem(7 , 0, QTableWidgetItem(""))
        #self.tableWidget.setItem(12 , 0, QTableWidgetItem(""))

        # hbox = QHBoxLayout()
        # hbox.setContentsMargins(0, 0, 0, 0)
        # hbox.setSpacing(0)
        # self.line_output = QLineEdit()
        # self.line_output.setFixedHeight(25)
        # self.saveButton = QPushButton("Scegli")
        # self.saveButton.setFixedHeight(25)
        # hbox.addWidget(self.line_output)
        # hbox.addWidget(self.saveButton)
        # cellWidget = QWidget()
        # cellWidget.setLayout(hbox)

        #self.tableWidget.setCellWidget(11,0, cellWidget)

        hbox1 = QHBoxLayout()
        hbox1.setContentsMargins(0, 0, 0, 0)
        hbox1.setSpacing(0)
        self.line_folder = QLineEdit()
        self.line_folder.setFixedHeight(25)
        self.saveFolder = QPushButton("Scegli")
        self.saveFolder.setFixedHeight(25)
        hbox1.addWidget(self.line_folder)
        hbox1.addWidget(self.saveFolder)
        cellWidget1 = QWidget()
        cellWidget1.setLayout(hbox1)

        self.tableWidget.setCellWidget(8,0, cellWidget1)


        self.spinRes=QSpinBox()
        self.spinRes.setValue(20)

        self.tableWidget.setCellWidget(9,0, self.spinRes)




        self.tableWidget.setVerticalHeaderLabels((u'Vettoriale sorgente*', u'Vettoriale confine*',u'Frequenza* (Khz)',u'Profondità* (km)',
                                                  u'Salinità* (ppt)',u'Acidità (pH)* ',u'Temperatura* (C°)',u'Output filename',u'Working folder',u'Risoluzione'))



        self.tableWidget.resizeRowsToContents();

        self.tabWidget.removeTab(1)

        self.menufile=self.menuFile_6





        #self.saveButton.clicked.connect(lambda: self.scegli_file("salvaraster"))
        self.saveFolder.clicked.connect(lambda: self.scegli_file("folder"))
        self.saveFreqlist.clicked.connect(lambda: self.scegli_file("csv"))
        self.reset_field_button.clicked.connect(self.reset_fields)
        self.actionCredits.triggered.connect(self.about)
        self.actionManuale.triggered.connect(self.help)


        #self.actionSetting.triggered.connect(self.configuration)


        self.buttonBox.accepted.connect(self.run_spread)

        self.clear_out_button.clicked.connect(self.reset_output)
        self.save_out_button.clicked.connect(self.esporta_output)

        self.label_status.setText("In attesa di dati")
        self.label_status.setStyleSheet('color : green; font-weight:bold')

        self.popolacombo()




    def run1(self):
        # fix_print_with_import
        # print("run effettuato")
        pass


    def help(self):
        #self.credits = u"Università della Tuscia\n Viterbo - Italy\nRaffaele Pelorosso, Federica Gobattoni\nDeveloper: Francesco Geri"
        #QMessageBox.about(self.dlg,"Credits", self.credits )
        if platform.uname()[0]=="Windows":
            os.system("start "+os.path.dirname(__file__)+"/../tutorial/manuale_envifate_rumore_in_acqua.pdf")
        if platform.uname()[0]=="Linux":
            os.system("xdg-open "+os.path.dirname(__file__)+"/../tutorial/manuale_envifate_rumore_in_acqua.pdf")
        else:
            os.system("open "+os.path.dirname(__file__)+"/../tutorial/manuale_envifate_rumore_in_acqua.pdf")



    # def configuration(self):
    #     d = do_setting.Dialog(self.iface)
    #     d.show()
    #     d.exec_()


    # def checkfields(self):
    #     self.combofield_lc.clear()
    #     self.popolafields(self.combo_lc,self.combofield_lc)

    # def popolafields(self,combo_in,combo_out):
    #     vect_source_text=combo_in.currentText()
    #     if vect_source_text!="":
    #         #vfields = self.allLayers[mainvect].pendingFields()
    #         mainvect = QgsProject.instance().mapLayersByName( vect_source_text )[0]
    #         vfields = mainvect.fields()
    #         # combo_out.addItem("No field")
    #         for field in vfields:
    #             combo_out.addItem(field.name())

    def popolacombo(self):
        self.progressBar.setValue(0)
        self.combo_bound.clear()
        self.combo_source.clear()






        self.allLayers = self.canvas.layers()
        self.listalayers=dict()
        #elementovuoto="No required"
        for i in self.allLayers:
            if i.type() == QgsMapLayer.VectorLayer:
                self.listalayers[i.name()]=i
                self.combo_source.addItem(str(i.name()))
                self.combo_bound.addItem(str(i.name()))



    def reset_fields(self):


        self.console.clear()
        self.tableWidget.setItem(2 , 0, QTableWidgetItem(""))
        self.tableWidget.setItem(3 , 0, QTableWidgetItem(""))
        self.tableWidget.setItem(4 , 0, QTableWidgetItem(""))
        self.tableWidget.setItem(5 , 0, QTableWidgetItem(""))
        self.tableWidget.setItem(6 , 0, QTableWidgetItem(""))
        self.tableWidget.setItem(7 , 0, QTableWidgetItem(""))


        self.line_srid.setText("")



        # for i in range(5,11):
        #     self.tableWidget.setItem(i , 0, QTableWidgetItem(""))

        # self.popolacombo()


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

    def scegli_file(self,tipofile):
        if tipofile=="sqlite":
            self.fname = QFileDialog.getOpenFileName(None, 'Open file', '/home','sqlite3 files (*.sqlite);;all files (*.*)')
            self.dlg_conf.pathtodb.setText(self.fname)
        if tipofile=="csv":
            self.fname = QFileDialog.getOpenFileName(None, 'Open file', '/home','csv files (*.csv);;all files (*.*)')
            self.line_freqList.setText(self.fname[0])
        if tipofile=="tif":
            self.fname = QFileDialog.getOpenFileName(None, 'Open file', '/home','GeoTiff files (*.tif);;all files (*.*)')
            self.dlg_reclass.output_raster_class.setText(self.fname[0])
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


    def run_f1(self, S,T):
        # pyqtRemoveInputHook()
        # pdb.set_trace()

        return (math.sqrt(0.78*(S/35))*(math.exp(T/26)))

    def run_f2(self,T):
        return (42*math.exp(T/17))


    def run_spread(self):

        self.check_freq_list=False




        self.text_vector = str(self.combo_source.currentText())
        self.text_area = str(self.combo_bound.currentText())


        #self.text_line_srid=str(self.line_srid.text())
        self.line_depth=str(self.tableWidget.item(3,0).text())
        self.line_salinity=str(self.tableWidget.item(4,0).text())
        self.line_ph=str(self.tableWidget.item(5,0).text())
        self.line_temperature=str(self.tableWidget.item(6,0).text())
        self.line_output=str(self.tableWidget.item(7,0).text())






        self.source=self.listalayers[self.text_vector]

        if self.source.wkbType()!=1:
            QMessageBox.warning(self,"Warning", "Il vettoriale sorgente deve essere di geometria puntuale" )
            return

        self.areastudio=self.listalayers[self.text_area]


        if self.areastudio.wkbType()!=6:
            QMessageBox.warning(self,"Warning", "Il vettoriale area deve essere di geometria poligonale" )
            return



        try:
            self.depth=int(self.line_depth)
            #self.lf=float(self.text_line_acquifer_depth)
            #self.dz=float(self.text_line_sourcethick)
        except Exception as e:
            QMessageBox.warning(self,"Warning", u"Errore nel dato di profondità (km)" )
            return


        try:
            self.temperature=int(self.line_temperature)
            #self.lf=float(self.text_line_acquifer_depth)
            #self.dz=float(self.text_line_sourcethick)
        except Exception as e:
            QMessageBox.warning(self,"Warning", u"Errore nel dato di temperatura (intero C°)" )
            return

        try:
            self.salinity=int(self.line_salinity)
            #self.lf=float(self.text_line_acquifer_depth)
            #self.dz=float(self.text_line_sourcethick)
        except Exception as e:
            QMessageBox.warning(self,"Warning", u"Errore nel dato di salinità (ppm)" )
            return


        try:
            self.ph=int(self.line_ph)
            #self.lf=float(self.text_line_acquifer_depth)
            #self.dz=float(self.text_line_sourcethick)
        except Exception as e:
            QMessageBox.warning(self,"Warning", u"Errore nel dato di ph" )
            return

        try:
            self.freq=int(self.line_freqList.text())
        except Exception as e:
            if os.path.exists(self.line_freqList.text()):
                self.check_freq_list=True
            else:
                QMessageBox.warning(self,"Warning", "Errore nella frequenza" )
                return

        self.misdist=15 #default value of mueasure distance


        if self.areastudio.crs().authid()!=self.source.crs().authid():
            QMessageBox.warning(self,"Warning", "Errore: i sistemi di riferimento non sono uniformi. Impossibile continuare con l'analisi." )
            return


        self.refsys=self.source.crs().authid().split(':')[1]

        self.path_working=self.line_folder.text()
        if self.path_working=="":
            self.path_working=os.path.dirname(__file__)

        self.path_temp_lc=self.path_working+"/temp_lc.tif"



        self.res=int(self.spinRes.text())

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
