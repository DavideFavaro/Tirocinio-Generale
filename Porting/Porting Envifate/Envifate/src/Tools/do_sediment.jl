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

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar

import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D


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

import functions,sediment, do_setting

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

        self.label_title.setText("Analisi sedimentazione marina")
        self.label_title.setStyleSheet('background-color : qlineargradient(spread:pad, x1:0, y1:0, x2:1, y2:0, stop:0 #9b4003, stop:1 rgba(0, 0, 0, 0)); color : white')
        self.tableWidget.setRowCount(15)
        self.tableWidget.horizontalHeader().setStretchLastSection(True)
        #self.tableWidget.horizontalHeaderItem(0).setText("newHeader")
        self.combo_bound = QComboBox()
        self.combo_source = QComboBox()
        self.combo_maindirwind = QComboBox()
        self.tableWidget.setCellWidget(0,0, self.combo_source)
        self.tableWidget.setCellWidget(1,0, self.combo_bound)
        self.tableWidget.setCellWidget(2,0, self.combo_maindirwind)

        self.tableWidget.setItem(3 , 0, QTableWidgetItem(""))
        self.tableWidget.setItem(4 , 0, QTableWidgetItem(""))
        self.tableWidget.setItem(5 , 0, QTableWidgetItem(""))
        self.tableWidget.setItem(6 , 0, QTableWidgetItem(""))
        self.tableWidget.setItem(7 , 0, QTableWidgetItem(""))
        self.tableWidget.setItem(8 , 0, QTableWidgetItem(""))
        self.tableWidget.setItem(9 , 0, QTableWidgetItem(""))
        self.tableWidget.setItem(10 , 0, QTableWidgetItem(""))
        self.tableWidget.setItem(11 , 0, QTableWidgetItem(""))
        self.tableWidget.setItem(12 , 0, QTableWidgetItem(""))



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
        self.tableWidget.setCellWidget(13,0, cellWidget)

        self.spinRes=QSpinBox()
        self.spinRes.setValue(25)

        self.tableWidget.setCellWidget(14,0, self.spinRes)


        self.tableWidget.resizeRowsToContents();

        # rowPosition = self.tableWidget.rowCount()
        # self.tableWidget.insertRow(rowPosition)
        # self.tableWidget.setItem(rowPosition , 0, QtGui.QTableWidgetItem("text1"))
        self.tableWidget.setVerticalHeaderLabels((u'Vettoriale sorgente*', u'Vettoriale confine*', u'Direzione della corrente (°)',u'Velocità media corrente (m/s)',
                                                  u'Profondità media (m)',u'Coefficiente dispersione X',u'Coefficiente dispersione Y',
                                                  u'Quantità di sedimento dragata (Kg/s)',u'Velocità di sedimentazione media (m/s)',
                                                  u'Periodo di analisi (min)',u'Intervallo di analisi (min)',
                                                  u'Ampiezza oscillazione corrente',u'Ciclo di marea (h)',u'Output file',u'Risoluzione (m)'))



        self.label_status.setText("In attesa di dati")
        self.label_status.setStyleSheet('color : green; font-weight:bold')

        self.clear_out_button.clicked.connect(self.reset_output)
        self.save_out_button.clicked.connect(self.esporta_output)

        self.popolacombo()


        self.saveButton.clicked.connect(lambda: self.scegli_file("salvaraster"))
        self.reset_field_button.clicked.connect(self.reset_fields)
        self.buttonBox.accepted.connect(self.run_sediment)
        self.actionManuale.triggered.connect(self.help)
        self.actionCredits.triggered.connect(self.about)
        self.actionSetting.triggered.connect(self.configuration)


        self.classiwind={}
        self.classiwind['N']=0
        self.classiwind['NE']=45
        self.classiwind['E']=90
        self.classiwind['SE']=135
        self.classiwind['S']=180
        self.classiwind['SW']=225
        self.classiwind['W']=270
        self.classiwind['NW']=315


        # demo data
        #da eliminare
        self.tableWidget.item(3,0).setText("1") #V
        self.tableWidget.item(4,0).setText("13") #h
        self.tableWidget.item(5,0).setText("1") #dx
        self.tableWidget.item(6,0).setText("10") #dy
        self.tableWidget.item(7,0).setText("4") #q
        self.tableWidget.item(8,0).setText("0.0359") #w
        self.tableWidget.item(9,0).setText("100") #final t
        self.tableWidget.item(10,0).setText("1") #dt

        #fine demo data



        #self.tabWidget.removeTab(1)


        self.figure = plt.figure()
        self.canvas_mat = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas_mat, self)
        self.layout_mat.addWidget(self.toolbar)
        self.layout_mat_2.addWidget(self.canvas_mat)

        #self.list_srid=[3003,3004,32632,32633,3857,4326]

    def run1(self):
        # fix_print_with_import
        print("run effettuato")


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




    def reset_fields(self):

        self.console.clear()

        for i in range(3,11):
            self.tableWidget.setItem(i , 0, QTableWidgetItem(""))


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
