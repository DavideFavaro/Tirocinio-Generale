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

        self.label_title.setText("Analisi inquimnento termico corsi d'acqua")
        self.label_title.setStyleSheet('background-color : qlineargradient(spread:pad, x1:0, y1:0, x2:1, y2:0, stop:0 #e00d1f, stop:1 rgba(0, 0, 0, 0)); color : cyan')
        self.tableWidget.setRowCount(6)
        self.tableWidget.horizontalHeader().setStretchLastSection(True)
        #self.tableWidget.horizontalHeaderItem(0).setText("newHeader")
        self.combo_source = QComboBox()
        self.combo_dem = QComboBox()
        self.combo_source_t = QComboBox()
        self.combo_source_q = QComboBox()
        self.combo_river = QComboBox()
        self.combo_river_t = QComboBox()
        self.combo_river_q = QComboBox()
        self.tableWidget.setCellWidget(0,0, self.combo_source)

        hbox_lc = QHBoxLayout()
        hbox_lc.setContentsMargins(0, 0, 0, 0)
        hbox_lc.setSpacing(0)
        self.line_soil = QLineEdit()
        self.combo_soil = QComboBox()
        #self.line_freqList.setFixedHeight(25)
        hbox_lc.addWidget(self.combo_source_t)
        hbox_lc.addWidget(self.combo_source_q)
        cellWidget_lc = QWidget()
        cellWidget_lc.setLayout(hbox_lc)
        self.tableWidget.setCellWidget(1,0, cellWidget_lc)



        self.tableWidget.setCellWidget(2,0, self.combo_river)


        hbox_lc = QHBoxLayout()
        hbox_lc.setContentsMargins(0, 0, 0, 0)
        hbox_lc.setSpacing(0)
        self.line_soil = QLineEdit()
        self.combo_soil = QComboBox()
        #self.line_freqList.setFixedHeight(25)
        hbox_lc.addWidget(self.combo_river_t)
        hbox_lc.addWidget(self.combo_river_q)
        cellWidget_lc = QWidget()
        cellWidget_lc.setLayout(hbox_lc)

        self.tableWidget.setCellWidget(3,0, cellWidget_lc)


        self.tableWidget.setCellWidget(4,0, self.combo_dem)



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


        self.tableWidget.setCellWidget(5,0, cellWidget)



        self.tableWidget.resizeRowsToContents();

        # rowPosition = self.tableWidget.rowCount()
        # self.tableWidget.insertRow(rowPosition)
        # self.tableWidget.setItem(rowPosition , 0, QtGui.QTableWidgetItem("text1"))
        self.tableWidget.setVerticalHeaderLabels((u'Vettoriale sorgente*', u'Sorgente campi portata/temperatura*', u'Vettoriale fiume*',u'Fiume campi portata/temperatura*',
                                                  u'DEM*',u'Working folder'))



        self.label_status.setText("In attesa di dati")
        self.label_status.setStyleSheet('color : green; font-weight:bold')

        self.clear_out_button.clicked.connect(self.reset_output)
        self.save_out_button.clicked.connect(self.esporta_output)

        # self.web = QWebView()
        # self.web.load(QUrl("https://grass.osgeo.org/grass70/manuals/addons/r.green.biomassfor.theoretical.html"))
        # self.web_layout.addWidget(self.web)

        self.popolacombo()

        self.combo_source.currentIndexChanged[str].connect(self.checkfield_source)
        self.combo_river.currentIndexChanged[str].connect(self.checkfield_river)

        self.saveButton.clicked.connect(lambda: self.scegli_file("folder"))
        self.reset_field_button.clicked.connect(self.reset_fields)
        self.buttonBox.accepted.connect(self.run_thermic)
        self.actionManuale.triggered.connect(self.help)
        self.actionCredits.triggered.connect(self.about)
        self.actionSetting.triggered.connect(self.configuration)



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


    def checkfield_source(self):
        self.combo_source_t.clear()
        self.combo_source_q.clear()
        self.popolafields(self.combo_source,self.combo_source_t)
        self.popolafields(self.combo_source,self.combo_source_q)

    def checkfield_river(self):
        self.combo_river_t.clear()
        self.combo_river_q.clear()
        self.popolafields(self.combo_river,self.combo_river_t)
        self.popolafields(self.combo_river,self.combo_river_q)


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
        self.combo_source.clear()
        self.combo_river.clear()
        self.combo_source_t.clear()
        self.combo_source_q.clear()
        self.combo_river_t.clear()
        self.combo_river_q.clear()
        self.line_folder.clear()
        self.progressBar.setValue(0)


        self.allLayers = self.canvas.layers()
        self.listalayers=dict()
        #elementovuoto="No required"
        for i in self.allLayers:
            if i.type() == QgsMapLayer.VectorLayer:
                self.listalayers[i.name()]=i
                self.combo_source.addItem(str(i.name()))
                self.combo_river.addItem(str(i.name()))
            if i.type()==QgsMapLayer.RasterLayer:
                self.listalayers[i.name()]=i
                self.combo_dem.addItem(str(i.name()))


        self.popolafields(self.combo_source,self.combo_source_t)
        self.popolafields(self.combo_source,self.combo_source_q)
        self.popolafields(self.combo_river,self.combo_river_t)
        self.popolafields(self.combo_river,self.combo_river_q)



    def reset_fields(self):
        self.console.clear()
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
        if tipofile=="folder":
            self.folder = QFileDialog.getExistingDirectory(self, "Select Directory")
            self.line_folder.setText(self.folder)
        # if self.tipofile=="salvacsv":
        #     self.fname = QFileDialog.getSaveFileName(None, 'Save file', '/home','csv files (*.csv);;all files (*.*)')
        #     self.dlg.lineEdit_csv.setText(self.fname)


    def about(self):
        QMessageBox.about(self, "Credits EnviFate",u"""<p>EnviFate: Open source tool for environmental risk analysis<br />Release 1.0<br />13-1-2017<br />License: GPL v. 3<br /><a href='https://bitbucket.org/fragit/envifate'>Home page plugin</a></p><hr><p>Lavoro svolto nell’ambito del  Progetto  di   ricerca   scientifica  “Definizione  di   metodi   standard     e  di strumenti applicativi   informatici per   il calcolo degli effetti dei fattori di perturbazione   ai sensi della decisione  2011/484/Ue,  da impiegarsi  nell’ambito  della valutazione di incidenza” finanziato dalla Regione Veneto. Partner principale è il DICAM, Dipartimento di Ingegneria Civile Ambientale e Meccanica dell’Università di Trento (Italia).</p><hr><p>Autori: Francesco Geri, Marco Ciolli</p><p>Universita' di Trento, Trento - Dipartimento di Ingegneria Civile Ambientale e Meccanica (DICAM) <a href="http://www.dicam.unitn.it/">www.dicam.unitn.it/</a></p><hr><p>Consulenti: Paolo Zatelli, Oscar Cainelli</p>""")





    def run_thermic(self):
        self.text_vector = str(self.combo_source.currentText())
        self.text_river = str(self.combo_river.currentText())
        self.text_dem = str(self.combo_dem.currentText())
        self.text_source_q = str(self.combo_source_q.currentText())
        self.text_source_t = str(self.combo_source_t.currentText())
        self.text_river_q = str(self.combo_river_q.currentText())
        self.text_river_t = str(self.combo_river_t.currentText())



        # wkbType: 1:point, 6:multipolygon, 2: Linestring

        self.dem=self.listalayers[self.text_dem]

        if not self.dem.isValid():
            QMessageBox.warning(self,"Warning", "The dem file is not valid" )
            return

        self.source=self.listalayers[self.text_vector]

        if self.source.wkbType()!=1:
            QMessageBox.warning(self,"Warning", "The source file must have point geometry" )
            return

        self.river=self.listalayers[self.text_river]



        # pyqtRemoveInputHook()
        # pdb.set_trace()



        if self.river.wkbType()!=2 and self.river.wkbType()!=5:
            QMessageBox.warning(self,"Warning", "Il vettoriale dei fiumi deve avere geometria lineare" )
            return

        # self.path_output=self.line_output.text()
        # if self.path_output=="":
        #     self.path_output=os.path.dirname(__file__)+"/runoff.tif"

        if self.river.crs().authid()!=self.source.crs().authid() or self.dem.crs().authid()!=self.source.crs().authid():
            QMessageBox.warning(self,"Warning", "Errore: i sistemi di riferimento non sono uniformi. Impossibile continuare con l'analisi." )
            return


        self.refsys=self.source.crs().authid().split(':')[1]


        self.path_working=self.line_folder.text()
        if self.path_working=="":
            self.path_working=os.path.dirname(__file__)




        messaggio="Inizio elaborazione analisi inquinamento termico dei corsi d'acqua\n"
        messaggio+="---------------------------\n\n"
        messaggio+="FILE DI INPUT:\n"
        messaggio+="Vettoriale sorgente: "+str(self.text_vector)+"\n"
        messaggio+="Vettoriale fiume: "+str(self.text_river)+"\n"
        messaggio+="DTM: "+str(self.text_dem)+"\n\n"
        messaggio+="---------------------------\n\n"
        self.console.appendPlainText(messaggio)


        self.label_status.setText("Preparazione dati")
        self.label_status.setStyleSheet('color : #e8b445;font-weight:bold')

        self.path_working=self.line_folder.text()
        if self.path_working=="":
            self.path_working=os.path.dirname(__file__)


        start_time = time.time()

        # idxsq = self.source.fields().indexFromName(self.text_source_q)
        # idxst = self.source.fields().indexFromName(self.text_source_t)
        #
        # idxrq = self.source.fields().indexFromName(self.text_river_q)
        # idxrt = self.source.fields().indexFromName(self.text_river_t)


        buffer_sorgente=processing.run('native:buffer', {"INPUT": self.source, "DISTANCE": 5, "OUTPUT":'memory:'})
        sorgente = buffer_sorgente['OUTPUT']


        for f in sorgente.getFeatures():
            for a in self.river.getFeatures():
                if a.geometry().intersects(f.geometry()):
                    intersection = a.geometry().intersection(f.geometry())
                    nuova_portata = { a.fieldNameIndex(self.text_river_q): f[self.text_source_q]}
                    nuova_temp = { a.fieldNameIndex(self.text_river_t): f[self.text_source_t]}
                    self.river.dataProvider().changeAttributeValues({a.id(): nuova_portata})
                    self.river.dataProvider().changeAttributeValues({a.id(): nuova_temp})
                    break #only one or less intersection are possible



        param_point_intersect_all={ 'INPUT_FIELDS' : [], 'OUTPUT' : 'memory:', 'INTERSECT' : self.river, 'INTERSECT_FIELDS' : [], 'INPUT' : self.river }
        point_intersect_all_result=processing.run('native:lineintersections',param_point_intersect_all)
        point_intersect_all=point_intersect_all_result['OUTPUT']

        point_intersect_result=processing.run('qgis:deleteduplicategeometries',{ 'OUTPUT' : 'memory:', 'INPUT' : point_intersect_all })
        point_intersect=point_intersect_result['OUTPUT']

        features = point_intersect.getFeatures()





        count=1
        qmain=0
        diz={}

        for f in features:
            x = f.geometry().asPoint().x()
            y = f.geometry().asPoint().y()
            z = self.dem.dataProvider().identify(QgsPointXY(x, y),QgsRaster.IdentifyFormatValue)
            zresult=z.results()
            diz[zresult[1]]=[f[self.text_river_q],f[self.text_river_t],f[self.text_river_q+"_2"],f[self.text_river_q+"_2"]]

        #print('dizionario')
        #print(diz)

        for key in sorted(diz.keys(),reverse=True):
            if count==1:
                if diz[key][0]>=diz[key][2]:
                    qmain=diz[key][0]
                    q=diz[key][0]
                    t=diz[key][1]
                else:
                    qmain=diz[key][2]
                    q=diz[key][2]
                    t=diz[key][3]

            if diz[key][0]==qmain:
                q1=q
                t1=t
                q2=diz[key][2]
                t2=diz[key][3]
            else:
                q2=q
                t2=t
                q1=diz[key][0]
                t1=diz[key][1]


            q=q1+q2
            t=((q1*t1)+(q2*t2))/(q1+q2)
            count+=1


        messaggio="Temperatura stimata al nodo finale "+str(t)+"\n"
        self.console.appendPlainText(messaggio)

        tempoanalisi=time.time() - start_time
        tempostimato=time.strftime("%H:%M:%S", time.gmtime(tempoanalisi))
        messaggio="---------------------------------\n"
        messaggio+="Fine modellazione\n"
        # messaggio+="\nTempo di analisi: "+tempostimato+"\n"
        messaggio+="---------------------------------\n\n"
        self.console.appendPlainText(messaggio)

        self.label_status.setText("In attesa di dati")
        self.label_status.setStyleSheet('color : green; font-weight:bold')


        # tempoanalisi=time.time() - start_time
        # tempostimato=time.strftime("%H:%M:%S", time.gmtime(tempoanalisi))
        # messaggio="---------------------------------\n"
        # messaggio+="Fine modellazione\n"
        # messaggio+="\nTempo di analisi: "+tempostimato+"\n"
        # messaggio+="---------------------------------\n\n"
        # self.console.appendPlainText(messaggio)

        # self.label_status.setText("In attesa di dati")
        # self.label_status.setStyleSheet('color : green; font-weight:bold')
        # self.progressBar.setValue(max_progress)
