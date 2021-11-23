# -*- coding: utf-8 -*-
"""
/***************************************************************************
 Envifate
                                 A QGIS plugin
 Envifate: Open source tool for environmental risk analysis
                              -------------------
        begin                : 2016-07-15
        git sha              : $Format:%H$
        copyright            : (C) 2016 by Francesco Geri
        email                : francescogeri@tim.com
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
from PyQt5.QtGui import QIcon, QStandardItem,QStandardItemModel
from PyQt5.QtWidgets import QAction, QDialog, QFormLayout, QMenu, QComboBox, QTableWidgetItem, QHBoxLayout, QLineEdit, QPushButton, QWidget, QSpinBox, QTableWidget, QTableWidgetItem, QMessageBox, QFileDialog,QListView

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
import struct
import time, datetime
import platform


from osgeo import gdal,ogr,osr

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar

import matplotlib.pyplot as plt


import pdb


from envifate_dialog import EnviDialog

from configuration_dialog import ConfigurationDialog

sys.path.append( os.path.dirname(__file__)+"/../library" )

import functions, leaching, daf, do_setting

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

        self.label_title.setText("Dispersione in falda")
        self.label_title.setStyleSheet('background-color : qlineargradient(spread:pad, x1:0, y1:0, x2:1, y2:0, stop:0 #CC7D31, stop:1 rgba(0, 0, 0, 0)); color : black')
        self.tableWidget.setRowCount(17)
        self.tableWidget.horizontalHeader().setStretchLastSection(True)
        #self.tableWidget.horizontalHeaderItem(0).setText("newHeader")
        self.combo_bound = QComboBox()
        self.combo_source = QComboBox()
        self.combo_contaminant = QComboBox()
        self.combo_texture = QComboBox()
        self.combo_alg = QComboBox()
        self.tableWidget.setCellWidget(0,0, self.combo_source)
        self.tableWidget.setCellWidget(1,0, self.combo_bound)
        #self.tableWidget.setCellWidget(2,0, self.combo_contaminant)

        hbox_subs = QHBoxLayout()
        hbox_subs.setContentsMargins(0, 0, 0, 0)
        hbox_subs.setSpacing(0)
        self.line_freqList = QLineEdit()
        self.combo_contaminant = QComboBox()
        #self.line_freqList.setFixedHeight(25)
        self.buttonmultiple = QPushButton("TAB")
        self.buttonmultiple.setFixedHeight(25)
        hbox_subs.addWidget(self.combo_contaminant)
        hbox_subs.addWidget(self.buttonmultiple)
        cellWidget_subs = QWidget()
        cellWidget_subs.setLayout(hbox_subs)


        self.tableWidget.setCellWidget(2,0, cellWidget_subs)

        self.tableWidget.setCellWidget(3,0, self.combo_texture)
        self.tableWidget.setCellWidget(4,0, self.combo_alg)

        self.tableWidget.setItem(5 , 0, QTableWidgetItem(""))
        self.tableWidget.setItem(6 , 0, QTableWidgetItem(""))
        self.tableWidget.setItem(7 , 0, QTableWidgetItem(""))
        self.tableWidget.setItem(8 , 0, QTableWidgetItem(""))

        self.tableWidget.setItem(9 , 0, QTableWidgetItem(""))
        self.tableWidget.setItem(10 , 0, QTableWidgetItem(""))
        self.tableWidget.setItem(11 , 0, QTableWidgetItem(""))
        self.tableWidget.setItem(12 , 0, QTableWidgetItem(""))

        self.tableWidget.setItem(13 , 0, QTableWidgetItem(""))
        self.tableWidget.setItem(14 , 0, QTableWidgetItem(""))
        # self.tableWidget.setItem(15 , 0, QTableWidgetItem(""))
        # self.tableWidget.setItem(16 , 0, QTableWidgetItem(""))


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


        self.tableWidget.setCellWidget(15,0, cellWidget)

        #self.spinTime=QtGui.QSpinBox()
        self.spinRes=QSpinBox()

        #self.tableWidget.setCellWidget(17,0, self.spinTime)
        self.tableWidget.setCellWidget(16,0, self.spinRes)

        self.spinRes.setValue(10)
        #self.spinTime.setValue(10)


        self.countmatrix_substance=0;
        self.multiplesubstance_control=False

        #lista sostanze
        self.listasostanze=[]

        #lista concentrazioni
        self.listaconc=[]

        self.tableWidget.resizeRowsToContents();

        # rowPosition = self.tableWidget.rowCount()
        # self.tableWidget.insertRow(rowPosition)
        # self.tableWidget.setItem(rowPosition , 0, QtGui.QTableWidgetItem("text1"))
        # self.tableWidget.setVerticalHeaderLabels((u'Vettoriale sorgente*', u'Vettoriale area*', u'Contaminante',
        #                                           u'Tessitura*',u'Metodo*',u'Estensione (m)',u'Densità (g/cm3)',u'Profondità sorgente secondaria* (cm)',
        #                                           u'Concentrazione inquinante* (kg/l)',u'Profondità falda* (cm)',u'Velocità di Darcy* (m/h)',
        #                                           u'Profondità miscelazione* (cm)',u'Indice dispersione X',u'Indice dispersione Y',
        #                                           u'Indice di 1° decadimento',u'Direzione falda* (°)',u'Precipitazioni medie* (mm/y)',
        #                                           u'Output file',u'Risoluzione*'))
        self.tableWidget.setVerticalHeaderLabels((u'Vettoriale sorgente*', u'Vettoriale confine*', u'Contaminante',
                                                  u'Tessitura*',u'Metodo*',u'Estensione (m)',u'Densità suolo (g/cm3)',u'Profondità sorgente secondaria* (cm)',
                                                  u'Concentrazione inquinante* (kg/l)',u'Profondità falda* (cm)',u'Velocità di Darcy (m/h)',
                                                  u'Profondità miscelazione* (cm)',u'Coeff. di decadimento 1° ordine',u'Direzione falda* (°)',u'Infiltrazione efficace* (mm/a)',
                                                  u'Output file',u'Risoluzione (m)*',u'SRID'))




        self.label_status.setText("In attesa di dati")
        self.label_status.setStyleSheet('color : green; font-weight:bold')

        self.popolacombo()

        # pyqtRemoveInputHook()
        # pdb.set_trace()

        self.figure = plt.figure()
        self.canvas_mat = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas_mat, self)
        self.layout_mat.addWidget(self.toolbar)
        self.layout_mat_2.addWidget(self.canvas_mat)



        self.saveButton.clicked.connect(lambda: self.scegli_file("salvaraster"))
        self.reset_field_button.clicked.connect(self.reset_fields)
        self.buttonBox.accepted.connect(self.leach)
        self.actionManuale.triggered.connect(self.help)
        self.actionCredits.triggered.connect(self.about)

        self.actionSetting.triggered.connect(self.configuration)

        #self.actionInput.triggered.connect(self.matrice_iter_daf)

        self.clear_out_button.clicked.connect(self.reset_output)
        self.save_out_button.clicked.connect(self.esporta_output)
        self.buttonmultiple.clicked.connect(self.matrice_iter_daf)

        ##### per la scelta multipla del modello decommentare le righe qui di seguito
        ##### decommentare anche la funzione self_combo_opz
        ##### vedi anche riga 304 con commento "in caso di scelta tra l'algoritmo fickian o domenico"
        #self.combo_opz.addItem("impulso")
        #self.combo_opz.addItem("continuo")
        #self.connect(self.combo_alg, QtCore.SIGNAL("currentIndexChanged(const QString&)"), self.set_combo_opz)



        self.combo_alg.addItem("Domenico/Schwartz")
        #self.combo_alg.addItem("Fickiano")


        ####### fine


        #self.button_run.clicked.connect(self.run1)


        self.classiwind={}
        self.classiwind['S']=0
        self.classiwind['SW']=45
        self.classiwind['W']=90
        self.classiwind['NW']=135
        self.classiwind['N']=180
        self.classiwind['NE']=225
        self.classiwind['E']=270
        self.classiwind['SE']=315

        self.classiwind2={}
        self.classiwind2[0]=90
        self.classiwind2[45]=135
        self.classiwind2[90]=180
        self.classiwind2[135]=225
        self.classiwind2[180]=270
        self.classiwind2[225]=315
        self.classiwind2[270]=0
        self.classiwind2[315]=45

        self.list_srid=[3003,3004,32632,32633,3857,4326]


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


    def run1(self):
        # fix_print_with_import
        print("run effettuato")


    # def set_combo_opz(self):
    #     check_alg = str(self.combo_alg.currentText())
    #     if check_alg=="Domenico/Schwartz":
    #         self.combo_opz.setEnabled(False)

    #     else:
    #         self.combo_opz.setEnabled(True)



    def help(self):
        #self.credits = u"Università della Tuscia\n Viterbo - Italy\nRaffaele Pelorosso, Federica Gobattoni\nDeveloper: Francesco Geri"
        #QMessageBox.about(self.dlg,"Credits", self.credits )
        if platform.uname()[0]=="Windows":
            os.system("start "+os.path.dirname(__file__)+"/../tutorial/modello_dispersione_in_falda.pdf")
        if platform.uname()[0]=="Linux":
            os.system("xdg-open "+os.path.dirname(__file__)+"/../tutorial/modello_dispersione_in_falda.pdf")
        else:
            os.system("open "+os.path.dirname(__file__)+"/../tutorial/modello_dispersione_in_falda.pdf")

    def configuration(self):
        self.d = do_setting.Dialog(self.iface)

        self.listView=QListView(self)
        self.d.horizontalLayout.addWidget(self.listView)
        self.modelConf = QStandardItemModel()

        self.listView.setModel(self.modelConf)
        servizi=['https://idt2.regione.veneto.it/geoportal/csw','https://idt2-geoserver.regione.veneto.it/geoserver/ows','https://idt2.regione.veneto.it/gwc/service/wmts'];

        for i in servizi:
            item = QStandardItem(i)
            item.setCheckable(True)
            self.modelConf.appendRow(item)

        self.AddButton = QPushButton('Add', self)

        self.d.VerticalLayout.addWidget(self.AddButton)

        self.d.show()
        self.d.exec_()


    def matrice_iter_daf(self):
        self.dmatrix = do_setting.Dialog(self.iface)

        self.tablesubstance = QTableWidget()

        self.tablesubstance.setRowCount(1)
        self.tablesubstance.setColumnCount(2)

        self.dmatrix.horizontalLayout.addWidget(self.tablesubstance)

        self.AddSubstance = QPushButton('+', self)
        self.RemSubstance = QPushButton('-', self)
        self.ReadSubstance = QPushButton('Leggi', self)

        self.AddSubstance.clicked.connect(self.aggiungiriga)
        self.RemSubstance.clicked.connect(self.togliriga)
        self.ReadSubstance.clicked.connect(self.leggisostanza)




        self.dmatrix.VerticalLayout.addWidget(self.AddSubstance)
        self.dmatrix.VerticalLayout.addWidget(self.RemSubstance)
        self.dmatrix.VerticalLayout.addWidget(self.ReadSubstance)
        combo_inq_matrix = QComboBox()
        for row in self.sql_fetch_inq:
            combo_inq_matrix.addItem(row[1])

        #self.tablesubstance.insertRow(self.tablesubstance.rowCount())

        self.tablesubstance.setCellWidget(self.countmatrix_substance,0, combo_inq_matrix)
        self.countmatrix_substance+=1


        self.tablesubstance.setHorizontalHeaderLabels((u'Sostanza', u'Concentrazione'))

        self.tablesubstance.horizontalHeader().setStretchLastSection(True)

        self.combo_contaminant.setEnabled(False)
        self.multiplesubstance_control=True

        self.tableWidget.item(8,0).setFlags( Qt.ItemIsSelectable |  Qt.ItemIsEnabled )
        self.tableWidget.item(8,0).setBackground(QtGui.QColor(166,166,166))

        # pyqtRemoveInputHook()
        # pdb.set_trace()

        self.dmatrix.show()
        self.dmatrix.exec_()



    def aggiungiriga(self):
        combo_inq_matrix = QComboBox()
        for row in self.sql_fetch_inq:
            combo_inq_matrix.addItem(row[1])
        self.tablesubstance.insertRow(self.tablesubstance.rowCount())
        self.tablesubstance.setCellWidget(self.countmatrix_substance,0, combo_inq_matrix)
        self.countmatrix_substance+=1
        # pyqtRemoveInputHook()
        # pdb.set_trace()


    def togliriga(self):
        if self.countmatrix_substance>1:
            self.tablesubstance.removeRow(1)
            self.countmatrix_substance-=1


    def leggisostanza(self):

        self.listasostanze=[]
        for nsostanza in range(0,self.tablesubstance.rowCount()):
            #self.diz_perm[self.table.item(self.c,0).text()]=self.table.item(self.c,1).text()
            try:
                id_substance=self.inquinanti[str(self.tablesubstance.cellWidget(nsostanza,0).currentText())]
                self.listasostanze.append(id_substance)
                self.listaconc.append(float(self.tablesubstance.item(nsostanza,1).text())*1000)
            except:
                QMessageBox.warning(self,"Warning", "Errore nella lista delle concentrazioni" )
                return



        self.dmatrix.close()


    def about(self):
        QMessageBox.about(self, "Credits EnviFate",u"""<p>EnviFate: Open source tool for environmental risk analysis<br />Release 1.0<br />13-1-2017<br />License: GPL v. 3<br /><a href='https://bitbucket.org/fragit/envifate'>Home page plugin</a></p><hr><p>Lavoro svolto nell’ambito del  Progetto  di   ricerca   scientifica  “Definizione  di   metodi   standard     e  di strumenti applicativi   informatici per   il calcolo degli effetti dei fattori di perturbazione   ai sensi della decisione  2011/484/Ue,  da impiegarsi  nell’ambito  della valutazione di incidenza” finanziato dalla Regione Veneto. Partner principale è il DICAM, Dipartimento di Ingegneria Civile Ambientale e Meccanica dell’Università di Trento (Italia).</p><hr><p>Autori: Francesco Geri, Marco Ciolli</p><p>Universita' di Trento, Trento - Dipartimento di Ingegneria Civile Ambientale e Meccanica (DICAM) <a href="http://www.dicam.unitn.it/">www.dicam.unitn.it/</a></p><hr><p>Consulenti: Paolo Zatelli, Oscar Cainelli</p>""")


    def popolafields(self,combo_in,combo_out):
        vect_source_text=combo_in.currentText()
        if vect_source_text!="":
            #vfields = self.allLayers[mainvect].pendingFields()
            mainvect = QgsProject.instance().mapLayersByName( vect_source_text )[0]
            vfields = mainvect.pendingFields()
            combo_out.addItem("No field")
            for field in vfields:
                combo_out.addItem(field.name())

    def popolacombo(self):
        self.progressBar.setValue(0)
        self.combo_source.clear()
        self.combo_bound.clear()
        # self.combo_source_field.clear()
        self.combo_contaminant.clear()
        self.combo_texture.clear()
        self.line_output.clear()

        self.allLayers = self.canvas.layers()
        self.listalayers=dict()
        #combo_bound
        #elementovuoto="No required"
        for i in self.allLayers:
            if i.type() == QgsMapLayer.VectorLayer:
                self.listalayers[i.name()]=i
                self.combo_source.addItem(str(i.name()))
                self.combo_bound.addItem(str(i.name()))

        # self.popolafields(self.combo_source,self.combo_source_field)
        # pyqtRemoveInputHook()
        # pdb.set_trace()
        conn = sqlite3.connect(os.path.dirname(__file__)+"/../library/substance.db")
        cursor=conn.cursor()
        query_substance="select id,nome from substance"
        cursor.execute(query_substance)
        self.sql_fetch_inq=cursor.fetchall()

        self.inquinanti=dict()
        for row in self.sql_fetch_inq:
            self.inquinanti[row[1]]=row[0]
            self.combo_contaminant.addItem(row[1])

        query_texture="select id,nome from texture"
        cursor.execute(query_texture)
        sql_texture=cursor.fetchall()
        #pyqtRemoveInputHook()
        #pdb.set_trace()
        self.texture=dict()
        for rowt in sql_texture:
            self.texture[rowt[1]]=rowt[0]
            self.combo_texture.addItem(rowt[1])


        conn.close()



    def reset_fields(self):

        self.console.clear()


        for i in range(5,17):
            self.tableWidget.setItem(i , 0, QTableWidgetItem(""))

        self.popolacombo()
        self.combo_contaminant.setEnabled(True)
        self.multiplesubstance_control=False
        self.listasostanze=[]
        self.listaconc=[]




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


    def cosdir_azim(self,azim):
        az = math.radians(azim)
        cosa = math.sin(az)
        cosb = math.cos(az)
        return cosa,cosb

    def leach(self):

        # inizio recupero e controllo dati da form



        self.text_line_sw=str(self.tableWidget.item(5,0).text())
        self.text_line_soild=str(self.tableWidget.item(6,0).text())
        self.text_line_sourcethick=str(self.tableWidget.item(7,0).text())
        self.text_line_conc=str(self.tableWidget.item(8,0).text())
        self.text_line_acquifer_depth=str(self.tableWidget.item(9,0).text())
        self.text_line_darcy=str(self.tableWidget.item(10,0).text())
        self.text_line_deltagw=str(self.tableWidget.item(11,0).text())
        # self.text_line_alfax=str(self.tableWidget.item(12,0).text())
        # self.text_line_alfay=str(self.tableWidget.item(13,0).text())
        self.text_line_1dec=str(self.tableWidget.item(12,0).text())
        self.text_line_dirwind=str(self.tableWidget.item(13,0).text())
        self.text_line_p=str(self.tableWidget.item(14,0).text())



        self.text_vector = str(self.combo_source.currentText())
        self.text_area = str(self.combo_bound.currentText())
        self.text_contaminant = str(self.combo_contaminant.currentText())
        self.text_texture = str(self.combo_texture.currentText())

        self.algoritmo = str(self.combo_alg.currentText())

        ########## in caso di scelta tra l'algoritmo fickian o domenico #######
        ########## va modificata anche l'interfaccia ui rinomincando ui_daf_double in ui_daf #######
        ########## commentata la prima riga qui di seguito e decommentate le altre #######
        ########## vedi anche riga 98 ########
        self.opzione="continuo"
        # self.opzione = str(self.combo_opz.currentText())

        #### fine

        self.res=int(self.spinRes.text())

        self.time=1



        ####controllo campi input


        if self.text_line_soild!='':
            try:
                self.ro=float(self.text_line_soild)
            except Exception as e:
                QMessageBox.warning(self,"Warning", "Errore nella densità" )
                return
        else:
            self.ro=1.70

        if self.text_line_sourcethick!='':
            try:
                self.dz=float(self.text_line_sourcethick)
            except Exception as e:
                QMessageBox.warning(self,"Warning", u"Errore nello spessore della sorgente" )
                return
        else:
            self.dz=1

        if self.text_line_darcy!='':
            try:
                self.ve=float(self.text_line_darcy)
            except Exception as e:
                QMessageBox.warning(self,"Warning", u"Errore nella velictà di Darcy" )
                return
        else:
            self.ve=0.000025

        if self.text_line_deltagw!='':
            try:
                self.dgw=float(self.text_line_deltagw)
            except Exception as e:
                QMessageBox.warning(self,"Warning", u"Errore nella profondità della mixed zone" )
                return
        else:
            self.dgw=1


        if self.text_line_sw!='':
            try:
                self.sw=float(self.text_line_sw)*1000
            except Exception as e:
                QMessageBox.warning(self,"Warning", "Errore nell'estensione ortogonale al flusso di falda della sorgente" )
                return
        else:
            self.sw=10000.0


        if self.text_line_1dec!='':
            try:
                self.lambda1=float(self.text_line_1dec)
            except Exception as e:
                QMessageBox.warning(self,"Warning", u"Errore nel coefficiente di decadimento di primo ordine" )
                return
        else:
            self.lambda1=0.0


        try:
            self.lf=float(self.text_line_acquifer_depth)
        except Exception as e:
            QMessageBox.warning(self,"Warning", u"La profondità dell'acquifero è obbligatoria" )
            return


        # attenzione che la concentrazione è in grammi e non in kg correggere anche sotto alla linea 379

        if self.multiplesubstance_control==False:
            try:
                self.C_kg=float(self.text_line_conc)
            except Exception as e:
                QMessageBox.warning(self,"Warning", u"La concentrazione iniziale è obbligatoria" )
                return


        try:
            self.pioggia=float(self.text_line_p)
        except Exception as e:
            QMessageBox.warning(self,"Warning", u"L'infiltrazione efficace è obbligatoria" )
            return

        try:
            self.dirfalda=int(self.text_line_dirwind)
        except Exception as e:
            QMessageBox.warning(self,"Warning", u"Indicare una direzione in gradi del flusso di falda" )
            return

        ####fine controllo campi input

        #pyqtRemoveInputHook()

        self.x_w=self.dirfalda-90

        self.azimut=int(self.dirfalda-90)

        # attenzione in realtà la massa è in grammi non in kg: correggere anche sopra alla linea 356
        #self.C=self.C_kg*1000

        self.source=self.listalayers[self.text_vector]

        # import pdb
        # pdb.set_trace()


        # wkbType: 1:point, 6:multipolygon, 2: Linestring

        if self.source.wkbType()!=1:
            QMessageBox.warning(self,"Warning", "La sorgente deve avere geometria puntuale" )
            return

        self.areastudio=self.listalayers[self.text_area]


        if self.areastudio.wkbType()!=6:
            QMessageBox.warning(self,"Warning", u"L'area di studio deve avere geometria poligonale" )
            return


        if self.areastudio.crs().authid()!=self.source.crs().authid():
            QMessageBox.warning(self,"Warning", "Errore: i sistemi di riferimento non sono uniformi. Impossibile continuare con l'analisi." )
            return

        self.refsys=self.source.crs().authid().split(':')[1]


        if self.multiplesubstance_control==False:
            self.listasostanze.append(self.inquinanti[self.text_contaminant])
            self.listaconc.append(self.C_kg)


        contatore_sostanza=0
        for sostanza in self.listasostanze:


            self.path_output=self.line_output.text()
            if self.path_output=="":
                self.path_output=os.path.dirname(__file__)+"/output_model"+str(sostanza)+".tif"

            #recupero dati database

            self.id_substance=sostanza
            self.id_texture=self.texture[self.text_texture]

            lst_fields=['c_henry','koc_kd']
            res_fields=functions.substance_extract(self.id_substance,lst_fields,os.path.dirname(__file__)+"/../library/")

            self.h=res_fields[0]
            self.kd=res_fields[1]

            lst_fields_t=['tot_por','c_water_avg','ief','por_eff','grain']
            res_texture=functions.texture_extract(self.text_texture,lst_fields_t,os.path.dirname(__file__)+"/../library/")

            self.tera_a=res_texture[0]
            self.tera_w=res_texture[1]
            #self.ief=res_texture[2]*self.p
            self.ief=res_texture[2]*math.pow(float(self.pioggia)/10,2)
            self.tera_e=res_texture[3]
            self.grain=res_texture[4]


            messaggio="Inizio elaborazione Acquifero Envifate\n"
            messaggio+="---------------------------\n\n"
            messaggio+="FILE DI INPUT:\n"
            messaggio+="Vettoriale sorgente: "+str(self.text_vector)+"\n"
            messaggio+="Vettoriale confine: "+str(self.text_area)+"\n"

            messaggio+="VARIABILI:\n"
            messaggio+="Sostanza inquinante: "+str(sostanza)+"\n"
            messaggio+=u"Concentrazione inquinante: "+str(self.listaconc[contatore_sostanza])+"\n"
            messaggio+="Tessitura del suolo: "+str(self.text_texture)+"\n\n"
            messaggio+=u"Densità del suolo: "+str(self.ro)+"\n"
            messaggio+="Spessore sorgente: "+str(self.dz)+"\n"
            messaggio+=u"Velocità Darcy: "+str(self.ve)+"\n"
            messaggio+=u"Profondità mixed zone: "+str(self.dgw)+"\n"
            messaggio+=u"Estensione orizzontale: "+str(self.sw)+"\n"
            messaggio+="Coefficiente di decadimento: "+str(self.lambda1)+"\n"
            messaggio+=u"Profondità acquifero: "+str(self.lf)+"\n"

            messaggio+=u"Precipitazione media annua: "+str(self.pioggia)+"\n"
            messaggio+=u"Direzione corrente falda: "+str(self.dirfalda)+"\n"
            messaggio+="Risoluzione: "+str(self.res)+"\n\n"
            messaggio+="Nome mappa prodotta: "+str(os.path.basename(self.path_output+str(sostanza)))+"\n\n"
            messaggio+="ALGORITMI UTILIZZATI: \n"
            # pyqtRemoveInputHook()
            # pdb.set_trace()
            messaggio+='Leaching: Calcolo del fattore di Lisciviazione (Agenzia per la Protezione dell’Ambiente. "Criteri metodologici per l\'applicazione dell\'analisi assoluta di rischio ai siti contaminati." (2008).\n\n'
            messaggio+='DAF: '+str(self.algoritmo)+' (Domenico P.A. e Schwartz F.W. (1998), Physical and Chemical Hydrogeology, John Wiley and Sons, New York)\n\n'
            messaggio+="---------------------------\n\n"


            self.console.appendPlainText(messaggio)

            #ve2=float(self.ve)*100*86400*365

            self.label_status.setText("Preparazione dati")
            self.label_status.setStyleSheet('color : #e8b445;font-weight:bold')


            #controllo se i valori sono presenti
            if '' in res_fields or '' in res_texture:
                messaggio="Analisi non effettuata - modificare i parametri"
                self.console.appendPlainText(messaggio)
                return

            element=leaching.leach(self.h,self.tera_w,self.tera_a,self.kd,self.ief,self.ro,self.dz,
                          self.lf,self.ve,self.dgw,self.sw)

            self.kw=element.calc_kw()
            self.ldf=element.calc_ldf()
            self.sam=element.calc_sam()

            self.LF=element.calc_LF()

            risultato="RISULTATI INTERMEDI:\n\n"
            risultato+="kw = %f \n" % self.kw

            risultato+="ldf = %f \n" % self.ldf

            #risultato+="sam = %f \n" % self.sam

            self.c0=self.listaconc[contatore_sostanza]*self.LF

            risultato+="fattore di lisciviazione : %f \n" % self.LF
            risultato+="concentrazione sorgente secondaria : %f mg/l \n" % self.c0

            self.console.appendPlainText(risultato)



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
                                   'modulo':'Dispersione in falda',
                                   'descrizione':'Simulazione di dispersione inquinante in falda',
                                   'srs':self.source.crs().authid(),
                                   'data':datetime.datetime.now().strftime("%d-%m-%y")})

            # if self.srid!="":
            #     srs = osr.SpatialReference()
            #     srs.ImportFromEPSG(self.srid)
            #     target_ds.SetProjection( srs.ExportToWkt() )
            # else:
            #     target_ds.SetProjection(projectionfrom)

            band = target_ds.GetRasterBand(1)
            band.SetNoDataValue(float(NoData_value))
            band.Fill(NoData_value)

            xsize = band.XSize
            ysize = band.YSize

            outData = np.array(band.ReadAsArray(0, 0, xsize,ysize).astype(np.float))

            polygons = [feature for feature in self.areastudio.getFeatures()]

            feature = next(self.source.getFeatures())
            geom = feature.geometry().asPoint()
            x_source=geom[0]
            y_source=geom[1]


            rows=ysize-1
            cols=xsize-1

            intervallo=self.res


            ########## da eliminare ##########
            self.max=1000
            ##################################

            max_progress=self.max/intervallo
            self.progressBar.setMaximum(max_progress)
            start_time = time.time()


            index_progress=0
            controllo=1



            length = intervallo

            # vl = QgsVectorLayer("Point", "Concentration", "memory")
            # pr = vl.dataProvider()
            # prfield=pr.addAttributes( [ QgsField("distance", QVariant.Int) ] )
            # prfield1=pr.addAttributes( [ QgsField("conc", QVariant.Double) ] )
            # vl.updateFields()

            original_point = QgsPoint(x_source,y_source)

            # feats = []

            self.list_result=[]

            self.sw=self.sw/1000




            while length <= self.max:
                index_progress+=1
                self.progressBar.setValue(index_progress)
                cosa, cosb = self.cosdir_azim(self.azimut)
                end_point = QgsPoint(original_point.x()+(length*cosa), original_point.y()+(length*cosb))


                xanalisi=length*100

                calcolo_daf=daf.class_daf(self.c0,xanalisi,0,0,0,
                                self.lambda1,self.ve,self.kd,self.ro,self.tera_e,1,self.time)
                daf_p=calcolo_daf.calc_DAF_ispra()
                cfinal_p=self.c0*daf_p
                #pdb.set_trace()
                # if cfinal_p>0:
                #     pdb.set_trace()
                # fet = QgsFeature()
                # fet.initAttributes(2)
                # fet.setAttribute(0,length)
                # fet.setAttribute(1,cfinal_p)
                # vl.updateFeature(fet)
                # fet.setGeometry(QgsGeometry.fromPoint(end_point))
                # feats.append(fet)

                self.list_result.append(cfinal_p)

                length = length+ intervallo

            # pr.addFeatures(feats)
            # vl.updateFields()
            # QgsMapLayerRegistry.instance().addMapLayer(vl)

            self.progressBar.setValue(0)
            max_progress=rows*cols
            self.progressBar.setMaximum(max_progress)
            start_time = time.time()


            self.label_status.setText("Processing data")
            self.label_status.setStyleSheet('color : #e8b445;font-weight:bold')

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
                            # pyqtRemoveInputHook()
                            # pdb.set_trace()
                            if poly.contains(punto_controllo):

                                #z= self.dem.dataProvider().identify(QgsPoint(x, y),QgsRaster.IdentifyFormatValue)
                                deltax=x-x_source
                                deltay=y-y_source
                                dist=math.sqrt(math.pow(deltay,2)+math.pow(deltax,2))

                                xvero=deltax*math.cos(math.radians(self.azimut))-deltay*math.sin(math.radians(self.azimut))
                                yvero=deltax*math.sin(math.radians(self.azimut))+deltay*math.cos(math.radians(self.azimut))

                                if xvero>0:

                                    element_daf=daf.class_daf(self.c0,xvero,yvero,0,0,self.lambda1,self.ve,self.kd,self.ro,self.tera_e,self.sw,self.time)
                                    if self.algoritmo=="Fickiano":

                                        if self.opzione=="impulso":
                                            cfinal=element_daf.calc_DAF()
                                        elif self.opzione=="continuo":
                                            cfinal=element_daf.calc_DAF_c()
                                    else:
                                        cfinal=self.c0*element_daf.calc_DAF_ispra()

                                    outData[row,col]=cfinal
                                else:
                                    outData[row,col]=0

                self.label_status.setText("Preparazione output")
                self.label_status.setStyleSheet('color : #e8b445;font-weight:bold')

                outData_raster=outData[::-1]
                band.WriteArray(outData_raster)
                astats=band.GetStatistics(0, 1)


            contatore_sostanza+=1

            band= None
            target_ds = None

            base_raster_name=os.path.basename(self.path_output+str(sostanza))
            raster_name=os.path.splitext(base_raster_name)[0]
            self.iface.addRasterLayer(self.path_output, raster_name)

            contatore_sostanza=0
            self.list_result=[]

            tempoanalisi=time.time() - start_time
            tempostimato=time.strftime("%H:%M:%S", time.gmtime(tempoanalisi))
            messaggio="---------------------------------\n"
            messaggio+="Fine modellazione\n"
            messaggio+="\nTempo di analisi: "+tempostimato+"\n"
            messaggio+="---------------------------------\n\n"
            messaggio+="ANALISI STATSTICHE DI BASE\nvalore minimo: "+str(astats[0])+"\n"+"valore massimo: "+str(astats[1])+"\n"+"valore medio: "+str(astats[2])+"\n"+"deviazione standard: "+str(astats[3])
            self.console.appendPlainText(messaggio)

        if self.multiplesubstance_control==False:
            ax1f1 = self.figure.add_subplot(111)
            ax1f1.plot(self.list_result)

            self.canvas_mat.draw()

            self.label_status.setText("In attesa di dati")
            self.label_status.setStyleSheet('color : green; font-weight:bold')
