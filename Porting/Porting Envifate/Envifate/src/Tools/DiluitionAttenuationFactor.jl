module DiluitionAttenuationfactor

# -*- coding: utf-8 -*-
#=
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
=#

#= imports
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
=#

using ArchGDAL
using Dates

include("../Library/Daf.jl")
include("../Library/Leaching.jl")

#=
    mutable struct Dialog
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
=#


#=
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
    end

    function esporta_output(ef::Dialog)
        resultmodel = ef.console.toPlainText()
        name = QFileDialog.getSaveFileName(ef, 'Save File')
        file = open(name[0],'w')
        file.write(resultmodel)
        file.close()

    function reset_output(ef::Dialog)
        ret = QMessageBox.warning(ef,"Attenzione", "Vuoi davvero eliminare i risultati del modello?",QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        if ret== QMessageBox.Yes:
            ef.console.clear()
        else:
            return False

    function run1(ef::Dialog)
        # fix_print_with_import
        print("run effettuato")

    function set_combo_opz(ef::Dialog)
        check_alg = str(ef.combo_alg.currentText())
        if check_alg=="Domenico/Schwartz":
            ef.combo_opz.setEnabled(False)
        else:
            ef.combo_opz.setEnabled(True)

    function help(ef::Dialog)
        #ef.credits = u"Università della Tuscia\n Viterbo - Italy\nRaffaele Pelorosso, Federica Gobattoni\nDeveloper: Francesco Geri"
        #QMessageBox.about(ef.dlg,"Credits", ef.credits )
        if platform.uname()[0]=="Windows":
            os.system("start "+os.path.dirname(__file__)+"/../tutorial/modello_dispersione_in_falda.pdf")
        if platform.uname()[0]=="Linux":
            os.system("xdg-open "+os.path.dirname(__file__)+"/../tutorial/modello_dispersione_in_falda.pdf")
        else:
            os.system("open "+os.path.dirname(__file__)+"/../tutorial/modello_dispersione_in_falda.pdf")

    function configuration(ef::Dialog)
        ef.d = do_setting.Dialog(ef.iface)

        ef.listView=QListView(ef)
        ef.d.horizontalLayout.addWidget(ef.listView)
        ef.modelConf = QStandardItemModel()

        ef.listView.setModel(ef.modelConf)
        servizi=['https://idt2.regione.veneto.it/geoportal/csw','https://idt2-geoserver.regione.veneto.it/geoserver/ows','https://idt2.regione.veneto.it/gwc/service/wmts'];

        for i in servizi:
            item = QStandardItem(i)
            item.setCheckable(True)
            ef.modelConf.appendRow(item)

        ef.AddButton = QPushButton('Add', ef)

        ef.d.VerticalLayout.addWidget(ef.AddButton)

        ef.d.show()
        ef.d.exec_()


    function matrice_iter_daf(ef::Dialog)
        ef.dmatrix = do_setting.Dialog(ef.iface)

        ef.tablesubstance = QTableWidget()

        ef.tablesubstance.setRowCount(1)
        ef.tablesubstance.setColumnCount(2)

        ef.dmatrix.horizontalLayout.addWidget(ef.tablesubstance)

        ef.AddSubstance = QPushButton('+', ef)
        ef.RemSubstance = QPushButton('-', ef)
        ef.ReadSubstance = QPushButton('Leggi', ef)

        ef.AddSubstance.clicked.connect(ef.aggiungiriga)
        ef.RemSubstance.clicked.connect(ef.togliriga)
        ef.ReadSubstance.clicked.connect(ef.leggisostanza)




        ef.dmatrix.VerticalLayout.addWidget(ef.AddSubstance)
        ef.dmatrix.VerticalLayout.addWidget(ef.RemSubstance)
        ef.dmatrix.VerticalLayout.addWidget(ef.ReadSubstance)
        combo_inq_matrix = QComboBox()
        for row in ef.sql_fetch_inq:
            combo_inq_matrix.addItem(row[1])

        #ef.tablesubstance.insertRow(ef.tablesubstance.rowCount())

        ef.tablesubstance.setCellWidget(ef.countmatrix_substance,0, combo_inq_matrix)
        ef.countmatrix_substance+=1


        ef.tablesubstance.setHorizontalHeaderLabels((u'Sostanza', u'Concentrazione'))

        ef.tablesubstance.horizontalHeader().setStretchLastSection(True)

        ef.combo_contaminant.setEnabled(False)
        ef.multiplesubstance_control=True

        ef.tableWidget.item(8,0).setFlags( Qt.ItemIsSelectable |  Qt.ItemIsEnabled )
        ef.tableWidget.item(8,0).setBackground(QtGui.QColor(166,166,166))

        # pyqtRemoveInputHook()
        # pdb.set_trace()

        ef.dmatrix.show()
        ef.dmatrix.exec_()



    function aggiungiriga(ef::Dialog)
        combo_inq_matrix = QComboBox()
        for row in ef.sql_fetch_inq:
            combo_inq_matrix.addItem(row[1])
        ef.tablesubstance.insertRow(ef.tablesubstance.rowCount())
        ef.tablesubstance.setCellWidget(ef.countmatrix_substance,0, combo_inq_matrix)
        ef.countmatrix_substance+=1
        # pyqtRemoveInputHook()
        # pdb.set_trace()


    function togliriga(ef::Dialog)
        if ef.countmatrix_substance>1:
            ef.tablesubstance.removeRow(1)
            ef.countmatrix_substance-=1


    function leggisostanza(ef::Dialog)

        ef.listasostanze=[]
        for nsostanza in range(0,ef.tablesubstance.rowCount()):
            #ef.diz_perm[ef.table.item(ef.c,0).text()]=ef.table.item(ef.c,1).text()
            try:
                id_substance=ef.inquinanti[str(ef.tablesubstance.cellWidget(nsostanza,0).currentText())]
                ef.listasostanze.append(id_substance)
                ef.listaconc.append(float(ef.tablesubstance.item(nsostanza,1).text())*1000)
            except:
                QMessageBox.warning(ef,"Warning", "Errore nella lista delle concentrazioni" )
                return



        ef.dmatrix.close()


    function about(ef::Dialog)
        QMessageBox.about(ef, "Credits EnviFate",u"""<p>EnviFate: Open source tool for environmental risk analysis<br />Release 1.0<br />13-1-2017<br />License: GPL v. 3<br /><a href='https://bitbucket.org/fragit/envifate'>Home page plugin</a></p><hr><p>Lavoro svolto nell’ambito del  Progetto  di   ricerca   scientifica  “Definizione  di   metodi   standard     e  di strumenti applicativi   informatici per   il calcolo degli effetti dei fattori di perturbazione   ai sensi della decisione  2011/484/Ue,  da impiegarsi  nell’ambito  della valutazione di incidenza” finanziato dalla Regione Veneto. Partner principale è il DICAM, Dipartimento di Ingegneria Civile Ambientale e Meccanica dell’Università di Trento (Italia).</p><hr><p>Autori: Francesco Geri, Marco Ciolli</p><p>Universita' di Trento, Trento - Dipartimento di Ingegneria Civile Ambientale e Meccanica (DICAM) <a href="http://www.dicam.unitn.it/">www.dicam.unitn.it/</a></p><hr><p>Consulenti: Paolo Zatelli, Oscar Cainelli</p>""")


    function popolafields(ef,combo_in,combo_out):
        vect_source_text=combo_in.currentText()
        if vect_source_text!="":
            #vfields = ef.allLayers[mainvect].pendingFields()
            mainvect = QgsProject.instance().mapLayersByName( vect_source_text )[0]
            vfields = mainvect.pendingFields()
            combo_out.addItem("No field")
            for field in vfields:
                combo_out.addItem(field.name())

    function popolacombo(ef::Dialog)
        ef.progressBar.setValue(0)
        ef.combo_source.clear()
        ef.combo_bound.clear()
        # ef.combo_source_field.clear()
        ef.combo_contaminant.clear()
        ef.combo_texture.clear()
        ef.line_output.clear()

        ef.allLayers = ef.canvas.layers()
        ef.listalayers=dict()
        #combo_bound
        #elementovuoto="No required"
        for i in ef.allLayers:
            if i.type() == QgsMapLayer.VectorLayer:
                ef.listalayers[i.name()]=i
                ef.combo_source.addItem(str(i.name()))
                ef.combo_bound.addItem(str(i.name()))

        # ef.popolafields(ef.combo_source,ef.combo_source_field)
        # pyqtRemoveInputHook()
        # pdb.set_trace()
        conn = sqlite3.connect(os.path.dirname(__file__)+"/../library/substance.db")
        cursor=conn.cursor()
        query_substance="select id,nome from substance"
        cursor.execute(query_substance)
        ef.sql_fetch_inq=cursor.fetchall()

        ef.inquinanti=dict()
        for row in ef.sql_fetch_inq:
            ef.inquinanti[row[1]]=row[0]
            ef.combo_contaminant.addItem(row[1])

        query_texture="select id,nome from texture"
        cursor.execute(query_texture)
        sql_texture=cursor.fetchall()
        #pyqtRemoveInputHook()
        #pdb.set_trace()
        ef.texture=dict()
        for rowt in sql_texture:
            ef.texture[rowt[1]]=rowt[0]
            ef.combo_texture.addItem(rowt[1])


        conn.close()



    function reset_fields(ef::Dialog)

        ef.console.clear()


        for i in range(5,17):
            ef.tableWidget.setItem(i , 0, QTableWidgetItem(""))

        ef.popolacombo()
        ef.combo_contaminant.setEnabled(True)
        ef.multiplesubstance_control=False
        ef.listasostanze=[]
        ef.listaconc=[]




    function scegli_file(ef,tipofile)
        if tipofile == "sqlite"
            ef.fname = QFileDialog.getOpenFileName(None, 'Open file', '/home','sqlite3 files (*.sqlite);;all files (*.*)')
            ef.dlg_conf.pathtodb.setText(ef.fname)
        if tipofile=="csv":
            ef.fname = QFileDialog.getOpenFileName(None, 'Open file', '/home','csv files (*.csv);;all files (*.*)')
            ef.dlg_conf.path_to_kmean.setText(ef.fname)
        if tipofile=="tif":
            ef.fname = QFileDialog.getOpenFileName(None, 'Open file', '/home','GeoTiff files (*.tif);;all files (*.*)')
            ef.dlg_reclass.output_raster_class.setText(ef.fname)
        if tipofile=="tutti":
            ef.fname = QFileDialog.getOpenFileName(None, 'Open file', '/home','all files (*.*)')
            ef.dlg_reclass.input_reclass.setText(ef.fname)
        if tipofile=="salvaraster":
            ef.fname = QFileDialog.getSaveFileName(None, 'Save file', '/home','GeoTiff files (*.tif);;all files (*.*)')
            ef.line_output.setText(ef.fname[0])
        # if ef.tipofile=="salvacsv":
        #     ef.fname = QFileDialog.getSaveFileName(None, 'Save file', '/home','csv files (*.csv);;all files (*.*)')
        #     ef.dlg.lineEdit_csv.setText(ef.fname)


    function cosdir_azim(ef,azim):
        az = math.radians(azim)
        cosa = math.sin(az)
        cosb = math.cos(az)
        return cosa,cosb
=#


#=
function leach(ef::Dialog)

    # inizio recupero e controllo dati da form



    ef.text_line_sw=str(ef.tableWidget.item(5,0).text())
    ef.text_line_soild=str(ef.tableWidget.item(6,0).text())
    ef.text_line_sourcethick=str(ef.tableWidget.item(7,0).text())
    ef.text_line_conc=str(ef.tableWidget.item(8,0).text())
    ef.text_line_acquifer_depth=str(ef.tableWidget.item(9,0).text())
    ef.text_line_darcy=str(ef.tableWidget.item(10,0).text())
    ef.text_line_deltagw=str(ef.tableWidget.item(11,0).text())
    # ef.text_line_alfax=str(ef.tableWidget.item(12,0).text())
    # ef.text_line_alfay=str(ef.tableWidget.item(13,0).text())
    ef.text_line_1dec=str(ef.tableWidget.item(12,0).text())
    ef.text_line_dirwind=str(ef.tableWidget.item(13,0).text())
    ef.text_line_p=str(ef.tableWidget.item(14,0).text())



    ef.text_vector = str(ef.combo_source.currentText())
    ef.text_area = str(ef.combo_bound.currentText())
    ef.text_contaminant = str(ef.combo_contaminant.currentText())
    ef.text_texture = str(ef.combo_texture.currentText())

    ef.algoritmo = str(ef.combo_alg.currentText())

    ########## in caso di scelta tra l'algoritmo fickian o domenico #######
    ########## va modificata anche l'interfaccia ui rinomincando ui_daf_double in ui_daf #######
    ########## commentata la prima riga qui di seguito e decommentate le altre #######
    ########## vedi anche riga 98 ########
    ef.opzione="continuo"
    # ef.opzione = str(ef.combo_opz.currentText())

    #### fine

    ef.res=int(ef.spinRes.text())

    ef.time=1



    ####controllo campi input


    if ef.text_line_soild!='':
        try:
            ef.ro=float(ef.text_line_soild)
        except Exception as e:
            QMessageBox.warning(ef,"Warning", "Errore nella densità" )
            return
    else:
        ef.ro=1.70

    if ef.text_line_sourcethick!='':
        try:
            ef.dz=float(ef.text_line_sourcethick)
        except Exception as e:
            QMessageBox.warning(ef,"Warning", u"Errore nello spessore della sorgente" )
            return
    else:
        ef.dz=1

    if ef.text_line_darcy!='':
        try:
            ef.ve=float(ef.text_line_darcy)
        except Exception as e:
            QMessageBox.warning(ef,"Warning", u"Errore nella velictà di Darcy" )
            return
    else:
        ef.ve=0.000025

    if ef.text_line_deltagw!='':
        try:
            ef.dgw=float(ef.text_line_deltagw)
        except Exception as e:
            QMessageBox.warning(ef,"Warning", u"Errore nella profondità della mixed zone" )
            return
    else:
        ef.dgw=1


    if ef.text_line_sw!='':
        try:
            ef.sw=float(ef.text_line_sw)*1000
        except Exception as e:
            QMessageBox.warning(ef,"Warning", "Errore nell'estensione ortogonale al flusso di falda della sorgente" )
            return
    else:
        ef.sw=10000.0


    if ef.text_line_1dec!='':
        try:
            ef.lambda1=float(ef.text_line_1dec)
        except Exception as e:
            QMessageBox.warning(ef,"Warning", u"Errore nel coefficiente di decadimento di primo ordine" )
            return
    else:
        ef.lambda1=0.0


    try:
        ef.lf=float(ef.text_line_acquifer_depth)
    except Exception as e:
        QMessageBox.warning(ef,"Warning", u"La profondità dell'acquifero è obbligatoria" )
        return


    # attenzione che la concentrazione è in grammi e non in kg correggere anche sotto alla linea 379

    if ef.multiplesubstance_control==False:
        try:
            ef.C_kg=float(ef.text_line_conc)
        except Exception as e:
            QMessageBox.warning(ef,"Warning", u"La concentrazione iniziale è obbligatoria" )
            return


    try:
        ef.pioggia=float(ef.text_line_p)
    except Exception as e:
        QMessageBox.warning(ef,"Warning", u"L'infiltrazione efficace è obbligatoria" )
        return

    try:
        ef.dirfalda=int(ef.text_line_dirwind)
    except Exception as e:
        QMessageBox.warning(ef,"Warning", u"Indicare una direzione in gradi del flusso di falda" )
        return

    ####fine controllo campi input

    #pyqtRemoveInputHook()

    ef.x_w=ef.dirfalda-90

    ef.azimut=int(ef.dirfalda-90)

    # attenzione in realtà la massa è in grammi non in kg: correggere anche sopra alla linea 356
    #ef.C=ef.C_kg*1000

    ef.source=ef.listalayers[ef.text_vector]

    # import pdb
    # pdb.set_trace()


    # wkbType: 1:point, 6:multipolygon, 2: Linestring

    if ef.source.wkbType()!=1:
        QMessageBox.warning(ef,"Warning", "La sorgente deve avere geometria puntuale" )
        return

    ef.areastudio=ef.listalayers[ef.text_area]


    if ef.areastudio.wkbType()!=6:
        QMessageBox.warning(ef,"Warning", u"L'area di studio deve avere geometria poligonale" )
        return


    if ef.areastudio.crs().authid()!=ef.source.crs().authid():
        QMessageBox.warning(ef,"Warning", "Errore: i sistemi di riferimento non sono uniformi. Impossibile continuare con l'analisi." )
        return

    ef.refsys=ef.source.crs().authid().split(':')[1]


    if ef.multiplesubstance_control==False:
        ef.listasostanze.append(ef.inquinanti[ef.text_contaminant])
        ef.listaconc.append(ef.C_kg)


    contatore_sostanza=0
    for sostanza in ef.listasostanze:


        ef.path_output=ef.line_output.text()
        if ef.path_output=="":
            ef.path_output=os.path.dirname(__file__)+"/output_model"+str(sostanza)+".tif"

        #recupero dati database

        ef.id_substance=sostanza
        ef.id_texture=ef.texture[ef.text_texture]

        lst_fields=['c_henry','koc_kd']
        res_fields=functions.substance_extract(ef.id_substance,lst_fields,os.path.dirname(__file__)+"/../library/")

        ef.h=res_fields[0]
        ef.kd=res_fields[1]

        lst_fields_t=['tot_por','c_water_avg','ief','por_eff','grain']
        res_texture=functions.texture_extract(ef.text_texture,lst_fields_t,os.path.dirname(__file__)+"/../library/")

        ef.tera_a=res_texture[0]
        ef.tera_w=res_texture[1]
        #ef.ief=res_texture[2]*ef.p
        ef.ief=res_texture[2]*math.pow(float(ef.pioggia)/10,2)
        ef.tera_e=res_texture[3]
        ef.grain=res_texture[4]


        messaggio="Inizio elaborazione Acquifero Envifate\n"
        messaggio+="---------------------------\n\n"
        messaggio+="FILE DI INPUT:\n"
        messaggio+="Vettoriale sorgente: "+str(ef.text_vector)+"\n"
        messaggio+="Vettoriale confine: "+str(ef.text_area)+"\n"

        messaggio+="VARIABILI:\n"
        messaggio+="Sostanza inquinante: "+str(sostanza)+"\n"
        messaggio+=u"Concentrazione inquinante: "+str(ef.listaconc[contatore_sostanza])+"\n"
        messaggio+="Tessitura del suolo: "+str(ef.text_texture)+"\n\n"
        messaggio+=u"Densità del suolo: "+str(ef.ro)+"\n"
        messaggio+="Spessore sorgente: "+str(ef.dz)+"\n"
        messaggio+=u"Velocità Darcy: "+str(ef.ve)+"\n"
        messaggio+=u"Profondità mixed zone: "+str(ef.dgw)+"\n"
        messaggio+=u"Estensione orizzontale: "+str(ef.sw)+"\n"
        messaggio+="Coefficiente di decadimento: "+str(ef.lambda1)+"\n"
        messaggio+=u"Profondità acquifero: "+str(ef.lf)+"\n"

        messaggio+=u"Precipitazione media annua: "+str(ef.pioggia)+"\n"
        messaggio+=u"Direzione corrente falda: "+str(ef.dirfalda)+"\n"
        messaggio+="Risoluzione: "+str(ef.res)+"\n\n"
        messaggio+="Nome mappa prodotta: "+str(os.path.basename(ef.path_output+str(sostanza)))+"\n\n"
        messaggio+="ALGORITMI UTILIZZATI: \n"
        # pyqtRemoveInputHook()
        # pdb.set_trace()
        messaggio+='Leaching: Calcolo del fattore di Lisciviazione (Agenzia per la Protezione dell’Ambiente. "Criteri metodologici per l\'applicazione dell\'analisi assoluta di rischio ai siti contaminati." (2008).\n\n'
        messaggio+='DAF: '+str(ef.algoritmo)+' (Domenico P.A. e Schwartz F.W. (1998), Physical and Chemical Hydrogeology, John Wiley and Sons, New York)\n\n'
        messaggio+="---------------------------\n\n"


        ef.console.appendPlainText(messaggio)

        #ve2=float(ef.ve)*100*86400*365

        ef.label_status.setText("Preparazione dati")
        ef.label_status.setStyleSheet('color : #e8b445;font-weight:bold')


        #controllo se i valori sono presenti
        if '' in res_fields or '' in res_texture:
            messaggio="Analisi non effettuata - modificare i parametri"
            ef.console.appendPlainText(messaggio)
            return

        element=leaching.leach(ef.h,ef.tera_w,ef.tera_a,ef.kd,ef.ief,ef.ro,ef.dz,
                        ef.lf,ef.ve,ef.dgw,ef.sw)

        ef.kw=element.calc_kw()
        ef.ldf=element.calc_ldf()
        ef.sam=element.calc_sam()

        ef.LF=element.calc_LF()

        risultato="RISULTATI INTERMEDI:\n\n"
        risultato+="kw = %f \n" % ef.kw

        risultato+="ldf = %f \n" % ef.ldf

        #risultato+="sam = %f \n" % ef.sam

        ef.c0=ef.listaconc[contatore_sostanza]*ef.LF

        risultato+="fattore di lisciviazione : %f \n" % ef.LF
        risultato+="concentrazione sorgente secondaria : %f mg/l \n" % ef.c0

        ef.console.appendPlainText(risultato)



        path_layer=ef.areastudio.dataProvider().dataSourceUri()
        path=path_layer.split("|")
        source_ds = ogr.Open(path[0])
        area_layer = source_ds.GetLayer()
        x_min=int(area_layer.GetExtent()[0])
        y_min=int(area_layer.GetExtent()[2])
        x_max=int(area_layer.GetExtent()[1])
        y_max=int(area_layer.GetExtent()[3])

        drivermem = gdal.GetDriverByName('MEM')
        pixel_size = ef.res
        NoData_value = -9999



        # Create the destination data source
        x_res = (x_max - x_min) / pixel_size
        y_res = (y_max - y_min) / pixel_size

        target_ds = gdal.GetDriverByName('GTiff').Create(ef.path_output, int(x_res), int(y_res), 1, gdal.GDT_Float32)
        target_ds.SetGeoTransform((x_min, pixel_size, 0, y_max, 0, -pixel_size))
        projectionfrom = target_ds.GetProjection()

        srs = osr.SpatialReference()
        srs.ImportFromEPSG(int(ef.refsys))
        target_ds.SetProjection( srs.ExportToWkt() )

        target_ds.SetMetadata({'credits':'Envifate - Francesco Geri, Oscar Cainelli, Paolo Zatelli, Gianluca Salogni, Marco Ciolli - DICAM Università degli Studi di Trento - Regione Veneto',
                                'modulo':'Dispersione in falda',
                                'descrizione':'Simulazione di dispersione inquinante in falda',
                                'srs':ef.source.crs().authid(),
                                'data':datetime.datetime.now().strftime("%d-%m-%y")})

        # if ef.srid!="":
        #     srs = osr.SpatialReference()
        #     srs.ImportFromEPSG(ef.srid)
        #     target_ds.SetProjection( srs.ExportToWkt() )
        # else:
        #     target_ds.SetProjection(projectionfrom)

        band = target_ds.GetRasterBand(1)
        band.SetNoDataValue(float(NoData_value))
        band.Fill(NoData_value)

        xsize = band.XSize
        ysize = band.YSize

        outData = np.array(band.ReadAsArray(0, 0, xsize,ysize).astype(np.float))

        polygons = [feature for feature in ef.areastudio.getFeatures()]

        feature = next(ef.source.getFeatures())
        geom = feature.geometry().asPoint()
        x_source=geom[0]
        y_source=geom[1]


        rows=ysize-1
        cols=xsize-1

        intervallo=ef.res


        ########## da eliminare ##########
        ef.max=1000
        ##################################

        max_progress=ef.max/intervallo
        ef.progressBar.setMaximum(max_progress)
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

        ef.list_result=[]

        ef.sw=ef.sw/1000




        while length <= ef.max:
            index_progress+=1
            ef.progressBar.setValue(index_progress)
            cosa, cosb = ef.cosdir_azim(ef.azimut)
            end_point = QgsPoint(original_point.x()+(length*cosa), original_point.y()+(length*cosb))


            xanalisi=length*100

            calcolo_daf=daf.class_daf(ef.c0,xanalisi,0,0,0,
                            ef.lambda1,ef.ve,ef.kd,ef.ro,ef.tera_e,1,ef.time)
            daf_p=calcolo_daf.calc_DAF_ispra()
            cfinal_p=ef.c0*daf_p
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

            ef.list_result.append(cfinal_p)

            length = length+ intervallo

        # pr.addFeatures(feats)
        # vl.updateFields()
        # QgsMapLayerRegistry.instance().addMapLayer(vl)

        ef.progressBar.setValue(0)
        max_progress=rows*cols
        ef.progressBar.setMaximum(max_progress)
        start_time = time.time()


        ef.label_status.setText("Processing data")
        ef.label_status.setStyleSheet('color : #e8b445;font-weight:bold')

        if controllo==1:

            for row in range(rows):

                for col in range(cols):
                    index_progress+=1
                    ef.progressBar.setValue(index_progress)
                    x = col*pixel_size+x_min+(pixel_size/2)
                    y = row*pixel_size+y_min+(pixel_size/2)

                    punto_controllo = QgsPointXY(x,y)


                    for pol in polygons:
                        poly = pol.geometry()
                        # pyqtRemoveInputHook()
                        # pdb.set_trace()
                        if poly.contains(punto_controllo):

                            #z= ef.dem.dataProvider().identify(QgsPoint(x, y),QgsRaster.IdentifyFormatValue)
                            deltax=x-x_source
                            deltay=y-y_source
                            dist=math.sqrt(math.pow(deltay,2)+math.pow(deltax,2))

                            xvero=deltax*math.cos(math.radians(ef.azimut))-deltay*math.sin(math.radians(ef.azimut))
                            yvero=deltax*math.sin(math.radians(ef.azimut))+deltay*math.cos(math.radians(ef.azimut))

                            if xvero>0:

                                element_daf=daf.class_daf(ef.c0,xvero,yvero,0,0,ef.lambda1,ef.ve,ef.kd,ef.ro,ef.tera_e,ef.sw,ef.time)
                                if ef.algoritmo=="Fickiano":

                                    if ef.opzione=="impulso":
                                        cfinal=element_daf.calc_DAF()
                                    elif ef.opzione=="continuo":
                                        cfinal=element_daf.calc_DAF_c()
                                else:
                                    cfinal=ef.c0*element_daf.calc_DAF_ispra()

                                outData[row,col]=cfinal
                            else:
                                outData[row,col]=0

            ef.label_status.setText("Preparazione output")
            ef.label_status.setStyleSheet('color : #e8b445;font-weight:bold')

            outData_raster=outData[::-1]
            band.WriteArray(outData_raster)
            astats=band.GetStatistics(0, 1)


        contatore_sostanza+=1

        band= None
        target_ds = None

        base_raster_name=os.path.basename(ef.path_output+str(sostanza))
        raster_name=os.path.splitext(base_raster_name)[0]
        ef.iface.addRasterLayer(ef.path_output, raster_name)

        contatore_sostanza=0
        ef.list_result=[]

        tempoanalisi=time.time() - start_time
        tempostimato=time.strftime("%H:%M:%S", time.gmtime(tempoanalisi))
        messaggio="---------------------------------\n"
        messaggio+="Fine modellazione\n"
        messaggio+="\nTempo di analisi: "+tempostimato+"\n"
        messaggio+="---------------------------------\n\n"
        messaggio+="ANALISI STATSTICHE DI BASE\nvalore minimo: "+str(astats[0])+"\n"+"valore massimo: "+str(astats[1])+"\n"+"valore medio: "+str(astats[2])+"\n"+"deviazione standard: "+str(astats[3])
        ef.console.appendPlainText(messaggio)

    if ef.multiplesubstance_control==False:
        ax1f1 = ef.figure.add_subplot(111)
        ax1f1.plot(ef.list_result)

        ef.canvas_mat.draw()

        ef.label_status.setText("In attesa di dati")
        ef.label_status.setStyleSheet('color : green; font-weight:bold')
=#






function leach( sw::Real=10000.0, soild::Real=1.70, sourcethick::Int64=1, conc, aquifer, darcy::Real=0.000025, Δgw::Int64=1, λ1::Real=0.0, dirwind, p, source, area, contaminant, texture, algorithm )
    srid = [ 3003, 3004, 32632, 32633, 3857, 4326 ]
    classiwind = Dict(
        "N" => 0,
        "NE" => 45,
        "E" => 90,
        "SE" => 135,
        "S" => 180,
        "SW" => 225,
        "W" => 270,
        "NW" => 315
    )
    classiwind2 = Dict(
        0 => 180,
        45 => 225,
        90 => 270,
        135 => 315,
        180 => 0,
        225 => 45,
        270 => 90,
        315 => 135
    )
    # inizio recupero e controllo dati da form
    # ef.text_line_alfax=str(ef.tableWidget.item(12,0).text())
    # ef.text_line_alfay=str(ef.tableWidget.item(13,0).text())



    ########## in caso di scelta tra l'algoritmo fickian o domenico #######
    ########## va modificata anche l'interfaccia ui rinomincando ui_daf_double in ui_daf #######
    ########## commentata la prima riga qui di seguito e decommentate le altre #######
    ########## vedi anche riga 98 ########
    ef.opzione="continuo"
    # ef.opzione = str(ef.combo_opz.currentText())

    #### fine

    resolutiojn = tryparse( Int64, spinRes )
    time=1

    ####controllo campi input


    try:
        ef.lf=float(ef.text_line_acquifer_depth)
    except Exception as e:
        QMessageBox.warning(ef,"Warning", u"La profondità dell'acquifero è obbligatoria" )
        return


    # attenzione che la concentrazione è in grammi e non in kg correggere anche sotto alla linea 379

    if ef.multiplesubstance_control==False:
        try:
            ef.C_kg=float(ef.text_line_conc)
        except Exception as e:
            QMessageBox.warning(ef,"Warning", u"La concentrazione iniziale è obbligatoria" )
            return


    try:
        ef.pioggia=float(ef.text_line_p)
    except Exception as e:
        QMessageBox.warning(ef,"Warning", u"L'infiltrazione efficace è obbligatoria" )
        return

    try:
        ef.dirfalda=int(ef.text_line_dirwind)
    except Exception as e:
        QMessageBox.warning(ef,"Warning", u"Indicare una direzione in gradi del flusso di falda" )
        return

    ####fine controllo campi input

    #pyqtRemoveInputHook()

    ef.x_w=ef.dirfalda-90

    ef.azimut=int(ef.dirfalda-90)

    # attenzione in realtà la massa è in grammi non in kg: correggere anche sopra alla linea 356
    #ef.C=ef.C_kg*1000

    ef.source=ef.listalayers[ef.text_vector]

    # import pdb
    # pdb.set_trace()


    # wkbType: 1:point, 6:multipolygon, 2: Linestring

    if ef.source.wkbType()!=1:
        QMessageBox.warning(ef,"Warning", "La sorgente deve avere geometria puntuale" )
        return

    ef.areastudio=ef.listalayers[ef.text_area]


    if ef.areastudio.wkbType()!=6:
        QMessageBox.warning(ef,"Warning", u"L'area di studio deve avere geometria poligonale" )
        return


    if ef.areastudio.crs().authid()!=ef.source.crs().authid():
        QMessageBox.warning(ef,"Warning", "Errore: i sistemi di riferimento non sono uniformi. Impossibile continuare con l'analisi." )
        return

    ef.refsys=ef.source.crs().authid().split(':')[1]


    if ef.multiplesubstance_control==False:
        ef.listasostanze.append(ef.inquinanti[ef.text_contaminant])
        ef.listaconc.append(ef.C_kg)


    contatore_sostanza=0
    for sostanza in ef.listasostanze:


        ef.path_output=ef.line_output.text()
        if ef.path_output=="":
            ef.path_output=os.path.dirname(__file__)+"/output_model"+str(sostanza)+".tif"

        #recupero dati database

        ef.id_substance=sostanza
        ef.id_texture=ef.texture[ef.text_texture]

        lst_fields=['c_henry','koc_kd']
        res_fields=functions.substance_extract(ef.id_substance,lst_fields,os.path.dirname(__file__)+"/../library/")

        ef.h=res_fields[0]
        ef.kd=res_fields[1]

        lst_fields_t=['tot_por','c_water_avg','ief','por_eff','grain']
        res_texture=functions.texture_extract(ef.text_texture,lst_fields_t,os.path.dirname(__file__)+"/../library/")

        ef.tera_a=res_texture[0]
        ef.tera_w=res_texture[1]
        #ef.ief=res_texture[2]*ef.p
        ef.ief=res_texture[2]*math.pow(float(ef.pioggia)/10,2)
        ef.tera_e=res_texture[3]
        ef.grain=res_texture[4]


        messaggio="Inizio elaborazione Acquifero Envifate\n"
        messaggio+="---------------------------\n\n"
        messaggio+="FILE DI INPUT:\n"
        messaggio+="Vettoriale sorgente: "+str(ef.text_vector)+"\n"
        messaggio+="Vettoriale confine: "+str(ef.text_area)+"\n"

        messaggio+="VARIABILI:\n"
        messaggio+="Sostanza inquinante: "+str(sostanza)+"\n"
        messaggio+=u"Concentrazione inquinante: "+str(ef.listaconc[contatore_sostanza])+"\n"
        messaggio+="Tessitura del suolo: "+str(ef.text_texture)+"\n\n"
        messaggio+=u"Densità del suolo: "+str(ef.ro)+"\n"
        messaggio+="Spessore sorgente: "+str(ef.dz)+"\n"
        messaggio+=u"Velocità Darcy: "+str(ef.ve)+"\n"
        messaggio+=u"Profondità mixed zone: "+str(ef.dgw)+"\n"
        messaggio+=u"Estensione orizzontale: "+str(ef.sw)+"\n"
        messaggio+="Coefficiente di decadimento: "+str(ef.lambda1)+"\n"
        messaggio+=u"Profondità acquifero: "+str(ef.lf)+"\n"

        messaggio+=u"Precipitazione media annua: "+str(ef.pioggia)+"\n"
        messaggio+=u"Direzione corrente falda: "+str(ef.dirfalda)+"\n"
        messaggio+="Risoluzione: "+str(ef.res)+"\n\n"
        messaggio+="Nome mappa prodotta: "+str(os.path.basename(ef.path_output+str(sostanza)))+"\n\n"
        messaggio+="ALGORITMI UTILIZZATI: \n"
        # pyqtRemoveInputHook()
        # pdb.set_trace()
        messaggio+='Leaching: Calcolo del fattore di Lisciviazione (Agenzia per la Protezione dell’Ambiente. "Criteri metodologici per l\'applicazione dell\'analisi assoluta di rischio ai siti contaminati." (2008).\n\n'
        messaggio+='DAF: '+str(ef.algoritmo)+' (Domenico P.A. e Schwartz F.W. (1998), Physical and Chemical Hydrogeology, John Wiley and Sons, New York)\n\n'
        messaggio+="---------------------------\n\n"


        ef.console.appendPlainText(messaggio)

        #ve2=float(ef.ve)*100*86400*365

        ef.label_status.setText("Preparazione dati")
        ef.label_status.setStyleSheet('color : #e8b445;font-weight:bold')


        #controllo se i valori sono presenti
        if '' in res_fields or '' in res_texture:
            messaggio="Analisi non effettuata - modificare i parametri"
            ef.console.appendPlainText(messaggio)
            return

        element=leaching.leach(ef.h,ef.tera_w,ef.tera_a,ef.kd,ef.ief,ef.ro,ef.dz,
                        ef.lf,ef.ve,ef.dgw,ef.sw)

        ef.kw=element.calc_kw()
        ef.ldf=element.calc_ldf()
        ef.sam=element.calc_sam()

        ef.LF=element.calc_LF()

        risultato="RISULTATI INTERMEDI:\n\n"
        risultato+="kw = %f \n" % ef.kw

        risultato+="ldf = %f \n" % ef.ldf

        #risultato+="sam = %f \n" % ef.sam

        ef.c0=ef.listaconc[contatore_sostanza]*ef.LF

        risultato+="fattore di lisciviazione : %f \n" % ef.LF
        risultato+="concentrazione sorgente secondaria : %f mg/l \n" % ef.c0

        ef.console.appendPlainText(risultato)



        path_layer=ef.areastudio.dataProvider().dataSourceUri()
        path=path_layer.split("|")
        source_ds = ogr.Open(path[0])
        area_layer = source_ds.GetLayer()
        x_min=int(area_layer.GetExtent()[0])
        y_min=int(area_layer.GetExtent()[2])
        x_max=int(area_layer.GetExtent()[1])
        y_max=int(area_layer.GetExtent()[3])

        drivermem = gdal.GetDriverByName('MEM')
        pixel_size = ef.res
        NoData_value = -9999



        # Create the destination data source
        x_res = (x_max - x_min) / pixel_size
        y_res = (y_max - y_min) / pixel_size

        target_ds = gdal.GetDriverByName('GTiff').Create(ef.path_output, int(x_res), int(y_res), 1, gdal.GDT_Float32)
        target_ds.SetGeoTransform((x_min, pixel_size, 0, y_max, 0, -pixel_size))
        projectionfrom = target_ds.GetProjection()

        srs = osr.SpatialReference()
        srs.ImportFromEPSG(int(ef.refsys))
        target_ds.SetProjection( srs.ExportToWkt() )

        target_ds.SetMetadata({'credits':'Envifate - Francesco Geri, Oscar Cainelli, Paolo Zatelli, Gianluca Salogni, Marco Ciolli - DICAM Università degli Studi di Trento - Regione Veneto',
                                'modulo':'Dispersione in falda',
                                'descrizione':'Simulazione di dispersione inquinante in falda',
                                'srs':ef.source.crs().authid(),
                                'data':datetime.datetime.now().strftime("%d-%m-%y")})

        # if ef.srid!="":
        #     srs = osr.SpatialReference()
        #     srs.ImportFromEPSG(ef.srid)
        #     target_ds.SetProjection( srs.ExportToWkt() )
        # else:
        #     target_ds.SetProjection(projectionfrom)

        band = target_ds.GetRasterBand(1)
        band.SetNoDataValue(float(NoData_value))
        band.Fill(NoData_value)

        xsize = band.XSize
        ysize = band.YSize

        outData = np.array(band.ReadAsArray(0, 0, xsize,ysize).astype(np.float))

        polygons = [feature for feature in ef.areastudio.getFeatures()]

        feature = next(ef.source.getFeatures())
        geom = feature.geometry().asPoint()
        x_source=geom[0]
        y_source=geom[1]


        rows=ysize-1
        cols=xsize-1

        intervallo=ef.res


        ########## da eliminare ##########
        ef.max=1000
        ##################################

        max_progress=ef.max/intervallo
        ef.progressBar.setMaximum(max_progress)
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

        ef.list_result=[]

        ef.sw=ef.sw/1000




        while length <= ef.max:
            index_progress+=1
            ef.progressBar.setValue(index_progress)
            cosa, cosb = ef.cosdir_azim(ef.azimut)
            end_point = QgsPoint(original_point.x()+(length*cosa), original_point.y()+(length*cosb))


            xanalisi=length*100

            calcolo_daf=daf.class_daf(ef.c0,xanalisi,0,0,0,
                            ef.lambda1,ef.ve,ef.kd,ef.ro,ef.tera_e,1,ef.time)
            daf_p=calcolo_daf.calc_DAF_ispra()
            cfinal_p=ef.c0*daf_p
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

            ef.list_result.append(cfinal_p)

            length = length+ intervallo

        # pr.addFeatures(feats)
        # vl.updateFields()
        # QgsMapLayerRegistry.instance().addMapLayer(vl)

        ef.progressBar.setValue(0)
        max_progress=rows*cols
        ef.progressBar.setMaximum(max_progress)
        start_time = time.time()


        ef.label_status.setText("Processing data")
        ef.label_status.setStyleSheet('color : #e8b445;font-weight:bold')

        if controllo==1:

            for row in range(rows):

                for col in range(cols):
                    index_progress+=1
                    ef.progressBar.setValue(index_progress)
                    x = col*pixel_size+x_min+(pixel_size/2)
                    y = row*pixel_size+y_min+(pixel_size/2)

                    punto_controllo = QgsPointXY(x,y)


                    for pol in polygons:
                        poly = pol.geometry()
                        # pyqtRemoveInputHook()
                        # pdb.set_trace()
                        if poly.contains(punto_controllo):

                            #z= ef.dem.dataProvider().identify(QgsPoint(x, y),QgsRaster.IdentifyFormatValue)
                            deltax=x-x_source
                            deltay=y-y_source
                            dist=math.sqrt(math.pow(deltay,2)+math.pow(deltax,2))

                            xvero=deltax*math.cos(math.radians(ef.azimut))-deltay*math.sin(math.radians(ef.azimut))
                            yvero=deltax*math.sin(math.radians(ef.azimut))+deltay*math.cos(math.radians(ef.azimut))

                            if xvero>0:

                                element_daf=daf.class_daf(ef.c0,xvero,yvero,0,0,ef.lambda1,ef.ve,ef.kd,ef.ro,ef.tera_e,ef.sw,ef.time)
                                if ef.algoritmo=="Fickiano":

                                    if ef.opzione=="impulso":
                                        cfinal=element_daf.calc_DAF()
                                    elif ef.opzione=="continuo":
                                        cfinal=element_daf.calc_DAF_c()
                                else:
                                    cfinal=ef.c0*element_daf.calc_DAF_ispra()

                                outData[row,col]=cfinal
                            else:
                                outData[row,col]=0

            ef.label_status.setText("Preparazione output")
            ef.label_status.setStyleSheet('color : #e8b445;font-weight:bold')

            outData_raster=outData[::-1]
            band.WriteArray(outData_raster)
            astats=band.GetStatistics(0, 1)


        contatore_sostanza+=1

        band= None
        target_ds = None

        base_raster_name=os.path.basename(ef.path_output+str(sostanza))
        raster_name=os.path.splitext(base_raster_name)[0]
        ef.iface.addRasterLayer(ef.path_output, raster_name)

        contatore_sostanza=0
        ef.list_result=[]

        tempoanalisi=time.time() - start_time
        tempostimato=time.strftime("%H:%M:%S", time.gmtime(tempoanalisi))
        messaggio="---------------------------------\n"
        messaggio+="Fine modellazione\n"
        messaggio+="\nTempo di analisi: "+tempostimato+"\n"
        messaggio+="---------------------------------\n\n"
        messaggio+="ANALISI STATSTICHE DI BASE\nvalore minimo: "+str(astats[0])+"\n"+"valore massimo: "+str(astats[1])+"\n"+"valore medio: "+str(astats[2])+"\n"+"deviazione standard: "+str(astats[3])
        ef.console.appendPlainText(messaggio)

    if ef.multiplesubstance_control==False:
        ax1f1 = ef.figure.add_subplot(111)
        ax1f1.plot(ef.list_result)

        ef.canvas_mat.draw()

        ef.label_status.setText("In attesa di dati")
        ef.label_status.setStyleSheet('color : green; font-weight:bold')

end # module