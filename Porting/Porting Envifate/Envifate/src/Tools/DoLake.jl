
module DoLake

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

#= imports
    from __future__ import print_function

    from builtins import str
    from builtins import range
    from qgis.PyQt import QtCore, QtGui
    from PyQt5.QtCore import QSettings, QTranslator, QCoreApplication, Qt, QObject, pyqtSignal
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

    import pdb

    from envifate_dialog import EnviDialog

    from configuration_dialog import ConfigurationDialog

    sys.path.append( os.path.dirname(__file__)+"/../library" )

    import functions, lake, do_setting
=#


export Dialog, #struct
       esporta_output, reset_output, help, run1,
       popolafields, configuration, popolacombo, reset_fields,
       scegli_file, about, run_lake #functions


using ArchGDAL
using Sys

include("DoSettings.jl")
include("../Library/Lakes.jl")
include("../Library/Functions.jl")


mutable struct Dialog(EnviDialog)

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

        self.tab_2.setEnabled(false)

        self.label_title.setText("Analisi inquinamento luminoso")
        self.label_title.setStyleSheet('background-color : qlineargradient(spread:pad, x1:0, y1:0, x2:1, y2:0, stop:0 #76B888, stop:1 rgba(0, 0, 0, 0)); color : black')
        self.tableWidget.setRowCount(11)
        self.tableWidget.horizontalHeader().setStretchLastSection(true)
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


        self.tableWidget.setCellWidget(8,0, cellWidget)

        self.spintime=QSpinBox()
        self.spinRes=QSpinBox()

        self.tableWidget.setCellWidget(9,0, self.spintime)
        self.tableWidget.setCellWidget(10,0, self.spinRes)

        self.spinRes.setValue(25)
        self.spintime.setValue(10)


        self.tableWidget.resizeRowsToContents();

        # rowPosition = self.tableWidget.rowCount()
        # self.tableWidget.insertRow(rowPosition)
        # self.tableWidget.setItem(rowPosition , 0, QtGui.QTableWidgetItem("text1"))
        self.tableWidget.setVerticalHeaderLabels((u'Vettoriale sorgente*', u'Vettoriale confine*', u'Direzione corrente',
                                                  u'Massa inquinante* (g)',u'Velocità corrente* (m/sec)',u'Coeff. diffusione asse X',
                                                  u'Coeff. diffusione asse Y',u'Coeff. λ',u'Output file',
                                                  u'Tempo (h)',u'Risoluzione'))



        self.label_status.setText("In attesa di dati")
        self.label_status.setStyleSheet('color : green; font-weight:bold')

        self.clear_out_button.clicked.connect(self.reset_output)
        self.save_out_button.clicked.connect(self.esporta_output)

        # self.web = QWebView()
        # self.web.load(QUrl("https://grass.osgeo.org/grass70/manuals/addons/r.green.biomassfor.theoretical.html"))
        # self.web_layout.addWidget(self.web)

        self.popolacombo()

        self.saveButton.clicked.connect(λ: self.scegli_file("salvaraster"))
        self.reset_field_button.clicked.connect(self.reset_fields)
        self.buttonBox.accepted.connect(self.run_lake)
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

        self.classiwind2={}
        self.classiwind2[0]=180
        self.classiwind2[45]=225
        self.classiwind2[90]=270
        self.classiwind2[135]=315
        self.classiwind2[180]=0
        self.classiwind2[225]=45
        self.classiwind2[270]=90
        self.classiwind2[315]=135

        #pyqtRemoveInputHook()
        #pdb.set_trace()
        self.tabWidget.removeTab(1)

        self.list_srid=[3003,3004,32632,32633,3857,4326]
end




function esporta_output(dialog::Dialog)
    resultmodel = dialog.console.toPlainText()
    name = QFileDialog.getSaveFileName( dialog, "Save File" )
    file = open(name[0],'w')
    file.write(resultmodel)
    file.close()
end


function reset_output(dialog::Dialog)
    ret = QMessageBox.warning( dialog, "Attenzione", "Vuoi davvero eliminare i risultati del modello?", QMessageBox.Yes | QMessageBox.No, QMessageBox.No )
    if ret == QMessageBox.Yes
        dialog.console.clear()
    else
        return false
    end
end


function help()
    #self.credits = u"Università della Tuscia\n Viterbo - Italy\nRaffaele Pelorosso, Federica Gobattoni\nDeveloper: Francesco Geri"
    #QMessageBox.about(self.dlg,"Credits", self.credits )
    path = "$(@__DIR__)\\..\\tutorial\\manuale_envifate_laghi.pdf"
    cmd = Sys.iswindows() ? "start" : Sys.islinux() ? "xdg-open" : "open"
    run(`$cmd $path`)
end


function run1()
    # fix_print_with_import
    print("run effettuato")
end


function popolafields( dialog::Dialog, combo_in, combo_out )
    vect_source_text = combo_in.currentText()
    if !isempty(vect_source_text)
        #vfields = self.allLayers[mainvect].pendingFields()
        mainvect = dialog.registry.mapLayersByName( vect_source_text )[0]
        vfields = mainvect.pendingFields()
        combo_out.addItem("No field")
        for field in vfields
            combo_out.addItem(field.name())
        end
    end
end


function configuration(dialog::Dialog)
    d = do_setting.Dialog(dialog.iface)
    d.show()
    d.exec_()
end


function popolacombo(dialog::Dialog)
    dialog.combo_source.clear()
    dialog.combo_maindirwind.clear()
    dialog.combo_bound.clear()
    dialog.line_output.clear()
    dialog.progressBar.setValue(0)

    dialog.allLayers = dialog.canvas.layers()
    dialog.listalayers = Dict()
    #elementovuoto = "No required"
    for i in dialog.allLayers
        if i.type() == QgsMapLayer.VectorLayer
            dialog.listalayers[ i.name() ] = i
            dialog.combo_source.addItem( string( i.name() ) )
            dialog.combo_bound.addItem( string( i.name() ) )
        end
    end

    #self.popolafields(self.combo_source,self.combo_source_field)

    dialog.combo_maindirwind.addItem("N")
    dialog.combo_maindirwind.addItem("NE")
    dialog.combo_maindirwind.addItem("E")
    dialog.combo_maindirwind.addItem("SW")
    dialog.combo_maindirwind.addItem("S")
    dialog.combo_maindirwind.addItem("SW")
    dialog.combo_maindirwind.addItem("W")
    dialog.combo_maindirwind.addItem("NW")
end


function reset_fields(dialog::Dialog)
    dialog.console.clear()

    for i in range(3, 7, step=1 )
        dialog.tableWidget.setItem( i , 0, QTableWidgetItem("") )
    end

    dialog.popolacombo()

    ##### da eliminare #####
    # self.line_speed.setText("4")
    # self.line_flickianx.setText("1000")
    # self.line_flickiany.setText("1000")
    # self.line_conc.setText("2000")
    # self.line_lambda.setText("0")
    # self.spintime.setValue(20)
    ##### fine da eliminare #####
end


function scegli_file( dialog::Dialog, tipofile )
    if tipofile == "sqlite"
        dialog.fname = QFileDialog.getOpenFileName( nothing, 'Open file', '/home', 'sqlite3 files (*.sqlite);;all files (*.*)' )
        dialog.dlg_conf.pathtodb.setText(dialog.fname)
    end
    if tipofile == "csv"
        dialog.fname = QFileDialog.getOpenFileName( nothing, 'Open file', '/home', 'csv files (*.csv);;all files (*.*)' )
        dialog.dlg_conf.path_to_kmean.setText(dialog.fname)
    end
    if tipofile == "tif"
        dialog.fname = QFileDialog.getOpenFileName( nothing, 'Open file', '/home', 'GeoTiff files (*.tif);;all files (*.*)' )
        dialog.dlg_reclass.output_raster_class.setText(dialog.fname)
    end
    if tipofile == "tutti"
        dialog.fname = QFileDialog.getOpenFileName( nothing, 'Open file', '/home', 'all files (*.*)' )
        dialog.dlg_reclass.input_reclass.setText(dialog.fname)
    end
    if tipofile == "salvaraster"
        dialog.fname = QFileDialog.getSaveFileName( nothing, 'Save file', '/home', 'GeoTiff files (*.tif);;all files (*.*)' )
        dialog.line_output.setText(dialog.fname[0])
    end
    # if self.tipofile=="salvacsv":
    #     self.fname = QFileDialog.getSaveFileName(None, 'Save file', '/home','csv files (*.csv);;all files (*.*)')
    #     self.dlg.lineEdit_csv.setText(self.fname)
    #end
end


function about(dialog::Dialog)
    QMessageBox.about(dialog, "Credits EnviFate",u"""<p>EnviFate: Open source tool for environmental risk analysis<br />Release 1.0<br />13-1-2017<br />License: GPL v. 3<br /><a href='https://bitbucket.org/fragit/envifate'>Home page plugin</a></p><hr><p>Lavoro svolto nell’ambito del  Progetto  di   ricerca   scientifica  “Definizione  di   metodi   standard     e  di strumenti applicativi   informatici per   il calcolo degli effetti dei fattori di perturbazione   ai sensi della decisione  2011/484/Ue,  da impiegarsi  nell’ambito  della valutazione di incidenza” finanziato dalla Regione Veneto. Partner principale è il DICAM, Dipartimento di Ingegneria Civile Ambientale e Meccanica dell’Università di Trento (Italia).</p><hr><p>Autori: Francesco Geri, Marco Ciolli</p><p>Universita' di Trento, Trento - Dipartimento di Ingegneria Civile Ambientale e Meccanica (DICAM) <a href="http://www.dicam.unitn.it/">www.dicam.unitn.it/</a></p><hr><p>Consulenti: Paolo Zatelli, Oscar Cainelli</p>""")
end


function run_lake(dialog::Dialog)

    dialog.text_line_conc = string( dialog.tableWidget.item(3,0).text() )
    dialog.text_line_speed = string( dialog.tableWidget.item(4,0).text() )
    dialog.text_line_flickianx = string( dialog.tableWidget.item(5,0).text() )
    dialog.text_line_flickiany = string( dialog.tableWidget.item(6,0).text() )
    dialog.text_line_λ = string( dialog.tableWidget.item(7,0).text() )

    dialog.text_vector = string( dialog.combo_source.currentText() )
    dialog.text_area = string( dialog.combo_bound.currentText() )

    #self.x_w=self.classiwind[self.combo_maindirwind.currentText()]
    dialog.x_w = dialog.classiwind[ dialog.combo_maindirwind.currentText() ]

    dialog.res = parse( Int64, dialog.spinRes.text() )
    dialog.ore1 = parse( Int64, dialog.spintime.text() )

    dialog.ore = dialog.ore1 * 3600

    if !isempty( dialog.text_line_flickianx )
        try
            dialog.fickian_x = parse( Float64, dialog.text_line_flickianx )
        catch e
            QMessageBox.warning( dialog, "Warning", "Errore nel coefficiente di trasporto Flickian per l'asse X" )
            return
        end
    else
        dialog.fickian_x = 0.05
    end
    if !isempty( dialog.text_line_flickiany )
        try
            dialog.fickian_y = parse( Float64, dialog.text_line_flickiany )
        catch e
            QMessageBox.warning( dialog, "Warning", "Errore nel coefficiente di trasporto Flickian per l'asse Y" )
            return
        end
    else
        dialog.fickian_y = 0.05
    end
    if !isempty( dialog.text_line_λ )
        try
            dialog.λk = parse( Float64, dialog.text_line_λ )
        catch e
            QMessageBox.warning( dialog, "Warning", "Errore nel coefficiente di decadimento Lambda" )
            return
        end
    else
        dialog.λk = 0.0
    end
    try
        dialog.vmedia = parse( Float64, dialog.text_line_speed )
    catch e
        QMessageBox.warning( dialog, "Warning", "La velocità media della corrente è obbligatoria" )
        return
    end
    try
        dialog.C = parse( Float64, dialog.text_line_conc )
    catch e
        QMessageBox.warning( dialog, "Warning", "La massa inquinante è obbligatoria" )
        return
    end

    # wkbType: 1:point, 6:multipolygon, 2: Linestring

    dialog.source = dialog.listalayers[ dialog.text_vector ]

    if dialog.source.wkbType() != 1
        QMessageBox.warning( dialog, "Warning", "The source file must have point geometry" )
        return
    end
    dialog.areastudio = dialog.listalayers[ dialog.text_area ]

    if dialog.areastudio.wkbType() != 6
        QMessageBox.warning( dialog, "Warning", "The boundaries file must have polygon geometry" )
        return
    end
    dialog.path_output = dialog.line_output.text()

    if isempty(dialog.path_output)
        dialog.path_output = os.path.dirname(__file__)*"/output_model.tif"
    end

    if dialog.areastudio.crs().authid()!=dialog.source.crs().authid():
        QMessageBox.warning( dialog, "Warning", "Errore: i sistemi di riferimento non sono uniformi. Impossibile continuare con l'analisi." )
        return
    end
    dialog.refsys = dialog.source.crs().authid().split(":")[1]


    #recupero dati database
 
    messaggio = "Inizio elaborazione dispersione laghi e lagune Envifate\n" *
                    "---------------------------\n\n" *
                    "FILE DI INPUT:\n" *
                    "Vettoriale sorgente: "*string( dialog.text_vector ) * "\n" *
                    "Vettoriale confine: "*string( dialog.text_area ) * "\n" *
                    "VARIABILI:\n" *
                    "Direzione corrente: "*string( dialog.combo_source.currentText() )*"\n" *
                    "Coefficiente Fickian asse X: "*string( dialog.text_line_flickianx )*"\n" *
                    "Coefficiente Fickian asse Y: "*string( dialog.text_line_flickiany )*"\n" *
                    "Coefficiente decadimento lamda: "*tring( dialog.λk )*"\n" *
                    "Massa inquinante: "*string( dialog.C )*"\n" *
                    "Tempo dall'iniezione: "*string( dialog.spinRes.text() )*"\n" *
                    "Risoluzione: "*string( dialog.res )*"\n\n" *
                    "ALGORITMO UTILIZZATO: Fickian Mixing Process
                     (Hemond, Harold F., and Elizabeth J. Fechner. Chemical fate and transport in the environment. Elsevier, 2014.)\n\n" *
                    "---------------------------\n\n"
    dialog.console.appendPlainText(messaggio)

    dialog.label_status.setText("Preparazione dati")
    dialog.label_status.setStyleSheet("color : #e8b445;font-weight:bold")

    vx = round( dialog.vmedia * cos( deg2rad(dialog.x_w) ), digits=3 )
    vy = round( dialog.vmedia * sin( deg2rad(dialog.x_w) ), digits=3 )
    dialog.velocity_x = √(vx^2)
    dialog.velocity_y = √(vy^2)

    path_layer = dialog.areastudio.dataProvider().dataSourceUri()
    path = split( path_layer, "|" )
    source_ds = ogr.Open(path[0])
    area_layer = source_ds.GetLayer()
    x_min, y_min, x_max, y_max = round.( Int64, [ area_layer.GetExtent()[0], area_layer.GetExtent()[2], area_layer.GetExtent()[1], area_layer.GetExtent()[3] ] )

    drivermem = ArchGDAL.getdriver("MEM")
    pixel_size = dialog.res
    NoData_value = -9999

    # Create the destination data source
    x_res = (x_max - x_min) / pixel_size
    y_res = (y_max - y_min) / pixel_size

    target_ds = ArchGDAL.getdriver("GTiff").Create( dialog.path_output, round( Int64, x_res ), round( Int64, y_res ), 1, ArchGDAL.GDAL.GDT_Float32 )
    target_ds.SetGeoTransform( x_min, pixel_size, 0, y_max, 0, -pixel_size )
    projectionfrom = target_ds.GetProjection()

    srs = osr.SpatialReference()
    srs.ImportFromEPSG( parse( Int64, dialog.refsys ) )
    target_ds.SetProjection( srs.ExportToWkt() )
    # geotransform = target_ds.GetGeoTransform()
    target_ds.SetMetadata({ "credits" => "Envifate - Francesco Geri, Oscar Cainelli, Paolo Zatelli, Gianluca Salogni, Marco Ciolli - DICAM Università degli Studi di Trento - Regione Veneto",
                            "modulo" => "Dispersione in laghi e bacini",
                            "descrizione" => "Simulazione di dispersione inquinante all\'interno di un corpo idrico superficiale fermo",
                            "srs" => dialog.source.crs().authid(),
                            "data" => datetime.datetime.now().strftime("%d-%m-%y") } )


    band = target_ds.GetRasterBand(1)
    band.SetNoDataValue(convert( Float64, NoData_value))
    band.Fill(NoData_value)
    xsize = band.XSize
    ysize = band.YSize

    outData = np.array( band.ReadAsArray(0, 0, xsize,ysize).astype(np.float) )

    feature = next(dialog.source.getFeatures())
    geom = feature.geometry().asPoint()
    x_source = geom[0]
    y_source = geom[1]

    polygons = [ feature for feature in dialog.areastudio.getFeatures() ]

    rows = ysize - 1
    cols = xsize - 1

    max_progress = rows * cols
    dialog.progressBar.setMaximum(max_progress)
    start_time = time.time()

    dialog.label_status.setText("Processing data")
    dialog.label_status.setStyleSheet("color : #e8b445;font-weight:bold")

    index_progress = 0
    controllo = 1
    if controllo == 1
        for row in range( 1, length(rows), step=1 )
            for col in range( 1, cols, step=1 )
                index_progress += 1
                dialog.progressBar.setValue( index_progress )
                x = col * pixel_size + x_min + ( pixel_size / 2 )
                y = row * pixel_size + y_min + ( pixel_size / 2 )

                punto_controllo = QgsPointXY(x,y)

                # pyqtRemoveInputHook()
                # pdb.set_trace()
                for pol in polygons
                    poly = pol.geometry()
                    if poly.contains(punto_controllo)
                    # if geom_area.contains(punto_controllo):

                        Δx = x - x_source
                        Δy = y - y_source
                        xvero = Δx * cos( deg2rad(dialog.x_w) ) - Δy * sin( deg2rad(dialog.x_w) )
                        yvero = Δx * sin( deg2rad(dialog.x_w) ) + Δy * cos( deg2rad(dialog.x_w) )

                        element = lake.lake( dialog.C, dialog.ore, xvero, yvero, dialog.fickian_x, dialog.fickian_y,
                                             dialog.velocity_x, dialog.velocity_y, dialog.λk )

                        cfinal = element.calc_concentration()

                        outData[row,col] = cfinal
                    else
                        outData[row,col] = 0
                    end
                end
            end
        end

        dialog.label_status.setText("Preparazione output")
        dialog.label_status.setStyleSheet("color : #e8b445;font-weight:bold")

        outData_raster = outData[::-1]
        band.WriteArray(outData_raster)
        astats = band.GetStatistics(0, 1)
    end

    band = nothing
    target_ds = nothing

    base_raster_name = basename(dialog.path_output)
    raster_name = splitext(base_raster_name)[0]
    dialog.iface.addRasterLayer( dialog.path_output, raster_name )

    layer = nothing
    for lyr in collect( QgsProject.instance().mapLayers().values() )
        if lyr.name() == raster_name
            layer = lyr
        end
    end

    functions.applystyle( layer, "gr", 0.5 )

    tempoanalisi = time.time() - start_time
    tempostimato = time.strftime( "%H:%M:%S", time.gmtime(tempoanalisi) )
    messaggio = "---------------------------------\n" *
                    "Fine modellazione\n" *
                    "\nTempo di analisi: "*tempostimato*"\n" *
                    "---------------------------------\n\n" *
                    "ANALISI STATSTICHE DI BASE\nvalore minimo: "*string(round( Int64, astats[0], digits=4 ))*"\n" *
                    "valore massimo: "*string(round( Int64, astats[1], digits=4 ))*"\n" *
                    "valore medio: " * string(round( Int64, astats[2], digits=4 ))*"\n" *
                    "deviazione standard: "*string(round( Int64, astats[3], digits=4 ))
    dialog.console.appendPlainText(messaggio)

    dialog.label_status.setText("In attesa di dati")
    dialog.label_status.setStyleSheet('color : green; font-weight:bold')
end

end # module