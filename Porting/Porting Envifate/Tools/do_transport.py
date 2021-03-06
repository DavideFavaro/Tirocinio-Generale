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

#SET VARIABILI AMBIENTE
# Windows
from builtins import str
from builtins import range
grass7path = r'C:\OSGeo4W\apps\grass\grass-7.2.svn'
grass7bin_win = r'C:\OSGeo4W\bin\grass72svn.bat'
# Linux
#grass7bin_lin = 'grass72'
grass7bin_lin = 'grass'
# MacOSX
grass7bin_mac = '/Applications/GRASS/GRASS-7.1.app/'


from qgis.PyQt import QtCore, QtGui
from PyQt5.QtCore import QSettings, QTranslator, QCoreApplication, Qt, QObject, pyqtSignal 
from PyQt5.QtGui import QIcon 
from PyQt5.QtWidgets import QAction, QDialog, QFormLayout

import os,sys
import time

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

import subprocess
import shutil
import binascii
import tempfile

from osgeo import gdal,ogr


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

sys.path.append( os.path.dirname(__file__)+"/../library" )

import functions



if sys.platform.startswith('linux'):
    # we assume that the GRASS GIS start script is available and in the PATH
    # query GRASS 7 itself for its GISBASE
    grass7bin = grass7bin_lin
elif sys.platform.startswith('win'):
    grass7bin = grass7bin_win
else:
    OSError('Platform not configured.')

startcmd = grass7bin + ' --config path'

p = subprocess.Popen(startcmd, shell=True, 
                     stdout=subprocess.PIPE, stderr=subprocess.PIPE)
out, err = p.communicate()
if p.returncode != 0:
    # fix_print_with_import
    print('ERROR: %s' % err, file=sys.stderr)
    # fix_print_with_import
    print("ERROR: Cannot find GRASS GIS 7 start script (%s)" % startcmd, file=sys.stderr)
    sys.exit(-1)
if sys.platform.startswith('linux'):
    gisbase = out.strip('\n')
elif sys.platform.startswith('win'):
    if out.find("OSGEO4W home is") != -1:
        gisbase = out.strip().split('\n')[1]
    else:
        gisbase = out.strip('\n')
    os.environ['GRASS_SH'] = os.path.join(gisbase, 'msys', 'bin', 'sh.exe')

# Set GISBASE environment variable
os.environ['GISBASE'] = gisbase
# define GRASS-Python environment
gpydir = os.path.join(gisbase, "etc", "python")
sys.path.append(gpydir)
########
# define GRASS DATABASE
if sys.platform.startswith('win'):
    gisdb = os.path.join(os.getenv('APPDATA', 'grassdata'))
else:
    gisdb = os.path.join(os.getenv('HOME', 'grassdata'))

# override for now with TEMP dir
gisdb = os.path.join(tempfile.gettempdir(), 'grassdata')
try:
    os.stat(gisdb)
except:
    os.mkdir(gisdb)

# location/mapset: use random names for batch jobs
string_length = 16
location = binascii.hexlify(os.urandom(string_length))
mapset   = 'PERMANENT'
location_path = os.path.join(gisdb, location)

# Create new location (we assume that grass7bin is in the PATH)
#  from EPSG code:
#startcmd = grass7bin + ' -c epsg:' + str(self.myepsg) + ' -e ' + location_path
startcmd = grass7bin + ' -c epsg:4326 -e ' + location_path
#  from SHAPE or GeoTIFF file
#startcmd = grass7bin + ' -c ' + myfile + ' -e ' + location_path
# fix_print_with_import

# fix_print_with_import
print(startcmd)
p = subprocess.Popen(startcmd, shell=True, 
                     stdout=subprocess.PIPE, stderr=subprocess.PIPE)
out, err = p.communicate()
if p.returncode != 0:
    # fix_print_with_import
    print('ERROR: %s' % err, file=sys.stderr)
    # fix_print_with_import
    print('ERROR: Cannot generate location (%s)' % startcmd, file=sys.stderr)
    sys.exit(-1)
else:
    # fix_print_with_import
    print('Created location %s' % location_path)

# Now the location with PERMANENT mapset exists.

########
# Now we can use PyGRASS or GRASS Scripting library etc. after 
# having started the session with gsetup.init() etc

# Set GISDBASE environment variable
os.environ['GISDBASE'] = gisdb

# Linux: Set path to GRASS libs (TODO: NEEDED?)
path = os.getenv('LD_LIBRARY_PATH')
dir  = os.path.join(gisbase, 'lib')
if path:
    path = dir + os.pathsep + path
else:
    path = dir
os.environ['LD_LIBRARY_PATH'] = path

# language
os.environ['LANG'] = 'en_US'
os.environ['LOCALE'] = 'C'



import grass.script as grass
import grass.script.setup as gsetup
from grass.script.core import run_command, parser, overwrite, read_command, parse_command  


###########
# Launch session and do something
gsetup.init(gisbase, gisdb, location, mapset)   

class Dialog(EnviDialog):
    
    def __init__(self, iface):
        QDialog.__init__(self, iface.mainWindow())
        self.iface = iface
        self.canvas=self.iface.mapCanvas()
        self.registry = QgsMapLayerRegistry.instance()
        self.msgBar = self.iface.messageBar()
        # Set up the user interface from Designer.
        self.setupUi(self)

        self.tabWidget.setCurrentIndex(0)

        self.tabWidget.removeTab(2)

        self.label_title.setText("Simulazione di trasporto soluti in falda")
        self.label_title.setStyleSheet('background-color : qlineargradient(spread:pad, x1:0, y1:0, x2:1, y2:0, stop:0 #CC7D31, stop:1 rgba(0, 0, 0, 0)); color : black')
        self.tableWidget.setRowCount(25)
        self.tableWidget.horizontalHeader().setStretchLastSection(True)
        #self.tableWidget.horizontalHeaderItem(0).setText("newHeader")
        self.combo_bound = QtGui.QComboBox()
        self.combo_source = QtGui.QComboBox()
        self.combo_source.setObjectName("combo_source")
        self.combo_top = QtGui.QComboBox()
        self.combo_bottom = QtGui.QComboBox()
        self.combo_status = QtGui.QComboBox()
        self.combo_status.setObjectName("combo_status")
        self.combo_campostatus = QtGui.QComboBox()
        self.combo_campostatus.setObjectName("combo_campostatus")        
        self.combo_solver = QtGui.QComboBox()
        self.combo_conc = QtGui.QComboBox()
        self.combo_conc.setObjectName("combo_conc")
        self.combo_q = QtGui.QComboBox()
        self.combo_phead = QtGui.QComboBox()
        self.combo_phead.setObjectName("combo_phead")
        self.combo_campophead = QtGui.QComboBox()
        self.combo_campophead.setObjectName("combo_campophead")

        self.tableWidget.setCellWidget(0,0, self.combo_source)
        self.tableWidget.setCellWidget(1,0, self.combo_conc)
        self.tableWidget.setCellWidget(2,0, self.combo_bound)
        self.tableWidget.setCellWidget(3,0, self.combo_top)
        self.tableWidget.setItem(4 , 0, QtGui.QTableWidgetItem("20")) # value top

        self.tableWidget.setCellWidget(5,0, self.combo_bottom)
        self.tableWidget.setItem(6 , 0, QtGui.QTableWidgetItem("0")) # value bottom
        self.tableWidget.setCellWidget(7,0, self.combo_status)
        self.tableWidget.setCellWidget(8,0, self.combo_campostatus)
        self.tableWidget.setCellWidget(9,0, self.combo_solver)
        self.tableWidget.setCellWidget(10,0, self.combo_q)

        self.tableWidget.setItem(10 , 0, QtGui.QTableWidgetItem("")) #concentrazione
        self.tableWidget.setItem(11 , 0, QtGui.QTableWidgetItem("0.00005")) #tensore x
        self.tableWidget.setItem(12 , 0, QtGui.QTableWidgetItem("0.00005")) # tensore y
        self.tableWidget.setItem(13 , 0, QtGui.QTableWidgetItem("1")) #ritardo
        self.tableWidget.setItem(14 , 0, QtGui.QTableWidgetItem("0.17")) #porosit??
        self.tableWidget.setItem(15 , 0, QtGui.QTableWidgetItem("")) #error break
        self.tableWidget.setItem(16 , 0, QtGui.QTableWidgetItem("0.01")) #dispersitivit?? trasversla
        self.tableWidget.setItem(17 , 0, QtGui.QTableWidgetItem("0.1")) # dispersivit?? long
        self.tableWidget.setItem(18 , 0, QtGui.QTableWidgetItem("0.0001")) #storavit??
        

       
        self.tableWidget.setCellWidget(19,0, self.combo_phead)
        self.tableWidget.setCellWidget(20,0, self.combo_campophead)

        hbox = QtGui.QHBoxLayout()
        hbox.addStretch(1)
        self.line_output = QtGui.QLineEdit()
        self.line_output.setFixedHeight(25)
        self.saveButton = QtGui.QPushButton("Scegli")
        self.saveButton.setFixedHeight(25)
        hbox.addWidget(self.line_output)
        hbox.addWidget(self.saveButton)
        cellWidget = QtGui.QWidget()
        cellWidget.setLayout(hbox)

        self.spinTime=QtGui.QSpinBox()
        self.spinTime.setRange(1, 9999);
        self.spinTime.setValue(150)

        self.tableWidget.setCellWidget(21,0, self.spinTime)


        self.tableWidget.setCellWidget(22,0, cellWidget)

        self.spinRes=QtGui.QSpinBox()
        self.spinRes.setValue(1)

        self.tableWidget.setCellWidget(23,0, self.spinRes)

        self.tableWidget.setItem(24 , 0, QtGui.QTableWidgetItem("3044")) 


        self.tableWidget.resizeRowsToContents();

        # rowPosition = self.tableWidget.rowCount()
        # self.tableWidget.insertRow(rowPosition)
        # self.tableWidget.setItem(rowPosition , 0, QtGui.QTableWidgetItem("text1"))
        self.tableWidget.setVerticalHeaderLabels((u'Vettoriale sorgente', u'Campo concentrazione sorgente',u'Vettoriale confine*', u'Aquifer top surface map (m)',
                                                  u'Aquifer top surface value (m)',u'Aquifer bottom surface map (m)',u'Aquifer bottom surface value (m)',
                                                  u'Status',u'Campo status',u'Solver',u'Groundwater source and skins (m3/s)',
                                                  u'Tensore di conduttivit?? idraulica X*',
                                                  u'Tensore di conduttivit?? idraulica Y*',u'Fattore di ritardo',
                                                  u'Porosit??*',u'Error break criteria',u'Dispersivit?? trasversale*',
                                                  u'Dispersivit?? longitudinale*',u'Storativit??',u'Mappa altezza piezometrica (m)',
                                                  u'Campo altezza piezometrica',u'Tempo analisi (h)',u'Output file',
                                                  u'Risoluzione',u'Identificativo SRS'))



        self.label_status.setText("In attesa di dati")
        self.label_status.setStyleSheet('color : green; font-weight:bold') 

        self.clear_out_button.clicked.connect(self.reset_output)
        self.save_out_button.clicked.connect(self.esporta_output)

        self.popolacombo()


        self.saveButton.clicked.connect(lambda: self.scegli_file("folder"))
        self.reset_field_button.clicked.connect(self.reset_fields)
        self.buttonBox.accepted.connect(self.run_transport)
        self.actionManuale.activated.connect(self.help)
        self.actionCredits.activated.connect(self.about)
   

        # cbox_source=self.tableWidget.findChild(QComboBox,"combo_source")
        # cbox_conc=self.tableWidget.findChild(QComboBox,"combo_conc")
        # cbox_phead=self.tableWidget.findChild(QComboBox,"combo_phead")
        # cbox_conc=self.tableWidget.findChild(QComboBox,"combo_campophead")

        # pyqtRemoveInputHook()
        # pdb.set_trace()      

        self.combo_source.currentIndexChanged[str].connect(lambda: self.checkfields(self.tableWidget.findChild(QComboBox,"combo_source"),self.tableWidget.findChild(QComboBox,"combo_conc")))
        self.combo_phead.currentIndexChanged[str].connect(lambda: self.checkfields(self.tableWidget.findChild(QComboBox,"combo_phead"),self.tableWidget.findChild(QComboBox,"combo_campophead")))
        self.combo_status.currentIndexChanged[str].connect(lambda: self.checkfields(self.tableWidget.findChild(QComboBox,"combo_status"),self.tableWidget.findChild(QComboBox,"combo_campostatus")))
        #self.button_run.clicked.connect(self.run1)



        self.tabWidget.removeTab(1)

        self.menufile=self.menuFile_6

        self.tesimenu=QAction("Doc. originale",self)    
        self.menufile.addAction(self.tesimenu) 
        self.tesimenu.activated.connect(self.tesi_or) 




    def run1(self):
        # fix_print_with_import
        print("run effettuato")


    def help(self):         
        #self.credits = u"Universit?? della Tuscia\n Viterbo - Italy\nRaffaele Pelorosso, Federica Gobattoni\nDeveloper: Francesco Geri"
        #QMessageBox.about(self.dlg,"Credits", self.credits ) 
        if platform.uname()[0]=="Windows":
            os.system("start "+os.path.dirname(__file__)+"/../tutorial/manuale_envifate_solute.pdf")
        if platform.uname()[0]=="Linux":
            os.system("xdg-open "+os.path.dirname(__file__)+"/../tutorial/manuale_envifate_solute.pdf")
        else:
            os.system("open "+os.path.dirname(__file__)+"/../tutorial/manuale_envifate_solute.pdf")

    def tesi_or(self):         
        #self.credits = u"Universit?? della Tuscia\n Viterbo - Italy\nRaffaele Pelorosso, Federica Gobattoni\nDeveloper: Francesco Geri"
        #QMessageBox.about(self.dlg,"Credits", self.credits ) 
        if platform.uname()[0]=="Windows":
            os.system("start "+os.path.dirname(__file__)+"/../tutorial/gebbert2007_diplom_stroemung_grass_gis.pdf")
        if platform.uname()[0]=="Linux":
            os.system("xdg-open "+os.path.dirname(__file__)+"/../tutorial/gebbert2007_diplom_stroemung_grass_gis.pdf")
        else:
            os.system("open "+os.path.dirname(__file__)+"/../tutorial/gebbert2007_diplom_stroemung_grass_gis.pdf")    

        
    def esporta_output(self):
        try:
            resultmodel=self.console.toPlainText()
            self.f_res_export=open("resultmodel.txt","w")
            self.f_res_export.write(resultmodel.encode('utf8'))
            self.f_res_export.close()
            QMessageBox.information(self,"Info", "File resultmodel.txt exported in the working folder" )            
        except:
            QMessageBox.information(selfs,"Info", "Nothing to export" )  



    def reset_output(self):   
        ret = QMessageBox.warning(self,"Attenzione", "Vuoi davvero eliminare i risultati del modello?",QtGui.QMessageBox.Yes | QtGui.QMessageBox.Default,QtGui.QMessageBox.No,QtGui.QMessageBox.Cancel | QtGui.QMessageBox.Escape)  
        if ret== QtGui.QMessageBox.Yes:
            self.console.clear()
        else:
            return False


    def popolacombo(self):
        self.progressBar.setValue(0)
        self.combo_source.clear()
        self.combo_bound.clear()
        self.combo_top.clear()
        self.combo_phead.clear()
        self.combo_campophead.clear()
        self.combo_bottom.clear()
        self.combo_q.clear()
        self.combo_solver.clear()
        self.combo_status.clear()
        self.combo_campostatus.clear()
        #self.combo_maindirwind.clear()
        self.line_output.clear()

        self.allLayers = self.canvas.layers()
        self.listalayers=dict()
        elementovuoto="No file"


        self.combo_top.addItem(elementovuoto)
        self.combo_bottom.addItem(elementovuoto)
        self.combo_q.addItem(elementovuoto)
        #self.combo_status.addItem(elementovuoto)
        

        for i in self.allLayers:
            self.listalayers[i.name()]=i
            if i.type() == QgsMapLayer.VectorLayer:
                self.combo_source.addItem(str(i.name()))
                self.combo_bound.addItem(str(i.name()))
                self.combo_phead.addItem(str(i.name()))
                self.combo_status.addItem(str(i.name()))
            if i.type()==QgsMapLayer.RasterLayer:
                self.combo_top.addItem(str(i.name()))
                self.combo_bottom.addItem(str(i.name()))
                self.combo_q.addItem(str(i.name()))
                
                


        self.popolafields(self.combo_source,self.combo_conc)
        self.popolafields(self.combo_phead,self.combo_campophead)
        self.popolafields(self.combo_status,self.combo_campostatus)


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


        self.combo_solver.addItem("gauss")
        self.combo_solver.addItem("lu")
        self.combo_solver.addItem("jacobi")
        self.combo_solver.addItem("sor")
        self.combo_solver.addItem("bicgstab")


    

    def reset_fields(self):


        self.console.clear()


        for i in range(11,19):
            self.tableWidget.setItem(i , 0, QtGui.QTableWidgetItem(""))
        self.tableWidget.setItem(4 , 0, QtGui.QTableWidgetItem(""))
        self.tableWidget.setItem(6 , 0, QtGui.QTableWidgetItem(""))
        self.tableWidget.setItem(22 , 0, QtGui.QTableWidgetItem(""))

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
            self.line_output.setText(self.fname)  
        if tipofile=="folder":
            self.folder = QFileDialog.getExistingDirectory(self, "Select Directory")
            self.line_output.setText(self.folder)                 
        # if self.tipofile=="salvacsv":
        #     self.fname = QFileDialog.getSaveFileName(None, 'Save file', '/home','csv files (*.csv);;all files (*.*)')
        #     self.dlg.lineEdit_csv.setText(self.fname)

    def about(self):         
        QMessageBox.about(self, "Credits EnviFate",u"""<p>EnviFate: Open source tool for environmental risk analysis<br />Release 1.0<br />13-1-2017<br />License: GPL v. 3<br /><a href='https://bitbucket.org/fragit/envifate'>Home page plugin</a></p><hr><p>Lavoro svolto nell???ambito del  Progetto  di   ricerca   scientifica  ???Definizione  di   metodi   standard     e  di strumenti applicativi   informatici per   il calcolo degli effetti dei fattori di perturbazione   ai sensi della decisione  2011/484/Ue,  da impiegarsi  nell???ambito  della valutazione di incidenza??? finanziato dalla Regione Veneto. Partner principale ?? il DICAM, Dipartimento di Ingegneria Civile Ambientale e Meccanica dell???Universit?? di Trento (Italia).</p><hr><p>Autori: Francesco Geri, Marco Ciolli</p><p>Universita' di Trento, Trento - Dipartimento di Ingegneria Civile Ambientale e Meccanica (DICAM) <a href="http://www.dicam.unitn.it/">www.dicam.unitn.it/</a></p><hr><p>Consulenti: Paolo Zatelli, Oscar Cainelli</p>""")          


    def checkfields(self,campo_in,campo_out):
        campo_out.clear()
        # self.combo_campophead.clear()
        # self.popolafields(self.combo_source,self.combo_conc)
        self.popolafields(campo_in,campo_out)

    def popolafields(self,combo_in,combo_out):
        vect_source_text=combo_in.currentText()
        if vect_source_text!="":
            #vfields = self.allLayers[mainvect].pendingFields()
            mainvect = self.registry.mapLayersByName( vect_source_text )[0]
            vfields = mainvect.pendingFields()
            # combo_out.addItem("No field")
            for field in vfields:
                combo_out.addItem(field.name())


    def run_transport(self):


        self.text_valuetop=str(self.tableWidget.item(4,0).text())
        self.text_valuebottom=str(self.tableWidget.item(6,0).text())
        self.text_tensorx=str(self.tableWidget.item(11,0).text())
        self.text_tensory=str(self.tableWidget.item(12,0).text())
        self.text_ritardo=str(self.tableWidget.item(13,0).text())
        self.text_porosita=str(self.tableWidget.item(14,0).text())
        self.text_errore=str(self.tableWidget.item(15,0).text())
        self.text_dispt=str(self.tableWidget.item(16,0).text())
        self.text_displ=str(self.tableWidget.item(17,0).text())
        self.text_stor=str(self.tableWidget.item(18,0).text())


        self.text_vector = str(self.combo_source.currentText())
        self.text_area = str(self.combo_bound.currentText())
        self.text_top = str(self.combo_top.currentText())
        self.text_bottom = str(self.combo_bottom.currentText())
        self.text_phead = str(self.combo_phead.currentText())
        self.text_campophead = str(self.combo_campophead.currentText())

        self.text_status = str(self.combo_status.currentText())
        self.text_campostatus = str(self.combo_campostatus.currentText())
        self.text_solver = str(self.combo_solver.currentText())
        self.text_q = str(self.combo_q.currentText())

        self.text_conc = str(self.combo_conc.currentText())

        self.text_srs=str(self.tableWidget.item(24,0).text())

        # try:
        #     self.conc=float(self.text_conc)
        # except Exception, e:
        #     QMessageBox.warning(self,"Warning", u"Errore nella variabile concentrazione dell'inquinante" )
        #     return


        self.res=int(self.spinRes.text())
        self.time=int(self.spinTime.text())*3600


        try:
            self.tensorx=float(self.text_tensorx)
        except Exception as e:
            QMessageBox.warning(self,"Warning", "Entrambi i valori di conduttivit?? idraulica sono obbligatori" )
            return

        try:
            self.tensory=float(self.text_tensory)
        except Exception as e:
            QMessageBox.warning(self,"Warning", "Entrambi i valori di conduttivit?? idraulica sono obbligatori" )
            return


        try:
            self.porosita=float(self.text_porosita)
        except Exception as e:
            QMessageBox.warning(self,"Warning", u"La porosit?? ?? obbligatoria" )
            return

        try:
            self.dispt=float(self.text_dispt)
        except Exception as e:
            QMessageBox.warning(self,"Warning", u"Entrambi i valori di dispersivit?? sono obbligatori" )
            return

        try:
            self.displ=float(self.text_displ)
        except Exception as e:
            QMessageBox.warning(self,"Warning", u"Entrambi i valori di dispersivit?? sono obbligatori" )
            return

        if self.text_stor!='':
            try:
                self.stor=float(self.text_stor)
            except Exception as e:
                QMessageBox.warning(self,"Warning", u"Errore nella variabile storativit??" )
                return
        else:
            self.stor=0.0001

        if self.text_errore!='':
            try:
                self.errore=float(self.text_errore)
            except Exception as e:
                QMessageBox.warning(self,"Warning", "Errore nel break error" )
                return
        else:
            self.errore=0.000001

        if self.text_ritardo!='':
            try:
                self.ritardo=float(self.text_ritardo)
            except Exception as e:
                QMessageBox.warning(self,"Warning", "Errore nel fattore di ritardo" )
                return
        else:
            self.ritardo=1.0

        if self.text_valuetop=='' and self.text_top=='No file':
            QMessageBox.warning(self,"Warning", u"Acquifer top surface ?? obbligatorio" )
            return


        messaggio_top=''
        messaggio_mappatop=''
        if self.text_top!='No file':
            self.top=self.listalayers[self.text_top]
            messaggio_mappatop="Raster acquifer top surface: "+str(self.top)+"\n"
            if not self.top.isValid():
                QMessageBox.warning(self,"Warning", u"Il file raster top acquifer surface non ?? valido" )
                return
        else:
            try:
                self.top=float(self.text_valuetop)
                messaggio_top="Acquifer top surface: "+str(self.top)+"\n"
            except Exception as e:
                QMessageBox.warning(self,"Warning", "Errore nella variabile Acquifer top surface" )
                return



        if self.text_valuebottom=='' and self.text_bottom=='No file':
            QMessageBox.warning(self,"Warning", u"Acquifer bottom surface ?? obbligatorio" )
            return

        messaggio_bottom=''
        messaggio_mappabottom=''

        if self.text_bottom!='No file':
            self.bottom=self.listalayers[self.text_bottom]
            messaggio_mappabottom="Raster acquifer bottom surface: "+str(self.bottom)+"\n"
            if not self.bottom.isValid():
                QMessageBox.warning(self,"Warning", u"Il file raster bottom acquifer surface non ?? valido" )
                return
        else:
            try:
                self.bottom=float(self.text_valuebottom)
                messaggio_bottom="Acquifer bottom surface: "+str(self.bottom)+"\n"
            except Exception as e:
                QMessageBox.warning(self,"Warning", "Errore nella variabile Acquifer bottom surface" )
                return

        messaggio_q=''
        if self.text_q!='No file':
            self.q=self.listalayers[self.text_q]
            messaggio_q="Mappa concentration sources and skins: "+str(self.q)+"\n"
            if not self.q.isValid():
                QMessageBox.warning(self,"Warning", u"Il file raster concentration sources and skins non ?? valido" )
                return
        else:
            self.q=0.0


        self.phead=self.listalayers[self.text_phead]

        if not self.phead.isValid():
            QMessageBox.warning(self,"Warning", "La mappa raster della testa piezometrica non ?? valida" )
            return


        self.status=self.listalayers[self.text_status]

        if not self.status.isValid():
            QMessageBox.warning(self,"Warning", "Il file raster status non ?? valido" )
            return


        # messaggio_status=''
        # if self.text_status!='No file':
        #     self.status=self.listalayers[self.text_status]
        #     messaggio_status="Mappa status: "+str(self.q)+"\n"
        #     if not self.status.isValid():
        #         QMessageBox.warning(self,"Warning", u"Il file raster status non ?? valido" )
        #         return
        # else:
        #     self.status=1

        self.source=self.listalayers[self.text_vector]

        if self.source.wkbType()!=QGis.WKBPoint:
            QMessageBox.warning(self,"Warning", "The source file must have point geometry" )
            return

        self.areastudio=self.listalayers[self.text_area]
     

        if self.areastudio.wkbType()!=QGis.WKBPolygon:
            QMessageBox.warning(self,"Warning", "The boundaries file must have polygon geometry" )
            return


        self.path_output=self.line_output.text()
        if self.path_output=="":
            self.path_output=os.path.dirname(__file__)

        try:
            self.myepsg=float(self.text_srs)
        except Exception as e:
            QMessageBox.warning(self,"Warning", u"Il sistema di riferimento geografico (SRS) ?? obbligatorio" )
            return

        messaggio="Inizio elaborazione trasporto soluto in falda Envifate\n"
        messaggio+="---------------------------\n\n"
        messaggio+="FILE DI INPUT:\n"
        messaggio+="Vettoriale sorgente: "+str(self.text_vector)+"\n"
        messaggio+="Vettoriale confine: "+str(self.text_area)+"\n"
        messaggio+=messaggio_q
        messaggio+=messaggio_mappatop
        messaggio+=messaggio_mappabottom
        # messaggio+=messaggio_status
        messaggio+="VARIABILI:\n"
        messaggio+=u"Concentrazione inquinante: da shapefile"
        messaggio+=u"Solver: "+str(self.combo_solver.currentText())+"\n"
        messaggio+=u"Tensore di conduttivit?? X: "+str(self.text_tensorx)+"\n"
        messaggio+=u"Tensore di conduttivit?? Y: "+str(self.text_tensory)+"\n"
        messaggio+="Fattore di ritardo: "+str(self.text_ritardo)+"\n"
        messaggio+=u"Porosit??: "+str(self.text_porosita)+"\n"
        messaggio+="Error break criteria: "+str(self.text_errore)+"\n"
        messaggio+=u"Dispersivit?? trasversale: "+str(self.text_dispt)+"\n"
        messaggio+=u"Dispersivit?? longitudinale: "+str(self.text_displ)+"\n"
        messaggio+=u"Storavit??: "+str(self.text_stor)+"\n"
        messaggio+=messaggio_top
        messaggio+=messaggio_bottom
        messaggio+="Tempo analisi: "+str(self.spinTime.text())+"\n"
        messaggio+="Risoluzione: "+str(self.res)+"\n\n"
        messaggio+='ALGORITMO UTILIZZATO: modello realizzato sfruttando r.gwflow e r.solute.transport di Grass GIS v.7 (Neteler and Mitasova, 2008). Per i dettagli matematici degli algoritmi vedere https://grass.osgeo.org/gdp/hydrology/gebbert2007_diplom_stroemung_grass_gis.pdf\n\n'
        messaggio+="---------------------------\n\n"

        self.console.appendPlainText(messaggio)


        self.label_status.setText("Preparazione dati")
        self.label_status.setStyleSheet('color : #e8b445;font-weight:bold')    
        #rasterizzazione layer boundaries
        #raster_fn = 'test.tif'

        # pyqtRemoveInputHook()
        # pdb.set_trace() 
        self.progressBar.setValue(50)


        start_time = time.time()  

        self.label_status.setText("Processing data")
        self.label_status.setStyleSheet('color : #e8b445;font-weight:bold')



        # say hello


        self.areastudio_path=self.areastudio.dataProvider().dataSourceUri().split('|')
        self.phead_path=self.phead.dataProvider().dataSourceUri().split('|')
        self.source_path=self.source.dataProvider().dataSourceUri().split('|')
        self.status_path=self.status.dataProvider().dataSourceUri().split('|')

        #g.proj -c epsg=3003
        run_command("g.proj", epsg=self.myepsg, flags = 'c')

        run_command("v.in.ogr", input=self.areastudio_path[0], output="areastudio",overwrite=True, flags = 'o')
        run_command("g.region", res=self.res, vector="areastudio")
        run_command("v.to.rast", input="areastudio", output="areastudio",overwrite=True, use = 'val')


        run_command("v.in.ogr", input=self.phead_path[0], output="phead",overwrite=True, flags = 'o')
        run_command("v.to.rast", input="phead", output="phead",overwrite=True, use = 'attr', attribute_column=self.text_campophead)

        run_command("v.in.ogr", input=self.status_path[0], output="status",overwrite=True, flags = 'o')
        run_command("v.to.rast", input="status", output="status",overwrite=True, use = 'attr', attribute_column=self.text_campostatus)

        # pyqtRemoveInputHook()
        # pdb.set_trace() 


        run_command("v.in.ogr", input=self.source_path[0], output="source",overwrite=True, flags = 'o')
        run_command("v.to.rast", input="source", output="source",overwrite=True, use = 'attr', attribute_column=self.text_conc)
        run_command("r.null", map="source", null=0)




        if self.q==0.0:
            run_command("r.mapcalc", expression="well = 0",overwrite=True)
        else:
            q_path=self.q.dataProvider().dataSourceUri().split('|')
            run_command("r.in.gdal", input=q_path[0], output="well",overwrite=True,flags = 'o')            


        if self.status==1:
            run_command("r.mapcalc", expression="status = 1",overwrite=True)

        #run_command("r.mapcalc", expression="status = if(col() == 1 || col() == 200 , 2, 1)")
        #run_command("r.mapcalc", expression="well = 0")
        ex_hydx="hydcondx = "+'{0:.10f}'.format(self.tensorx)
        ex_hydy="hydcondy = "+'{0:.10f}'.format(self.tensory)
        run_command("r.mapcalc", expression=ex_hydx,overwrite=True)
        run_command("r.mapcalc", expression=ex_hydy,overwrite=True)
        run_command("r.mapcalc", expression="recharge = 0",overwrite=True)

        run_command("r.mapcalc", expression="top_conf = "+str(self.top),overwrite=True)
        run_command("r.mapcalc", expression="bottom = "+str(self.bottom),overwrite=True)

        run_command("r.mapcalc", expression="poros = "+str(self.porosita),overwrite=True)
        run_command("r.mapcalc", expression="syield = "+str(self.stor),overwrite=True)
        run_command("r.mapcalc", expression="null = 0.0",overwrite=True)
        run_command("r.mapcalc", expression="cs = 0.0",overwrite=True)
        run_command("r.mapcalc", expression="diff = 0.0000001",overwrite=True)
        run_command("r.mapcalc", expression="R = "+str(self.ritardo),overwrite=True)

        run_command("r.gwflow", solver="cg", top="top_conf", bottom="bottom", phead="phead",\
          status="status", hc_x="hydcondx", hc_y="hydcondy", q="well", s="syield",\
          recharge="recharge", output="gwresult_conf", dt=self.time, type="confined",overwrite=True)  

        run_command("r.solute.transport", solver=self.text_solver, top="top_conf",\
          bottom="bottom", phead="gwresult_conf", status="status", hc_x="hydcondx", hc_y="hydcondy",\
          rd="R", cs="cs", q="well", nf="poros", output="stresult_conf", dt=self.time, diff_x="diff",\
          diff_y="diff", c="source", al=self.displ, at=self.dispt,overwrite=True)     



        run_command("r.out.gdal", input="gwresult_conf", output=self.path_output+"/gwflow.tif",overwrite=True,format = 'GTiff')
        run_command("r.out.gdal", input="stresult_conf", output=self.path_output+"/solute.tif",overwrite=True,format = 'GTiff')


        base_gwflow_name=os.path.basename(self.path_output+"/gwflow.tif")
        gwflow_name=os.path.splitext(base_gwflow_name)[0]
        gwflow_file=self.iface.addRasterLayer(self.path_output+"/gwflow.tif", gwflow_name)   

        base_solute_name=os.path.basename(self.path_output+"/solute.tif")
        solute_name=os.path.splitext(base_solute_name)[0]
        solute_file=self.iface.addRasterLayer(self.path_output+"/solute.tif", solute_name)              
        # max_progress=rows*cols
        # self.progressBar.setMaximum(max_progress)  



        tempoanalisi=time.time() - start_time
        tempostimato=time.strftime("%H:%M:%S", time.gmtime(tempoanalisi))
        messaggio="---------------------------------\n"
        messaggio+="Fine modellazione\n"
        messaggio+="\nTempo di analisi: "+tempostimato+"\n"
        messaggio+="---------------------------------\n\n"
        self.progressBar.setValue(100)
        self.console.appendPlainText(messaggio)       

        self.label_status.setText("In attesa di dati")
        self.label_status.setStyleSheet('color : green; font-weight:bold')  

        #shutil.rmtree(location_path)        