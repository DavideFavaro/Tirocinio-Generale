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
from qgis.PyQt import QtCore, QtGui

from PyQt5.QtCore import QSettings, QTranslator, QCoreApplication, Qt, QObject, pyqtSignal, pyqtRemoveInputHook, QVariant
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

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar

import matplotlib.pyplot as plt

import random
import math
import time
import platform

import pdb

from envifate_dialog import EnviDialog 

from configuration_dialog import ConfigurationDialog

sys.path.append( os.path.dirname(__file__)+"/../library" )

import functions, river, do_setting

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

        self.label_title.setText("Dispersione fluviale")
        self.tableWidget.setRowCount(15)
        self.tableWidget.horizontalHeader().setStretchLastSection(True)
        #self.tableWidget.horizontalHeaderItem(0).setText("newHeader")
        self.combovector1 = QComboBox()
        self.combosource1 = QComboBox()
        self.comboslope1 = QComboBox()
        self.combodem1 = QComboBox()
        self.tableWidget.setCellWidget(0,0, self.combovector1)
        self.tableWidget.setCellWidget(1,0, self.combosource1)
        self.tableWidget.setCellWidget(2,0, self.comboslope1)
        self.tableWidget.setCellWidget(3,0, self.combodem1)
        
        self.tableWidget.setItem(4 , 0, QTableWidgetItem(""))
        self.tableWidget.setItem(5 , 0, QTableWidgetItem(""))
        self.tableWidget.setItem(6 , 0, QTableWidgetItem(""))
        self.tableWidget.setItem(7 , 0, QTableWidgetItem(""))
        self.tableWidget.setItem(8 , 0, QTableWidgetItem(""))
        self.tableWidget.setItem(9 , 0, QTableWidgetItem(""))

        self.spintime1=QSpinBox()
        self.spintime2=QSpinBox()
        self.spintime_step=QSpinBox()
        self.spinres1=QSpinBox()

        self.tableWidget.setCellWidget(10,0, self.spintime1)
        self.tableWidget.setCellWidget(11,0, self.spintime2)
        self.tableWidget.setCellWidget(12,0, self.spintime_step)
        self.tableWidget.setCellWidget(13,0, self.spinres1)

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

        # self.tableWidget.setCellWidget(14,0, cellWidget)
        self.tableWidget.setItem(14 , 0, QTableWidgetItem(""))


        self.spinres1.setValue(10)
        self.spintime1.setValue(5)
        self.spintime2.setValue(10)
        self.spintime_step.setValue(1)

        # rowPosition = self.tableWidget.rowCount()
        # self.tableWidget.insertRow(rowPosition)
        # self.tableWidget.setItem(rowPosition , 0, QtGui.QTableWidgetItem("text1"))
        self.tableWidget.setVerticalHeaderLabels((u'Vettoriale fiume*', u'Vettoriale sorgente*',u'Mappa pendenza',
                                                  u'DTM*',u'Massa inquinante* (g)',u'Contorno bagnato (m)',u'Raggio idraulico medio (m)',
                                                  u'Coff. di diffusione di Fick',u'Coeff. di scabrezza',u'Coeff. di decadimento',
                                                  u'Tempo start (min)',u'Tempo end (min)',u'Tempo intervallo (min)',u'Risoluzione (m)',u'Output filename'))
        

        self.label_status.setText("In attesa di dati")
        self.label_status.setStyleSheet('color : green; font-weight:bold') 

        self.figure = plt.figure()        
        self.canvas_mat = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas_mat, self)
        self.layout_mat.addWidget(self.toolbar)
        self.layout_mat_2.addWidget(self.canvas_mat) 

        self.clear_out_button.clicked.connect(self.reset_output)
        self.save_out_button.clicked.connect(self.esporta_output)

        self.popolacombo()

        #self.saveButton.clicked.connect(lambda: self.scegli_file("salvavector"))
        self.reset_field_button.clicked.connect(self.reset_fields)
        self.button_reset_draw.clicked.connect(self.reset_draw)
        self.buttonBox.accepted.connect(self.run_river)
        self.buttonBox_2.accepted.connect(self.run_temp_river)

        # self.web = QWebView()
        # self.web.load(QUrl("https://grass.osgeo.org/grass70/manuals/addons/r.green.biomassfor.theoretical.html"))
        # self.web_layout.addWidget(self.web)
        self.actionManuale.triggered.connect(self.help)
        self.actionCredits.triggered.connect(self.about)
        self.actionSetting.triggered.connect(self.configuration)

    def about(self):         
        QMessageBox.about(self, "Credits EnviFate",u"""<p>EnviFate: Open source tool for environmental risk analysis<br />Release 1.0<br />13-1-2017<br />License: GPL v. 3<br /><a href='https://bitbucket.org/fragit/envifate'>Home page plugin</a></p><hr><p>Lavoro svolto nell’ambito del  Progetto  di   ricerca   scientifica  “Definizione  di   metodi   standard     e  di strumenti applicativi   informatici per   il calcolo degli effetti dei fattori di perturbazione   ai sensi della decisione  2011/484/Ue,  da impiegarsi  nell’ambito  della valutazione di incidenza” finanziato dalla Regione Veneto. Partner principale è il DICAM, Dipartimento di Ingegneria Civile Ambientale e Meccanica dell’Università di Trento (Italia).</p><hr><p>Autori: Francesco Geri, Marco Ciolli</p><p>Universita' di Trento, Trento - Dipartimento di Ingegneria Civile Ambientale e Meccanica (DICAM) <a href="http://www.dicam.unitn.it/">www.dicam.unitn.it/</a></p><hr><p>Consulenti: Paolo Zatelli, Oscar Cainelli</p>""")          


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


    def configuration(self):
        d = do_setting.Dialog(self.iface)
        d.show()
        d.exec_()

    def help(self):         
        #self.credits = u"Università della Tuscia\n Viterbo - Italy\nRaffaele Pelorosso, Federica Gobattoni\nDeveloper: Francesco Geri"
        #QMessageBox.about(self.dlg,"Credits", self.credits ) 
        if platform.uname()[0]=="Windows":
            os.system("start "+os.path.dirname(__file__)+"/../tutorial/manuale_envifate_dispersione_fluviale.pdf")
        if platform.uname()[0]=="Linux":
            os.system("xdg-open "+os.path.dirname(__file__)+"/../tutorial/manuale_envifate_dispersione_fluviale.pdf")
        else:
            os.system("open "+os.path.dirname(__file__)+"/../tutorial/manuale_envifate_dispersione_fluviale.pdf")


    def run1(self):
        # fix_print_with_import
        print("run effettuato")
        
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
        self.progressBar_2.setValue(0)
        # self.spinstart.setValue(0)
        # self.spinend.setValue(0)
        # self.spindistance.setValue(0)

        self.combovector1.clear()
        self.combosource1.clear()
        self.combodem1.clear()
        self.comboslope1.clear()

        self.allLayers = self.canvas.layers()
        self.listalayers=dict()
        #elementovuoto="No required"
        for i in self.allLayers:
            self.listalayers[i.name()]=i
            if i.type() == QgsMapLayer.VectorLayer:
                self.combovector1.addItem(str(i.name()))
                self.combosource1.addItem(str(i.name()))
            if i.type()==QgsMapLayer.RasterLayer:
                self.comboslope1.addItem(str(i.name()))
                self.combodem1.addItem(str(i.name()))




    def reset_fields(self):

        self.tableWidget.setItem(4 , 0, QTableWidgetItem(""))
        self.tableWidget.setItem(5 , 0, QTableWidgetItem(""))
        self.tableWidget.setItem(6 , 0, QTableWidgetItem(""))
        self.tableWidget.setItem(7 , 0, QTableWidgetItem(""))
        self.tableWidget.setItem(8 , 0, QTableWidgetItem(""))
        self.tableWidget.setItem(9 , 0, QTableWidgetItem(""))
        self.tableWidget.setItem(14 , 0, QTableWidgetItem(""))
        self.tableWidget.setItem(15 , 0, QTableWidgetItem(""))


        self.popolacombo()

        # pyqtRemoveInputHook()
        # pdb.set_trace() 
        act_layers = [layer for layer in QgsProject.instance().mapLayers().values()]
        for act_layer in act_layers:
            if str(act_layer.name())=="Concentration":
                self.tab_2.setEnabled(True)



    def reset_draw(self):
        pass


    def scegli_file(self,tipofile):
        if tipofile=="salvaraster":
            self.fname = QFileDialog.getSaveFileName(None, 'Save file', '/home','all files (*.*)')
            self.line_output.setText(self.fname)
        if tipofile=="salvavector":
            self.fname = QFileDialog.getSaveFileName(None, 'Save file', '/home','shp files (*.shp)')
            self.line_output.setText(self.fname[0]+".shp")      
        if tipofile=="salvacsv":
            self.fname = QFileDialog.getSaveFileName(None, 'Save file', '/home','csv files (*.csv);;all files (*.*)')
            self.dlg.lineEdit_csv.setText(self.fname)

    def get_variabili(self):
        # inizio recupero e controllo dati da form

        # pyqtRemoveInputHook()
        # pdb.set_trace() 

        self.text_line_conc=str(self.tableWidget.item(4,0).text())
        self.text_line_sez=str(self.tableWidget.item(5,0).text())
        self.text_line_radius=str(self.tableWidget.item(6,0).text())
        self.text_line_flickian=str(self.tableWidget.item(7,0).text())
        self.text_line_manning=str(self.tableWidget.item(8,0).text())
        self.text_line_k=str(self.tableWidget.item(9,0).text())

        self.outputname=str(self.tableWidget.item(14,0).text())

        #self.srid=str(self.tableWidget.item(15,0).text())
        self.text_river= str(self.combovector1.currentText())
        self.text_source = str(self.combosource1.currentText())
        self.text_slope = str(self.comboslope1.currentText())
        self.text_dem = str(self.combodem1.currentText())

        self.res=int(self.spinres1.text())
        self.min=int(self.spintime1.text())
        self.min_end=int(self.spintime2.text())
        self.min_int=int(self.spintime_step.text())

        self.sec=self.min*60
        self.sec_end=self.min_end*60
        self.sec_int=self.min_int*60



        if self.text_line_flickian!='':
            try:
                self.fickian_x=float(self.text_line_flickian)
            except Exception as e:
                QMessageBox.warning(self,"Warning", u"Errore nel coefficiente di trasporto Flickian" )
                return False
        else:
            self.fickian_x=0.05

        if self.text_line_sez!='':
            try:
                self.w=float(self.text_line_sez)
            except Exception as e:
                QMessageBox.warning(self,"Warning", u"Errore nel sezione idraulica" )
                return False
        else:
            self.w=1


        if self.text_line_k!='':
            try:
                self.k=float(self.text_line_k)
            except Exception as e:
                QMessageBox.warning(self,"Warning", u"Errore nel coefficiente di decadimento" )
                return False
        else:
            self.k=0

        if self.text_line_manning!='':
            try:
                self.manning=float(self.text_line_manning)
            except Exception as e:
                QMessageBox.warning(self,"Warning", u"Errore nel coefficiente di Manning" )
                return False
        else:
            self.manning=0.050



        try:
            self.radius=float(self.text_line_radius)
            #self.lf=float(self.text_line_acquifer_depth)
            #self.dz=float(self.text_line_sourcethick)
        except Exception as e:
            QMessageBox.warning(self,"Warning", u"Il raggio idraulico medio è obbligatorio" )
            return False

        try:
            self.C=float(self.text_line_conc)
            #self.lf=float(self.text_line_acquifer_depth)
            #self.dz=float(self.text_line_sourcethick)
        except Exception as e:
            QMessageBox.warning(self,"Warning", u"Concentration is mandatory" )
            return False


        self.river=self.listalayers[self.text_river]
        # pyqtRemoveInputHook()
        # pdb.set_trace() 
        if self.river.wkbType()!=5:
            QMessageBox.warning(self,"Warning", u"The river file must have line geometry" )
            return False


        self.source=self.listalayers[self.text_source]

        if self.source.wkbType()!=1:
            QMessageBox.warning(self,"Warning", u"The source file must have point geometry" )
            return False

        self.slope=self.listalayers[self.text_slope]

        if not self.slope.isValid():
            QMessageBox.warning(self,"Warning", u"La mappa delle pendenza è mancante o non valida" )
            return False


        self.dem=self.listalayers[self.text_dem]

        if not self.slope.isValid():
            QMessageBox.warning(self,"Warning", u"Il DEM è mancante o non valido" )
            return False


        if self.min>=self.min_end:
            QMessageBox.warning(self,"Warning", u"Il tempo end non può essere minore o uguale al tempo start" )
            return False             



        if self.river.crs().authid()!=self.source.crs().authid() or self.slope.crs().authid()!=self.source.crs().authid() or self.dem.crs().authid()!=self.source.crs().authid():
            QMessageBox.warning(self,"Warning", "Errore: i sistemi di riferimento non sono uniformi. Impossibile continuare con l'analisi." )
            return

        self.refsys=self.source.crs().authid().split(':')[1]

    def euclidean_distance(self,x1,y1,x2,y2):
        return math.sqrt((x2-x1)**2 + (y2-y1)**2)

    def run_river(self):

        checkvariabili=self.get_variabili()
        if checkvariabili==False:
            return

        self.label_status.setText("Preparazione dati")
        self.label_status.setStyleSheet('color : #e8b445;font-weight:bold')

        messaggio="Inizio elaborazione dispersione fiumi Envifate\n"
        messaggio+="---------------------------\n\n"
        messaggio+="FILE DI INPUT:\n"
        messaggio+="Vettoriale sorgente: "+str(self.text_source)+"\n"
        messaggio+="DTM: "+str(self.text_dem)+"\n"
        messaggio+="Mappa pendenza: "+str(self.text_slope)+"\n\n"

        messaggio+="VARIABILI:\n"
        messaggio+="Sezione bagnata: "+str(self.text_line_sez)+"\n"
        messaggio+="Raggio idraulico medio: "+str(self.text_line_radius)+"\n"
        messaggio+="Coefficiente Fickian: "+str(self.fickian_x)+"\n"
        messaggio+="Coefficiente decadimento lamda: "+str(self.k)+"\n"
        messaggio+="Coefficiente di scabrezza: "+str(self.manning)+"\n"
        messaggio+="Massa inquinante: "+str(self.C)+"\n"
        messaggio+="Tempo iniziale dell'analisi: "+str(self.min)+"\n"
        messaggio+="Tempo finale dell'iniezione: "+str(self.min_end)+"\n"
        messaggio+="Intervallo temporale: "+str(self.min_int)+"\n"
        messaggio+="Risoluzione: "+str(self.res)+"\n\n"
        messaggio+='ALGORITMO UTILIZZATO: Fickian Mixing Process (Hemond, Harold F., and Elizabeth J. Fechner. Chemical fate and transport in the environment. Elsevier, 2014.)\n\n'
        messaggio+="---------------------------\n\n"
        self.console.appendPlainText(messaggio)


        #calcolo ciclo intervallo temporale di analisi
        cicli=(int(self.min_end)-int(self.min))/int(self.min_int)
        cicli+=1

        feature = next(self.river.getFeatures())
        geomfeature = feature.geometry()
        features = self.river.getFeatures()

        s_feature = next(self.source.getFeatures())
        s_geom =  QgsPoint(s_feature.geometry().asPoint())
        x_source=s_feature.geometry().asPoint()[0]
        y_source=s_feature.geometry().asPoint()[1]   


        firstpoint=geomfeature.interpolate(0)

        xfirst=firstpoint.asPoint()[0]
        yfirst=firstpoint.asPoint()[1]
 

        demfirstpoint=self.dem.dataProvider().identify(QgsPointXY(xfirst, yfirst),QgsRaster.IdentifyFormatValue).results()[1]
        demsource=self.dem.dataProvider().identify(QgsPointXY(x_source, y_source),QgsRaster.IdentifyFormatValue).results()[1]




        if demfirstpoint>=demsource:
            trend=1
            old_x=x_source
            old_y=y_source
        else:
            trend=0
            old_x=xfirst
            old_y=yfirst           

     
        #pyqtRemoveInputHook()

        
        for f in features:
            geom = f.geometry()

            length = geom.length()
            #currentdistance = self.res
            currentdistance=1
            avanzamento=1

            realdistance=0

            feats = [] 
            featlines=[] 

            self.list_result=[]
            # self.list_dist=[]

            max_progress=length/currentdistance
            self.progressBar.setMaximum(max_progress)  
            start_time = time.time()  

            index_progress=0 
            count_index=0
            # pyqtRemoveInputHook()
            # pdb.set_trace()

            if self.outputname=="":
                self.outputname="concentrazione"


            vl = QgsVectorLayer("Point?crs=EPSG:"+self.refsys,self.outputname, "memory")
            #vl.setCrs(QgsCoordinateReferenceSystem(3003))
            # vl.setCrs(crspoint)           
            pr = vl.dataProvider()  
            prfield=pr.addAttributes( [ QgsField("distance", QVariant.Int) ] )            
            prfield2=pr.addAttributes( [ QgsField("vmedia", QVariant.Double) ] )


            list_vmedia=[]

            vline = QgsVectorLayer("LineString?crs=EPSG:"+self.refsys, self.outputname, "memory")
            # crsline = vline.crs()
            # crsline.createFromId(self.srid)
            # vline.setCrs(crsline)
            prline = vline.dataProvider()
            # vline.startEditing()
            prlfield=prline.addAttributes( [ QgsField("distance", QVariant.Int) ] )            
            prlfield2=prline.addAttributes( [ QgsField("vmedia", QVariant.Double) ] )

            sec_cicli=self.sec
            checknumcampo=0
            while sec_cicli<=self.sec_end:  
                checknumcampo+=1
                nomecampo="conc"+str(int(sec_cicli)/60)              
                prfield1=pr.addAttributes( [ QgsField(nomecampo, QVariant.Double) ] )
                prlfield1=prline.addAttributes( [ QgsField(nomecampo, QVariant.Double) ] )
                sec_cicli=sec_cicli+self.sec_int


            

            vl.updateFields()
            vline.updateFields()
            fetline = QgsFeature()
            fetline.setGeometry( QgsGeometry.fromPolyline( [s_geom] ))
            prline.addFeatures( [ fetline ] )
            geomline = fetline.geometry()

            # max_value=0
            # max_value_dist=0

            if trend==1:
                controllo=0
            else:
                controllo=1

                
            self.label_status.setText("Processing data")
            self.label_status.setStyleSheet('color : #e8b445;font-weight:bold')


            while currentdistance < length: 



                index_progress+=1
                self.progressBar.setValue(index_progress)
                
                point = geom.interpolate(currentdistance)


                x=point.asPoint()[0]
                y=point.asPoint()[1]




                if trend==1:
                    if self.euclidean_distance(x_source,y_source,x,y)<=1:
                        controllo=1
                        # avanzamento=self.res

                if trend==0:
                    if self.euclidean_distance(x_source,y_source,x,y)<=1:
                        controllo=0


                if controllo==1:


                    if avanzamento==self.res:

                        realdistance=realdistance+self.res
                        count_index+=1

                        z=self.slope.dataProvider().identify(QgsPointXY(x, y),QgsRaster.IdentifyFormatValue)
                        zresult=z.results()[1]

                        v_inst=(math.pow(self.radius,0.666667)*math.sqrt(zresult/100))*(self.manning)
                        
                        list_vmedia.append(v_inst)

                        self.vmedia=sum(list_vmedia)/count_index

                        sec=self.sec

                        fet = QgsFeature()
                        fet.initAttributes(2+cicli)


                        fetline = QgsFeature()
                        fetline.initAttributes(2+cicli)

                        fet.setAttribute(0,realdistance)
                        fet.setAttribute(1,self.vmedia)

                        fetline.setAttribute(0,realdistance)
                        fetline.setAttribute(1,self.vmedia) 

                        ciclo=1

                        while sec<=self.sec_end:                            

                            element=river.river(self.C,sec,realdistance,self.fickian_x,self.vmedia,self.w,self.k)                    
                            Cfinal=element.calc_concentration()
                            if ciclo==1:
                                self.list_result.append(Cfinal)
                            fet.setAttribute(1+ciclo,Cfinal)
                            fetline.setAttribute(1+ciclo,Cfinal)
                            ciclo=ciclo+1
                            sec=sec+self.sec_int

                            # pdb.set_trace()

                            # if Cfinal>max_value:
                            #     max_value=Cfinal
                            #     max_value_dist=realdistance

                            
                            # self.list_dist.append(realdistance)

                        
                        vl.updateFeature(fet)
                        fet.setGeometry(point)

          
                        fetline.setGeometry( QgsGeometry.fromPolyline( [QgsPoint(old_x,old_y),QgsPoint(x,y)] ))
                        
                        vline.updateFeature(fetline)

                        feats.append(fet)
                        featlines.append(fetline)

                        old_x=x
                        old_y=y

                        avanzamento=0
                    avanzamento=avanzamento+1




                currentdistance = currentdistance + 1
                
            self.label_status.setText("Preparazione output")
            self.label_status.setStyleSheet('color : #e8b445;font-weight:bold')

            pr.addFeatures(feats)
            prline.addFeatures(featlines)
            #vl.updateExtents()
            vl.updateFields()
            vline.updateFields()

            QgsProject.instance().addMapLayer(vl)
            QgsProject.instance().addMapLayer(vline)

            tempoanalisi=time.time() - start_time
            tempostimato=time.strftime("%H:%M:%S", time.gmtime(tempoanalisi))
            messaggio="---------------------------------\n"
            messaggio+="Fine modellazione\n"
            messaggio+="\nTempo di analisi: "+tempostimato+"\n"
            messaggio+="---------------------------------\n\n"
            self.console.appendPlainText(messaggio)             

            self.tab_2.setEnabled(True)

            self.run_temp_river()

        self.label_status.setText("In attesa di dati")
        self.label_status.setStyleSheet('color : green; font-weight:bold') 


    def run_temp_river(self):
        


        ax1f1 = self.figure.add_subplot(111)
        ax1f1.set_xlabel('distanza (cella * risoluzione)')
        ax1f1.set_ylabel('concentrazione')
        plt.tight_layout()

        ax1f1.plot(self.list_result)


        self.canvas_mat.draw() 