module Transport

# -*- coding: utf-8 -*-
#=
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
=#
#=
    def help(self):         
        #self.credits = u"Università della Tuscia\n Viterbo - Italy\nRaffaele Pelorosso, Federica Gobattoni\nDeveloper: Francesco Geri"
        #QMessageBox.about(self.dlg,"Credits", self.credits ) 
        if platform.uname()[0]=="Windows":
            os.system("start "+os.path.dirname(__file__)+"/../tutorial/manuale_envifate_solute.pdf")
        if platform.uname()[0]=="Linux":
            os.system("xdg-open "+os.path.dirname(__file__)+"/../tutorial/manuale_envifate_solute.pdf")
        else:
            os.system("open "+os.path.dirname(__file__)+"/../tutorial/manuale_envifate_solute.pdf")

    def tesi_or(self):         
        #self.credits = u"Università della Tuscia\n Viterbo - Italy\nRaffaele Pelorosso, Federica Gobattoni\nDeveloper: Francesco Geri"
        #QMessageBox.about(self.dlg,"Credits", self.credits ) 
        if platform.uname()[0]=="Windows":
            os.system("start "+os.path.dirname(__file__)+"/../tutorial/gebbert2007_diplom_stroemung_grass_gis.pdf")
        if platform.uname()[0]=="Linux":
            os.system("xdg-open "+os.path.dirname(__file__)+"/../tutorial/gebbert2007_diplom_stroemung_grass_gis.pdf")
        else:
            os.system("open "+os.path.dirname(__file__)+"/../tutorial/gebbert2007_diplom_stroemung_grass_gis.pdf")    


        self.popolacombo()
=#


function run_transport( top, tensorx::Real, tensory::Real, porosity::Real, dispersiveness_t::Real, dispersiveness_l::Real, storativity::Real=0.0001,
                        break_error::Real=0.000001, delay::Real=1.0 )




    if self.text_valuetop=='' and self.text_top=='No file':
        QMessageBox.warning(self,"Warning", u"Acquifer top surface è obbligatorio" )
        return


    messaggio_top=''
    messaggio_mappatop=''
    if self.text_top!='No file':
        self.top=self.listalayers[self.text_top]
        messaggio_mappatop="Raster acquifer top surface: "+str(self.top)+"\n"
        if not self.top.isValid():
            QMessageBox.warning(self,"Warning", u"Il file raster top acquifer surface non è valido" )
            return
    else:
        try:
            self.top=float(self.text_valuetop)
            messaggio_top="Acquifer top surface: "+str(self.top)+"\n"
        except Exception as e:
            QMessageBox.warning(self,"Warning", "Errore nella variabile Acquifer top surface" )
            return



    if self.text_valuebottom=='' and self.text_bottom=='No file':
        QMessageBox.warning(self,"Warning", u"Acquifer bottom surface è obbligatorio" )
        return

    messaggio_bottom=''
    messaggio_mappabottom=''

    if self.text_bottom!='No file':
        self.bottom=self.listalayers[self.text_bottom]
        messaggio_mappabottom="Raster acquifer bottom surface: "+str(self.bottom)+"\n"
        if not self.bottom.isValid():
            QMessageBox.warning(self,"Warning", u"Il file raster bottom acquifer surface non è valido" )
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
            QMessageBox.warning(self,"Warning", u"Il file raster concentration sources and skins non è valido" )
            return
    else:
        self.q=0.0


    self.phead=self.listalayers[self.text_phead]

    if not self.phead.isValid():
        QMessageBox.warning(self,"Warning", "La mappa raster della testa piezometrica non è valida" )
        return


    self.status=self.listalayers[self.text_status]

    if not self.status.isValid():
        QMessageBox.warning(self,"Warning", "Il file raster status non è valido" )
        return



 

    if self.source.wkbType()!=QGis.WKBPoint:
        QMessageBox.warning(self,"Warning", "The source file must have point geometry" )
        return

    

    if self.areastudio.wkbType()!=QGis.WKBPolygon:
        QMessageBox.warning(self,"Warning", "The boundaries file must have polygon geometry" )
        return


    self.path_output=self.line_output.text()
    if self.path_output=="":
        self.path_output=os.path.dirname(__file__)

    try:
        self.myepsg=float(self.text_srs)
    except Exception as e:
        QMessageBox.warning(self,"Warning", u"Il sistema di riferimento geografico (SRS) è obbligatorio" )
        return










    self.text_valuetop=str(self.tableWidget.item(4,0).text())
    self.text_valuebottom=str(self.tableWidget.item(6,0).text())

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

    self.res=int(self.spinRes.text())
    self.time=int(self.spinTime.text())*3600




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
    messaggio+=u"Tensore di conduttività X: "+str(self.text_tensorx)+"\n"
    messaggio+=u"Tensore di conduttività Y: "+str(self.text_tensory)+"\n"
    messaggio+="Fattore di delay: "+str(self.text_delay)+"\n"
    messaggio+=u"Porosità: "+str(self.text_porosity)+"\n"
    messaggio+="Error break criteria: "+str(self.text_break_error)+"\n"
    messaggio+=u"Dispersività trasversale: "+str(self.text_dispersiveness_t)+"\n"
    messaggio+=u"Dispersività longitudinale: "+str(self.text_dispersiveness_l)+"\n"
    messaggio+=u"Storavità: "+str(self.text_storativity)+"\n"
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

    run_command("r.mapcalc", expression="poros = "+str(self.porosity),overwrite=True)
    run_command("r.mapcalc", expression="syield = "+str(self.storativity),overwrite=True)
    run_command("r.mapcalc", expression="null = 0.0",overwrite=True)
    run_command("r.mapcalc", expression="cs = 0.0",overwrite=True)
    run_command("r.mapcalc", expression="diff = 0.0000001",overwrite=True)
    run_command("r.mapcalc", expression="R = "+str(self.delay),overwrite=True)

    run_command("r.gwflow", solver="cg", top="top_conf", bottom="bottom", phead="phead",\
      status="status", hc_x="hydcondx", hc_y="hydcondy", q="well", s="syield",\
      recharge="recharge", output="gwresult_conf", dt=self.time, type="confined",overwrite=True)  

    run_command("r.solute.transport", solver=self.text_solver, top="top_conf",\
      bottom="bottom", phead="gwresult_conf", status="status", hc_x="hydcondx", hc_y="hydcondy",\
      rd="R", cs="cs", q="well", nf="poros", output="stresult_conf", dt=self.time, diff_x="diff",\
      diff_y="diff", c="source", al=self.dispersiveness_l, at=self.dispersiveness_t,overwrite=True)     



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