module Thermic

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
            os.system("start "+os.path.dirname(__file__)+"/../tutorial/manuale_envifate_ruscellamento.pdf")
        if platform.uname()[0]=="Linux":
            os.system("xdg-open "+os.path.dirname(__file__)+"/../tutorial/manuale_envifate_ruscellamento.pdf")
        else:
            os.system("open "+os.path.dirname(__file__)+"/../tutorial/manuale_envifate_ruscellamento.pdf")
=#
#=
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
=#



import ArchGDAL as agd

function run_thermic( dem, source, river,  output_path::AbstractString=".\\")
    self.text_vector = str(self.combo_source.currentText())
    self.text_river = str(self.combo_river.currentText())
    self.text_dem = str(self.combo_dem.currentText())
    self.text_source_q = str(self.combo_source_q.currentText())
    self.text_source_t = str(self.combo_source_t.currentText())
    self.text_river_q = str(self.combo_river_q.currentText())
    self.text_river_t = str(self.combo_river_t.currentText())



    if agd.geomdim(source) != 0
        throw(DomainError(source, "`source` must be a point"))
    end

    if agd.geomdim(river) != 1
        throw(DomainError(source, "`river` must be a line"))
    end

    layers = agd.getlayers(river)
    refsys = agd.getspatialref(source)

    if agd.importWKT(agd.getproj(dtm)) !=  refsys ||  agd.getspatialref(layers) != refsys
        throw(DomainError("The reference systems are not uniform. Aborting analysis." ))
    end



    #   start_time = time.time()

 """ NON SO DA DOVE PRENDA I DATI """
    buffer_sorgente = processing.run("native:buffer", {"INPUT": self.source, "DISTANCE": 5, "OUTPUT":"memory:"})
    sorgente = buffer_sorgente["OUTPUT"]
 """                              """

    for f in agd.getfeature(sorgente)
        for a in agd.getfeature(river)
            geom_a = agd.getgeom(a)
            geom_f = agd.getgeom(f)
            if agd.intersects(geom_a, geom_f)
                intersection = intersection(geom_a, geom_f)
             """ """
                nuova_portata = { a.fieldNameIndex(self.text_river_q): f[self.text_source_q]}
                nuova_temp = { a.fieldNameIndex(self.text_river_t): f[self.text_source_t]}
                self.river.dataProvider().changeAttributeValues({a.id(): nuova_portata})
                self.river.dataProvider().changeAttributeValues({a.id(): nuova_temp})
                break #only one or less intersection are possible
             """ """
            end
        end
    end


 """"""
    param_point_intersect_all={ 'INPUT_FIELDS' : [], 'OUTPUT' : 'memory:', 'INTERSECT' : self.river, 'INTERSECT_FIELDS' : [], 'INPUT' : self.river }
    point_intersect_all_result=processing.run('native:lineintersections',param_point_intersect_all)
    point_intersect_all=point_intersect_all_result['OUTPUT']

    point_intersect_result=processing.run('qgis:deleteduplicategeometries',{ 'OUTPUT' : 'memory:', 'INPUT' : point_intersect_all })
    point_intersect=point_intersect_result['OUTPUT']

    features = point_intersect.getFeatures()
 """"""




    count=1
    qmain=0
    diz={}

    for f in features
        x = agd.getx(f, 0)
        y = agd.gety(f, 0)
        z = self.dem.dataProvider().identify(QgsPointXY(x, y),QgsRaster.IdentifyFormatValue)
        zresult=z.results()
        diz[zresult[1]]=[f[self.text_river_q],f[self.text_river_t],f[self.text_river_q+"_2"],f[self.text_river_q+"_2"]]



    for key in keys(diz)
        if count == 1
            if diz[key][0] >= diz[key][2]
                qmain = diz[key][0]
                q = diz[key][0]
                t = diz[key][1]
            else
                qmain = diz[key][2]
                q = diz[key][2]
                t = diz[key][3]
            end
        end
        if diz[key][0] == qmain
            q1 = q
            t1 = t
            q2 = diz[key][2]
            t2 = diz[key][3]
        else
            q2 = q
            t2 = t
            q1 = diz[key][0]
            t1 = diz[key][1]
        end

        q = q1 + q2
        t = ( ( q1 * t1 ) + ( q2 * t2 ) ) / ( q1 + q2 )
        count += 1
    end



    #    tempoanalisi = time.time() - start_time
end


end # module