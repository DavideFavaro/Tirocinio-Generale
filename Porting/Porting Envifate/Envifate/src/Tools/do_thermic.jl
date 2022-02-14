module Thermic

#=
/tutorial/manuale_envifate_ruscellamento.pdf")
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



"""
    run_thermic( dem::AbstractArray, source, river, source_temperature::Real, source_flow_rate::Real, river_temperature::Real, river_flow_rate::Real, output_path::AbstractString=".\\")


"""

function run_thermic( dem::AbstractArray, source, river, source_temperature::Real, source_flow_rate::Real, river_temperature::Real, river_flow_rate::Real, output_path::AbstractString=".\\")

    src_geom = agd.getgeom(collect(agd.getlayer(source, 0))[1])

    if agd.geomdim(src_geom) != 0
        throw(DomainError(source, "`source` must be a point"))
    end

    river_layer = agd.getlayer(river, 0)

 # IL CONTROLLO SULLA GEOMETRIA PUO' ESSERE FATTO SOLO SU FEATURES
    if agd.geomdim(river_layer) != 1
        throw(DomainError(source, "`river` must be a line"))
    end

    layer = agd.getlayer(river, 0)
    refsys = agd.getspatialref(source)

    if agd.importWKT(agd.getproj(dem)) !=  refsys ||  agd.getspatialref(layer) != refsys
        throw(DomainError("The reference systems are not uniform. Aborting analysis." ))
    end

    demband = agd.getband(dem, 1)

    #   start_time = time.time()
    
 """ CREDO STIA CONSIDERANDO UN'AREA DI 5 METRI(?) INTORNO ALLA SORGENTE
    buffer_sorgente = processing.run("native:buffer", {"INPUT": self.source, "DISTANCE": 5, "OUTPUT":"memory:"})
    sorgente = buffer_sorgente["OUTPUT"]
 """
    # Consider a 5m area around the source
    source_area = agd.buffer(source, 5.0)
    # Find the portion of the river that intersects said area; there could be no intersection
    for feature in layer
        segment = agd.getgeom(feature)
        if agd.intersects(source_area, segment)
            intersection = agd.intersection(source_area, segment)

         """ IN QUESTA PARTE SETTA I VALORI DI TEMPERATURA E PORTATA DEL FIUME A QUELLI DELLA FUNZIONANTE
            nuova_portata = { a.fieldNameIndex(self.text_river_q): f[self.text_source_q]}
            nuova_temp = { a.fieldNameIndex(self.text_river_t): f[self.text_source_t]}
            self.river.dataProvider().changeAttributeValues({a.id(): nuova_portata})
            self.river.dataProvider().changeAttributeValues({a.id(): nuova_temp})
         """
            # BISOGNA VEDERE COME GESTIRE L'INPUT DI source E river
            new_flow = Dict("Portata" => source_flow_rate) 
            new_temperature = Dict("Temperatura" => source_temperature)
            flow_idx, temp_idx = agd.findfieldindex( Ref(feature), [:portata, :temperatura] )
            agd.setfield!.( Ref(feature), [flow_idx, temp_idx], [new_flow, new_temperature] )
            break
        end
    end


 """ STA PRENDENDO L'INTERSEZIONE CON river DI QUALCOSA, MA NON SO COSA """
 # FORSE STA PRENDENDO LE CELLE DAL RASTER CHE CORRISPONDONO A river
    param_point_intersect_all = Dict( :Input_Fields => [], :Output => :memory, :Intersect => river, :Intersect_Fields => [], :Input => river )
    point_intersect_all_result = processing.run( "native:lineintersections", param_point_intersect_all )
    point_intersect_all = point_intersect_all_result[:Output]

    point_intersect_result = processing.run( "qgis:deleteduplicategeometries", Dict( :Output => :memory, :Input => point_intersect_all ) )
    point_intersect = point_intersect_result[:Output]
 """                                                                 """
    # NON E' CORRETTO, BISOGNA VEDERE COME OTTENIAMO point_intersect
    features = agd.getfeatures(point_intersect)




    count = 1
    frmain = 0
    dict = Dict()
    for feature in features
        geom = agd.getgeom(feature)
        x = agd.getx(geom, 0)
        y = agd.gety(geom, 0)
        r, c = toIndexes(dtm, x, y) 
        z = demband[r, c]

        res = agd.getfield.( Ref(feature), [:portata_river, :temperatura_river, :portata_river_2, :temperatura_river_2] )
        if z in keys(dict)
            dict[z] = res
        else
            push!( dict, z => res )
        end
    end

    flow_rate = 0.0
    temperature = 0.0
    for key in keys(dict)
        if count == 1
            if dict[key][1] >= dict[key][3]
                frmain = dict[key][1]
                flow_rate = dict[key][1]
                temperature = dict[key][2]
            else
                frmain = dict[key][3]
                flow_rate = diz[key][3]
                temperature = diz[key][4]
            end
        end
        if dict[key][1] == frmain
            flow_rate1 = flow_rate
            temperature1 = temperature
            flow_rate2 = dict[key][3]
            temperature2 = dict[key][4]
        else
            flow_rate1 = flow_rate
            temperature1 = temperature
            flow_rate2 = dict[key][1]
            temperature2 = dict[key][2]
        end

        flow_rate = flow_rate1 + flow_rate2
        temperature = ( ( flow_rate1 * temperature1 ) + ( flow_rate2 * temperature2 ) ) / ( flow_rate1 + flow_rate2 )
        count += 1
    end



    #    tempoanalisi = time.time() - start_time
end



end # module