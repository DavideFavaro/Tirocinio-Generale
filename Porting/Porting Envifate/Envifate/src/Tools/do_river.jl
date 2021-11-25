module Rivers

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


import ArchGDAL as agd
using ArgParse
using Dates

include("../Library/Functions.jl")



mutable struct River
  """docstring for element"""

  ma::Real
  t::Real
  x::Real
  dl::Real
  v::Real
  w::Real
  k::Real

  concentration

  # element=river(args.concentration,args.time,args.distance,args.fickian,args.velocity)
  River( ma, t, x, dl, v, w, k ) = new( float(ma), float(t), float(x), float(dl), float(v), float(w), float(k) )
end

function calc_concentration!( r::River )   
  c1 = r.x - (r.v * r.t)
  c1_1 = -(c1^2)

  c2 = c1_1 / ( 4 * r.dl * r.t )
  c2_1 = exp( -r.k * r.t )
  c3 = exp(c2) * c2_1

  c4 = ( r.ma / r.w ) / ( √( 4π * r.dl * r.t ) )

  r.concentration = c4 * c3
  return r.concentration
end

#=
    def help(self):         
        #self.credits = u"Università della Tuscia\n Viterbo - Italy\nRaffaele Pelorosso, Federica Gobattoni\nDeveloper: Francesco Geri"
        #QMessageBox.about(self.dlg,"Credits", self.credits ) 
        if platform.uname()[0]=="Windows":
            os.system("start "+os.path.dirname(__file__)+"/../tutorial/manuale_envifate_dispersione_fluviale.pdf")
        if platform.uname()[0]=="Linux":
            os.system("xdg-open "+os.path.dirname(__file__)+"/../tutorial/manuale_envifate_dispersione_fluviale.pdf")
        else:
            os.system("open "+os.path.dirname(__file__)+"/../tutorial/manuale_envifate_dispersione_fluviale.pdf")
=#


         #                                     min         min_end   min_int        C / conc             radius
function run_river( dem, river, slope, source, start_time, end_time, time_interval, resolution::Integer, concentration::Real, mean_hydraulic_radius::Real, fickian_x::Real=0.05,
                  # w / sez                     k                    manning
                    hydraulic_section::Real=1.0, decay_coeff::Real=0, manning_coeff::Real=0.05, output_path::AbstractString )
    
 """ CONTROLLO SUI RASTER `dem` E `slope`
    if not self.slope.isValid():
        QMessageBox.warning(self,"Warning", u"La mappa delle pendenza è mancante o non valida" )
        return False

    if not self.dem.isValid():
        QMessageBox.warning(self,"Warning", u"Il DEM è mancante o non valido" )
        return False
 """
    if agd.geomdim(source) != 0
        throw(DomainError(source, "`source` must be a point"))
    end

 # NON SONO SICURO CHE IL CONTROLLO SIA EQUIVALENTE A: `self.river.wkbType()!=5`
    if agd.geomdim(river) != 1
        throw(DomainError(source, "`river` must be a line"))
    end

    if agd.getspatialref(river) != agd.getspatialref(source) || agd.getspatialref(slope) != agd.getspatialref(source) || agd.getspatialref(dem) != agd.getspatialref(source)
        throw(DomainError("The reference systems are not uniform. Aborting analysis."))
    end

    if start_time >= end_time
        throw(DomainError("`end_time` must be greater than `start_time`"))
    end

    refsys = agd.importEPSG(agd.fromWKT(agd.getspatialref(source)))

 """ PRINT DI COSE
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
    messaggio+="Coefficiente decadimento lamda: "+str(self.decay_coeff)+"\n"
    messaggio+="Coefficiente di scabrezza: "+str(self.manning_coeff)+"\n"
    messaggio+="Massa inquinante: "+str(self.concentration)+"\n"
    messaggio+="Tempo iniziale dell'analisi: "+str(self.start_time)+"\n"
    messaggio+="Tempo finale dell'iniezione: "+str(self.min_end)+"\n"
    messaggio+="Intervallo temporale: "+str(self.min_int)+"\n"
    messaggio+="Risoluzione: "+str(self.res)+"\n\n"
    messaggio+='ALGORITMO UTILIZZATO: Fickian Mixing Process (Hemond, Harold F., and Elizabeth J. Fechner. Chemical fate and transport in the environment. Elsevier, 2014.)\n\n'
    messaggio+="---------------------------\n\n"
    self.console.appendPlainText(messaggio)
 """

    #calcolo ciclo intervallo temporale di analisi
    cicli = ( (end_time - start_time) / time_interval ) + 1
        
 # CREDO SIA L'EQUIVALENTE DI:
    #   feature = next(self.river.getFeatures())
    #   geomfeature = feature.geometry()
    #   features = self.river.getFeatures()
    features = getgeom.(collect(agd.getfeature(river)))

    src_feature = collect(agd.getfeature(source))
    src_geom = agd.getgeom(src_feature[1])
    x_source = agd.getx(src_geom, 0)
    y_source = agd.gety(src_geom, 0)

 # DEVO TROVARE L'EQUIVALENTE DI INTERPOLATE
    firstpoint = features[1].interpolate(0)

    x_first = agd.getx(firstpoint, 0)
    y_first = agd.gety(firstpoint, 0)
 
 # POTREBBE NON ESSERE COSI' SEMPLICE
    demfirstpoint = dem[x_first, y_first]
    demsource = dem[x_source, y_source]

    if demfirstpoint >= demsource
        trend = 1
        old_x = x_source
        old_y = y_source
    else
        trend = 0
        old_x = x_first
        old_y = y_first           
    end
 
    for f in features
        length = agd.geomlength(geom)

        currentdistance = 1
        avanzamento = 1
        realdistance = 0
        feats = [] 
        featlines=[] 
        list_result=[]
 
        #   start_time = time.time()  

     # NON SO SE VADA ELIMINATO
        if self.outputname == ""
            self.outputname = "concentrazione"
        end

        vl = QgsVectorLayer("Point?crs=EPSG:"+self.refsys,self.outputname, "memory")       
        pr = vl.dataProvider()  
        prfield = pr.addAttributes( [ QgsField("distance", QVariant.Int) ] )            
        prfield2 = pr.addAttributes( [ QgsField("vmedia", QVariant.Double) ] )


        list_vmedia=[]

        vline = QgsVectorLayer("LineString?crs=EPSG:"+self.refsys, self.outputname, "memory")
        prline = vline.dataProvider()
        prlfield=prline.addAttributes( [ QgsField("distance", QVariant.Int) ] )            
        prlfield2=prline.addAttributes( [ QgsField("vmedia", QVariant.Double) ] )

        sec_cicli = sec
        checknumcampo = 0
        while sec_cicli <= sec_end  
            checknumcampo += 1
            nomecampo = "conc $(sec_cicli/60)"              
            prfield1 = pr.addAttributes( [ QgsField(nomecampo, QVariant.Double) ] )
            prlfield1 = prline.addAttributes( [ QgsField(nomecampo, QVariant.Double) ] )
            sec_cicli = sec_cicli + sec_int
        end


        vl.updateFields()
        vline.updateFields()
        fetline = QgsFeature()
        fetline.setGeometry( QgsGeometry.fromPolyline( [s_geom] ))
        prline.addFeatures( [ fetline ] )
        geomline = fetline.geometry()

        if trend == 1
            controllo = 0
        else
            controllo=1
        end

        while currentdistance < length
            
            point = geom.interpolate(currentdistance)
            x = agd.getx(point, 0)
            y = agd.gety(point, 0)

            dist = √( (x - x_source)^2 + (y - y_source)^2 )
            if trend == 1 && dist <= 1
                controllo = 1
            end
            if trend == 0 && dist <= 1
                controllo = 0
            end

            if controllo == 1
                if avanzamento == resolution
                    realdistance = realdistance + ressolution
                    count_index += 1
                    z = slope[x, y][1]
                    v_inst = ( mean_hydraulic_radius^(2/3) * √(z/100) ) * manning_coeff
                    push!( list_vmedia, v_inst )
                    vmadia = sum(list_vmedia) / count_index

                 """ NON SO COSA RAPPRESENTINO QUESTE RIGHE
                    fet = QgsFeature()
                    fet.initAttributes(2+cicli)

                    fetline = QgsFeature()
                    fetline.initAttributes(2+cicli)

                    fet.setAttribute(0,realdistance)
                    fet.setAttribute(1,self.vmedia)

                    fetline.setAttribute(0,realdistance)
                    fetline.setAttribute(1,self.vmedia) 
                 """

                    ciclo = 1

                    while sec <= sec_end                            
                        element = River( concentration, sec, realdistance, fickian_x, vmedia, hydraulic_section, decay_coeff )                    
                        Cfinal = calc_concentration!(element)
                        if ciclo == 1
                            push!( list_result, Cfinal )
                        end
                     """ NON SO COSA STIA FACENDO
                        fet.setAttribute(1+ciclo,Cfinal)
                        fetline.setAttribute(1+ciclo,Cfinal)
                     """
                        ciclo = ciclo + 1
                        sec = sec + sec_int
                    end
                    
                    vl.updateFeature(fet)
                    fet.setGeometry(point)

        
                    fetline.setGeometry( QgsGeometry.fromPolyline( [QgsPoint(old_x,old_y),QgsPoint(x,y)] ))
                    
                    vline.updateFeature(fetline)

                    feats.append(fet)
                    featlines.append(fetline)

                    old_x=x
                    old_y=y

                    avanzamento=0
                end
                avanzamento=avanzamento+1
            end

            currentdistance = currentdistance + 1
        end    


        pr.addFeatures(feats)
        prline.addFeatures(featlines)
        vl.updateFields()
        vline.updateFields()

        QgsProject.instance().addMapLayer(vl)
        QgsProject.instance().addMapLayer(vline)

     """ PIRNT DI COSE
        tempoanalisi=time.time() - start_time
        tempostimato=time.strftime("%H:%M:%S", time.gmtime(tempoanalisi))
        messaggio="---------------------------------\n"
        messaggio+="Fine modellazione\n"
        messaggio+="\nTempo di analisi: "+tempostimato+"\n"
        messaggio+="---------------------------------\n\n"
        self.console.appendPlainText(messaggio)             
     """
        self.tab_2.setEnabled(True)

        self.run_temp_river()
    end
end

end # module