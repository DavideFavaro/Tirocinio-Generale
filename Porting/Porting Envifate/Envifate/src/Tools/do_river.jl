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
  River( ma, t, x, dl, v, w, k ) = new(ma,t,x,dl,v,w,k)
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
function run_river( dem, slope, river, source, start_time, end_time, time_interval, resolution::Integer, concentration::Real, mean_hydraulic_radius::Real, fickian_x::Real=0.05,
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

    refsys = agd.getspatialref(source)

    if agd.getspatialref(river) != refsys || agd.importWKT(agd.getproj(slope)) != refsys || agd.importWKT(agd.getproj(dem)) != refsys
        throw(DomainError("The reference systems are not uniform. Aborting analysis."))
    end

    if start_time >= end_time
        throw(DomainError("`end_time` must be greater than `start_time`"))
    end


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

    start_sec, end_sec, int_sec = 60( time_start, time_end, time_interval )
    #calcolo ciclo intervallo temporale di analisi
    cicles = ( (end_time - start_time) / time_interval ) + 1
        
 # CREDO SIA L'EQUIVALENTE DI:
    #   feature = next(self.river.getFeatures())
    #   geomfeature = feature.geometry()
    #   features = self.river.getFeatures()
    features = agd.getgeom.(collect(agd.getfeature(river)))

    src_feature = collect(agd.getfeature(source))
    src_geom = agd.getgeom(src_feature[1])
    x_source = agd.getx(src_geom, 0)
    y_source = agd.gety(src_geom, 0)
    r_source, c_source = toIndexes(dem, x_source, y_source)

 # DEVO TROVARE L'EQUIVALENTE DI INTERPOLATE
    firstpoint = features[1].interpolate(0)

    x_first = agd.getx(firstpoint, 0)
    y_first = agd.gety(firstpoint, 0)
    r_first ,c_first = toIndexes(dem, x_first, y_first) 
 
    slopeband = agd.getband(slope, 1)

    demband = agd.getband(dem, 1)
    demfirstpoint = demband[r_first, c_first]
    demsource = demband[r_source, c_source]

    trend = 0
    old_x = x_first
    old_r = r_first
    old_y = y_first  
    old_c = c_first  
    if demfirstpoint >= demsource
        trend = 1
        old_x = x_source
        old_r = r_source
        old_y = y_source
        old_c = c_source       
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

        sec_cicli = start_sec
        checknumcampo = 0
        while sec_cicli <= end_sec  
            checknumcampo += 1
            nomecampo = "conc $(sec_cicli/60)"              
            prfield1 = pr.addAttributes( [ QgsField(nomecampo, QVariant.Double) ] )
            prlfield1 = prline.addAttributes( [ QgsField(nomecampo, QVariant.Double) ] )
            sec_cicli = sec_cicli + int_sec
        end


        vl.updateFields()
        vline.updateFields()
        fetline = QgsFeature()
        fetline.setGeometry( QgsGeometry.fromPolyline( [s_geom] ))
        prline.addFeatures( [ fetline ] )
        geomline = fetline.geometry()



        controllo = 1
        if trend == 1
            controllo = 0
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
                    realdistance = realdistance + resolution
                    count_index += 1
                    z = slopeband[x, y]
                    v_inst = ( mean_hydraulic_radius^(2/3) * √(z/100) ) * manning_coeff
                    push!( list_vmedia, v_inst )
                    mean_v = sum(list_vmedia) / count_index

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
                    for (i,t) in enumerate(start_sec:int_sec:end_sec)
                        element = River( concentration, t, realdistance, fickian_x, mean_v, hydraulic_section, decay_coeff )                    
                        Cfinal = calc_concentration!(element)
                        if i == 1
                            push!( list_result, Cfinal )
                        end
                     """ NON SO COSA STIA FACENDO
                        fet.setAttribute(1+ciclo,Cfinal)
                        fetline.setAttribute(1+ciclo,Cfinal)
                     """   
                    end
                    
                    vl.updateFeature(fet)
                    fet.setGeometry(point)

        
                    fetline.setGeometry( QgsGeometry.fromPolyline( [QgsPoint(old_x,old_y),QgsPoint(x,y)] ))
                    
                    vline.updateFeature(fetline)

                    feats.append(fet)
                    featlines.append(fetline)

                    old_x = x
                    old_y = y
                    avanzamento = 0
                end
                avanzamento += 1
            end
            currentdistance += 1
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