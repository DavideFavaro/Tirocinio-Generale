module Rivers



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
   os.system("start "+os.path.dirname(__file__)+"/../tutorial/manuale_envifate_dispersione_fluviale.pdf")
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

    src_geom = agd.getgeom(collect(agd.getlayer(source, 0))[1])

    if agd.geomdim(geom) != 0
        throw(DomainError(source, "`source` must be a point"))
    end

   river_layer = agd.getlayer(river, 0)

 # IL CONTROLLO SI PUO' FARE SOLO SULLE SINGOLE FEATURES
   if agd.geomdim(river) != 1
      throw(DomainError(source, "`river` must be a line"))
   end

   refsys = agd.getspatialref(geom)

   if agd.getspatialref(river_layer) != refsys || agd.importWKT(agd.getproj(slope)) != refsys || agd.importWKT(agd.getproj(dem)) != refsys
      throw(DomainError("The reference systems are not uniform. Aborting analysis."))
   end

   if start_time >= end_time
      throw(DomainError("`end_time` must be greater than `start_time`"))
   end


 # messaggio+='ALGORITMO UTILIZZATO: Fickian Mixing Process (Hemond, Harold F., and Elizabeth J. Fechner. Chemical fate and transport in the environment. Elsevier, 2014.)\n\n'

    start_sec, end_sec, int_sec = 60 .* ( time_start, time_end, time_interval )
    #calcolo ciclo intervallo temporale di analisi
    cicles = ( (end_time - start_time) / time_interval ) + 1
        
 # CREDO SIA L'EQUIVALENTE DI:
    #   feature = next(self.river.getFeatures())
    #   geomfeature = feature.geometry()
    #   features = self.river.getFeatures()
    features = agd.getgeom.(collect(layer))
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
 
    for geom in features
        length = agd.geomlength(geom)
        avanzamento = 1
        realdistance = 0
        feats = [] 
        featlines=[] 
        list_result=[]
 
        #   start_time = time.time()  
        
        if outputname == ""
            outputname = "concentrazione"
        end

        
     """ CREA UN LAYER
        vl = QgsVectorLayer("Point?crs=EPSG:"+self.refsys,self.outputname, "memory")       
        pr = vl.dataProvider()  
        prfield = pr.addAttributes( [ QgsField("distance", QVariant.Int) ] )            
        prfield2 = pr.addAttributes( [ QgsField("vmedia", QVariant.Double) ] )
     """
     # NON FUNZIONA
        vl = agd.createlayer( outputname, agd.wkbPoint, refsys )


        list_vmedia = []



     """ CREA UN'ALTRO LAYER
        vline = QgsVectorLayer("LineString?crs=EPSG:"+self.refsys, self.outputname, "memory")
        prline = vline.dataProvider()
        prlfield=prline.addAttributes( [ QgsField("distance", QVariant.Int) ] )            
        prlfield2=prline.addAttributes( [ QgsField("vmedia", QVariant.Double) ] )
     """
     # NON FUNZIONA
        vline = agd.createlayer( outputname, agd.wkbLineString, refsys )



        sec_cicli = start_sec
        checknumcampo = 0
        for sec_cicli in start_sec:int_sec:end_sec 
            checknumcampo += 1
            fieldname = "conc $(sec_cicli/60)"
            
         """ AGGIUNGE CAMPI AI LAYER CREATI
            prfield1 = pr.addAttributes( [ QgsField(nomecampo, QVariant.Double) ] )
            prlfield1 = prline.addAttributes( [ QgsField(nomecampo, QVariant.Double) ] )
         """
            agd.addfielddefn!(vl, fieldname, agd.OFTReal )
            agd.addfielddefn!(vline, fieldname, agd.OFTReal )
        end


     """ CREA UNA NUOVA FEATURE E LA AGGIUNGE A vline
        fetline = QgsFeature()
        fetline.setGeometry( QgsGeometry.fromPolyline( [s_geom] ))
        prline.addFeatures( [ fetline ] )
        geomline = fetline.geometry()
     """
        fetline = agd.createfeature( x -> x, vline )
        agd.setgeom!(fetline, agd.wkbMultiLineString)
        agd.addfeature!(vline, fetline)
        geomline = agd.getgeom(fetline)



        controllo = 1
        if trend == 1
            controllo = 0
        end

        for currentdistance in 1:length
            #   point = geom.interpolate(currentdistance)
            point = agd.getpoint(geom, currentdistance)
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

                 """ STA INIZIALIZZANDO DELLE FEATURES
                    fet = QgsFeature()
                    fet.initAttributes(2+cicli)

                    fetline = QgsFeature()
                    fetline.initAttributes(2+cicli)

                    fet.setAttribute(0,realdistance)
                    fet.setAttribute(1,self.vmedia)

                    fetline.setAttribute(0,realdistance)
                    fetline.setAttribute(1,self.vmedia) 
                 """
                    fet = agd.createfeature(x->x, vl)
                    agd.fillunsetwithdefault!(fet)
                    agd.setfield!( Ref(fet), [0, 1], [realdistance, vmedia])

                    fetline = agd.createfeature(x->x, vline)
                    agd.fillunsetwithdefault!(fetline)
                    agd.setfield!( Ref(fetline), [0, 1], [realdistance, vmedia])

                    for (cicle, t) in enumerate(start_sec:int_sec:end_sec)
                        element = River( concentration, t, realdistance, fickian_x, mean_v, hydraulic_section, decay_coeff )                    
                        finalC = calc_concentration!(element)
                        if cicle == 1
                            push!( list_result, Cfinal )
                        end

                     """ AGGIORNA I VALORI DELLE FEATURES
                        fet.setAttribute(1+ciclo,Cfinal)
                        fetline.setAttribute(1+ciclo,Cfinal)
                     """
                        agd.setfield!(fet, cicle+1, finalC )
                        agd.setfield!(fetline, cicle+1, finalC )
                    end
                    
                 """ AGGIORNA I VALORI DELLA FEATURE(?) E NE CAMBIA LA GEOMETRIA 
                    vl.updateFeature(fet)
                    fet.setGeometry(point)
                 """
                    agd.setgeom!(fet, point)


                 """ COME SOPRA
                    fetline.setGeometry( QgsGeometry.fromPolyline( [QgsPoint(old_x,old_y),QgsPoint(x,y)] ))
                    vline.updateFeature(fetline)
                 """
                    line = agd.createlinestring( Flloat64.([old_x, old_y]), Flloat64.([x, y]) )
                    agd.setgeom!(fetline, line)

                    push!(feats, fet)
                    push!(featlines, fetline)

                    old_x = x
                    old_y = y
                    avanzamento = 0
                end
                avanzamento += 1
            end
        end    

        agd.addfeature!.(Ref(vl), feats)
        agd.addfeature!.(Ref(vline), featlines)

     """ AGGIUNGE I VETTORI APPENA CREATI
        QgsProject.instance().addMapLayer(vl)
        QgsProject.instance().addMapLayer(vline)
     """
     """ PIRNT DI COSE
        tempoanalisi=time.time() - start_time
        tempostimato=time.strftime("%H:%M:%S", time.gmtime(tempoanalisi))
        messaggio="---------------------------------\n"
        messaggio+="Fine modellazione\n"
        messaggio+="\nTempo di analisi: "+tempostimato+"\n"
        messaggio+="---------------------------------\n\n"
        self.console.appendPlainText(messaggio)             
     """
     """ NON SO COSA FACCIA QUESTA PARTE
        self.tab_2.setEnabled(True)

        self.run_temp_river()
     """
    end
end

end # module