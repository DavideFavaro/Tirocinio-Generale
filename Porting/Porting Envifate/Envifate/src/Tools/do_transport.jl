module Transport

#=
    def help(self):         
        #self.credits = u"Università della Tuscia\n Viterbo - Italy\nRaffaele Pelorosso, Federica Gobattoni\nDeveloper: Francesco Geri"
        #QMessageBox.about(self.dlg,"Credits", self.credits ) 
        if platform.uname()[0]=="Windows":
            os.system("start "+os.path.dirname(__file__)+"/../tutorial/manuale_envifate_solute.pdf")
=#



using ArchGDAL


const agd = ArchGDAL


                      # top          bottom          q                      phead
function run_transport( aquifer_top, aquifer_bottom, concentration_sources, piezomatic_head, status, source, srs, tensor_x::Real, tensor_y::Real, porosity::Real, dispersiveness_t::Real,
                        dispersiveness_l::Real, storativity::Real=0.0001, break_error::Real=0.000001, delay::Real=1.0, time::Integer, resolution::Integer, output_path::AbstractString=".\\" )


 #= E' UN RASTER MA LO PONE A ZERO
    messaggio_q=''
    if self.text_q!='No file':
        self.q=self.listalayers[self.text_q]
        messaggio_q="Mappa concentration sources and skins: "+str(self.q)+"\n"
        if not self.q.isValid():
            QMessageBox.warning(self,"Warning", u"Il file raster concentration sources and skins non è valido" )
            return
    else:
        self.q=0.0
 =#

  #= COSA SONO QUESTI PARAMETRI
    self.text_campophead = str(self.combo_campophead.currentText())
    self.text_campostatus = str(self.combo_campostatus.currentText())
    self.text_solver = str(self.combo_solver.currentText())
    self.text_conc = str(self.combo_conc.currentText())
 =#

    src_geom = agd.getgeom(collect(agd.getlayer(source, 0))[1])

    if agd.geomdim(source) != 0
        throw(DomainError(source, "`source` must be a point"))
    end
 

   time *= 3600 


 """ PRINT DI COSE
    messaggio+=messaggio_q
    messaggio+=messaggio_mappatop
    messaggio+=messaggio_mappabottom
    # messaggio+=messaggio_status
    messaggio+="VARIABILI:\n"
    messaggio+=u"Concentrazione inquinante: da shapefile"
    messaggio+=u"Solver: "+str(self.combo_solver.currentText())+"\n"
    messaggio+=messaggio_top
    messaggio+=messaggio_bottom
    messaggio+='ALGORITMO UTILIZZATO: modello realizzato sfruttando r.gwflow e r.solute.transport di Grass GIS v.7 (Neteler and Mitasova, 2008). Per i dettagli matematici degli algoritmi vedere https://grass.osgeo.org/gdp/hydrology/gebbert2007_diplom_stroemung_grass_gis.pdf\n\n'
 """

   # start_time = time.time()


   self.areastudio_path=self.areastudio.dataProvider().dataSourceUri().split('|')
   self.phead_path=self.phead.dataProvider().dataSourceUri().split('|')
   self.source_path=self.source.dataProvider().dataSourceUri().split('|')
   self.status_path=self.status.dataProvider().dataSourceUri().split('|')


   #g.proj -c epsg=3003
   run_command("g.proj", epsg=self.myepsg, flags = 'c')

   # Carica e trasforma in raster "area"
    # Inporta "areastudio" come vettoriale
   run_command("v.in.ogr", input=self.areastudio_path[0], output="areastudio",overwrite=True, flags = 'o')
    # Definisce la regione di "areastudio"
   run_command("g.region", res=self.res, vector="areastudio")
    # Trasforma "areastudio" in raster
   run_command("v.to.rast", input="areastudio", output="areastudio",overwrite=True, use = 'val')





 """ Carica e trasforma in raster il vettore di "top"
  # Inporta "phead" come vettoriale
    run_command("v.in.ogr", input=self.phead_path[0], output="phead",overwrite=True, flags = 'o')
  # Trasforma "phead" in raster
    run_command("v.to.rast", input="phead", output="phead",overwrite=True, use = 'attr', attribute_column=self.text_campophead)
 """
   pheadband = agd.getband(piezomatic_head, 1)

 """ Carica e trasforma in raster "status" 
    run_command("v.in.ogr", input=self.status_path[0], output="status",overwrite=True, flags = 'o')
    run_command("v.to.rast", input="status", output="status",overwrite=True, use = 'attr', attribute_column=self.text_campostatus)
 """
   statusband = agd.getband(status, 1)

 """ Carica source lo trasforma in raster e setta la cella corrispondente a 0/null
  # Inporta "source" come vettoriale
    run_command("v.in.ogr", input=self.source_path[0], output="source",overwrite=True, flags = 'o')
  # Trasforma "source" in raster
    run_command("v.to.rast", input="source", output="source",overwrite=True, use = 'attr', attribute_column=self.text_conc)
  # Rimpiazza i null value in "source" con 0
    run_command("r.null", map="source", null=0)
 """
   x_source = agd.getx(src_geom, 0)
   y_source = agd.gety(src_geom, 0)
   r_source, c_source = toIndexes(dtm, x_source, y_source)


 """ CREDO STIA SOLO SETTANDO DEI VALORI
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
 """


   # create solv les
    # setup 
   res = cg( phead, status )



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



    #   tempoanalisi=time.time() - start_time
     
end

end # module