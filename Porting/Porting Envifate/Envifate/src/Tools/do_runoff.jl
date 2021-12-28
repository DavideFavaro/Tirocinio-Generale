module Runoffs

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



function getFeatureByFid( features::Vector{agd.Feature}, fid )
    for f in features
        if agd.getfid(f) == fid
            return f
        end
    end
    return nothing
end


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


        d.exec_()

    def popolacombo(self):
        self.combo_source.clear()
        self.combo_dem.clear()
        self.combo_soil.clear()
        self.combo_bound.clear()
        self.combofield_lc.clear()
        self.combo_target.clear()
        self.combofield_soil.clear()
        self.combofield_target.clear()
        self.combo_fieldp.clear()
        self.combo_lc.clear()
        self.line_folder.clear()
        self.progressBar.setValue(0)


        self.allLayers = self.canvas.layers()
        self.listalayers=dict()
        #elementovuoto="No required"
        for i in self.allLayers:
            if i.type() == QgsMapLayer.VectorLayer:
                self.listalayers[i.name()]=i
                self.combo_source.addItem(str(i.name()))
                self.combo_bound.addItem(str(i.name()))
                self.combo_lc.addItem(str(i.name()))
                self.combo_target.addItem(str(i.name()))
            if i.type()==QgsMapLayer.RasterLayer:
                self.listalayers[i.name()]=i
                self.combo_dem.addItem(str(i.name()))

        self.combo_soil.addItem("Valore campo")
        self.combo_soil.addItem("A")
        self.combo_soil.addItem("B")
        self.combo_soil.addItem("C")
        self.combo_soil.addItem("D")

        self.popolafields(self.combo_lc,self.combofield_lc)
        self.popolafields(self.combo_lc,self.combofield_soil)
        self.popolafields(self.combo_source,self.combo_fieldp)
        self.popolafields(self.combo_target,self.combofield_target)


    def extract_values(self, raster,x,y):
        z=raster.dataProvider().identify(QgsPointXY(x, y),QgsRaster.IdentifyFormatValue)
        zresult=z.results()
        zvalue=zresult[1]
        return(zvalue)
=#
#=
    def run_runoff(self):
        self.text_vector = str(self.combo_source.currentText())
        self.text_area = str(self.combo_bound.currentText())
        self.text_dem = str(self.combo_dem.currentText())
        self.text_lc = str(self.combo_lc.currentText())
        self.text_lcfield = str(self.combofield_lc.currentText())
        self.text_p = str(self.combo_fieldp.currentText())
        self.text_target = str(self.combo_target.currentText())
        self.text_targetfield = str(self.combofield_target.currentText())
        self.text_lcfield = str(self.combofield_lc.currentText())
        self.text_soil = str(self.combo_soil.currentText())
        self.text_soilfield = str(self.combofield_soil.currentText())

        self.res=int(self.spinRes.text())






        # wkbType: 1:point, 6:multipolygon, 2: Linestring

        self.dem=self.listalayers[self.text_dem]

        if not self.dem.isValid():
            QMessageBox.warning(self,"Warning", "The dem file is not valid" )
            return

        self.source=self.listalayers[self.text_vector]

        if self.source.wkbType()!=1:
            QMessageBox.warning(self,"Warning", "The source file must have point geometry" )
            return

        self.areastudio=self.listalayers[self.text_area]


        if self.areastudio.wkbType()!=6:
            QMessageBox.warning(self,"Warning", "The boundaries file must have polygon geometry" )
            return

        self.target=self.listalayers[self.text_target]

        if self.target.wkbType()!=6:
            QMessageBox.warning(self,"Warning", "L'area target deve avere geometria poligonale" )
            return

        self.lc=self.listalayers[self.text_lc]


        if self.lc.wkbType()!=6:
            QMessageBox.warning(self,"Warning", "Not a valid landcover geometry" )
            return

        # self.path_output=self.line_output.text()
        # if self.path_output=="":
        #     self.path_output=os.path.dirname(__file__)+"/runoff.tif"

        if self.areastudio.crs().authid()!=self.source.crs().authid() or self.lc.crs().authid()!=self.source.crs().authid() or self.dem.crs().authid()!=self.source.crs().authid() or self.target.crs().authid()!=self.source.crs().authid():
            QMessageBox.warning(self,"Warning", "Errore: i sistemi di riferimento non sono uniformi. Impossibile continuare con l'analisi." )
            return


        self.refsys=self.source.crs().authid().split(':')[1]


        self.path_working=self.line_folder.text()
        if self.path_working=="":
            self.path_working=os.path.dirname(__file__)

        self.path_temp_lc=self.path_working+"/temp_lc.tif"
        self.path_temp_soil=self.path_working+"/temp_soil.tif"


        #recupero dati database

        listaclc=functions.cn_list_extract()


        controllo_soil=0

        messaggio="Inizio elaborazione analisi dispersione per ruscellamento\n"
        messaggio+="---------------------------\n\n"
        messaggio+="FILE DI INPUT:\n"
        messaggio+="Vettoriale sorgente: "+str(self.text_vector)+"\n"
        messaggio+="Vettoriale confine: "+str(self.text_area)+"\n"
        messaggio+="Vettoriale target: "+str(self.text_target)+"\n"
        messaggio+="DTM: "+str(self.text_dem)+"\n\n"

        messaggio+="VARIABILI:\n"
        messaggio+="Risoluzione: "+str(self.res)+"\n\n"
        messaggio+='ALGORITMO UTILIZZATO: calcolo della separazione delle componenti infiltrazione e ruscellamento tramite metodo SCS-CN; US Department of Agriculture Soil Conservation Service, 1972. National Engineering Handbook, Section 4, Hydrology. US Government Printing Office, Washington, DC, 544pp.\n\n'
        messaggio+="---------------------------\n\n"
        self.console.appendPlainText(messaggio)


        self.label_status.setText("Preparazione dati")
        self.label_status.setStyleSheet('color : #e8b445;font-weight:bold')

        self.path_working=self.line_folder.text()
        if self.path_working=="":
            self.path_working=os.path.dirname(__file__)

        self.path_output=self.path_working+"/runoff.tif"



        path_layer=self.areastudio.dataProvider().dataSourceUri()
        path=path_layer.split("|")
        source_ds = ogr.Open(path[0])
        area_layer = source_ds.GetLayer()
        x_min=int(area_layer.GetExtent()[0])
        y_min=int(area_layer.GetExtent()[2])
        x_max=int(area_layer.GetExtent()[1])
        y_max=int(area_layer.GetExtent()[3])

        drivermem = gdal.GetDriverByName('MEM')
        pixel_size = self.res
        NoData_value = -9999

        # Create the destination data source
        x_res = (x_max - x_min) / pixel_size
        y_res = (y_max - y_min) / pixel_size

        target_ds = gdal.GetDriverByName('GTiff').Create(self.path_output, int(x_res), int(y_res), 1, gdal.GDT_Float32)
        target_ds.SetGeoTransform((x_min, pixel_size, 0, y_max, 0, -pixel_size))
        projectionfrom = target_ds.GetProjection()

        srs = osr.SpatialReference()
        srs.ImportFromEPSG(int(self.refsys))
        target_ds.SetProjection( srs.ExportToWkt() )

        target_ds.SetMetadata({'credits':'Envifate - Francesco Geri, Oscar Cainelli, Paolo Zatelli, Gianluca Salogni, Marco Ciolli - DICAM Università degli Studi di Trento - Regione Veneto',
                               'modulo':'Analisi ruscellamento',
                               'descrizione':'Analisi di ruscellamento di un inquinante attraverso il metodo della separazione delle componenti',
                               'srs':self.source.crs().authid(),
                               'data':datetime.datetime.now().strftime("%d-%m-%y")})

        # geotransform = target_ds.GetGeoTransform()


        band = target_ds.GetRasterBand(1)
        band.SetNoDataValue(float(NoData_value))
        band.Fill(NoData_value)
        xsize = band.XSize
        ysize = band.YSize


        intervallo=int(self.dem.rasterUnitsPerPixelX())

        outData = np.array(band.ReadAsArray(0, 0, xsize,ysize).astype(np.float))

        lc_clip_proc = processing.run('qgis:clip', {'INPUT':self.lc, 'OVERLAY':self.areastudio, 'OUTPUT':self.path_working+'/clip.gpkg'})
        lc_clip=QgsVectorLayer(lc_clip_proc['OUTPUT'], 'lc_clip', 'ogr')
        lc_clip.setCrs(self.source.crs())


        #path__layer_lc=lc_clip['OUTPUT'].dataProvider().dataSourceUri()
        path__layer_lc=lc_clip.dataProvider().dataSourceUri()
        path_lc=path__layer_lc.split("|")
        source_ds_lc = ogr.Open(path_lc[0])
        lc_layer = source_ds_lc.GetLayer()
        lc_ds = gdal.GetDriverByName('GTiff').Create(self.path_temp_lc, int(x_res), int(y_res), 1, gdal.GDT_Float32)
        lc_ds.SetGeoTransform((x_min, 25, 0, y_max, 0, -25))
        lc_ds.SetProjection( srs.ExportToWkt() )
        band_lc = lc_ds.GetRasterBand(1)
        band_lc.SetNoDataValue(float(-9999))
        band_lc.Fill(-9999)
        xsize = band_lc.XSize
        ysize = band_lc.YSize
        gdal.RasterizeLayer(lc_ds, [1], lc_layer,options=["ATTRIBUTE="+self.text_lcfield])
        lc_ds=None

        lc_layer=QgsRasterLayer(self.path_temp_lc,"lc_layer")



        if self.text_soil=="Valore campo":
            controllo_soil=1

            source_ds_soil = ogr.Open(path_lc[0])
            soil_layer = source_ds_soil.GetLayer()
            soil_ds = gdal.GetDriverByName('GTiff').Create(self.path_temp_soil, int(x_res), int(y_res), 1, gdal.GDT_Float32)
            soil_ds.SetGeoTransform((x_min, 25, 0, y_max, 0, -25))
            soil_ds.SetProjection( srs.ExportToWkt() )
            band_soil = soil_ds.GetRasterBand(1)
            band_soil.SetNoDataValue(float(-9999))
            band_soil.Fill(-9999)
            gdal.RasterizeLayer(lc_ds, [1], soil_layer,options=["ATTRIBUTE="+self.text_soilfield])
            soil_ds=None

            soil_layer=QgsRasterLayer(self.path_temp_soil,"soil_layer")


        grass_area=str(x_min)+','+str(x_max)+','+str(y_min)+','+str(y_max)+' ['+str(self.areastudio.crs().authid())+']'
        # grass_coord=str(x_source)+','+str(y_source)+' ['+str(self.source.crs().authid())+']'

        namewatershed=self.path_working+'/watershed'
        namedrain=self.path_working+'/wshed.shp'

        params = { 'GRASS_RASTER_FORMAT_OPT' : '','GRASS_REGION_CELLSIZE_PARAMETER' : 0, 'GRASS_REGION_PARAMETER' :grass_area,
                   'GRASS_VECTOR_EXPORT_NOCAT' : False, '-a' : False, 'start_coordinates' : None,  '-n' : False,
                   'input' : self.dem.dataProvider().dataSourceUri(),'-c' : True, 'drain' : namedrain, 'GRASS_MIN_AREA_PARAMETER' : 0.0001,
                   'start_points' : self.source.dataProvider().dataSourceUri(), 'output' : namewatershed }



        waterwshed_proc = processing.run('grass7:r.drain', params)


        #aggiungo per controllo la viewshed alla toc
        #iface.addVectorLayer(namedrain,'watershed','ogr')
        #watershed=QgsProject.instance().mapLayersByName('watershed.shp')


        vdrain = QgsVectorLayer(namedrain, 'vdrain', 'ogr')

        idxlevel = self.source.fields().indexFromName(self.text_p)
        idxcat = vdrain.fields().indexFromName('cat')

        idxtargetname = self.target.fields().indexFromName(self.text_targetfield)



        # pyqtRemoveInputHook()
        # pdb.set_trace()

        features = vdrain.getFeatures()

        nfeat=0

        polygons_t = [feature_t for feature_t in self.target.getFeatures()]


        start_time = time.time()

        for f in features:
            geom = f.geometry()
            length = geom.length()
            currentdistance=intervallo
            nfeat+=1

            featlines=[]


            firstpoint=geom.interpolate(0)

            old_x=firstpoint.asPoint()[0]
            old_y=firstpoint.asPoint()[1]

            fileoutput=self.path_working+'/drain'+str(nfeat)+'.shp'


            max_progress=length/intervallo
            self.progressBar.setMaximum(max_progress)



            vline = QgsVectorLayer("LineString?crs=EPSG:"+self.refsys, "drain"+str(nfeat), "memory")

            prline = vline.dataProvider()
            prlfield=prline.addAttributes( [ QgsField("concentrazione", QVariant.Double) ] )

            self.label_status.setText("Processing data")
            self.label_status.setStyleSheet('color : #e8b445;font-weight:bold')

            idf=f.attributes()[idxcat]
            feat_drain = next(self.source.getFeatures(QgsFeatureRequest().setFilterFid(idf-1)))
            p00=feat_drain.attributes()[idxlevel]

            index_progress=0

            while currentdistance < length:
                if index_progress==0:
                    p0=p00
                else:
                    p0=pe


                point = geom.interpolate(currentdistance)
                x=point.asPoint()[0]
                y=point.asPoint()[1]
                clc=self.extract_values(lc_layer,x,y)
                if controllo_soil==1:
                    soil=self.extract_values(soil_layer,x,y)
                else:
                    soil=self.text_soil

                try:
                    cn=listaclc[clc][self.classisoil[soil]]
                    S=self.calc_s(int(cn))
                except Exception as e:
                    S=0

                pcheck=(math.pow((p0-(0.2*S)),2))/(p0-(0.2*S)+S)
                if pcheck>(0.2*S):
                    pe=pcheck
                    fetline = QgsFeature()
                    fetline.setGeometry( QgsGeometry.fromPolyline( [QgsPoint(old_x,old_y),QgsPoint(x,y)] ))
                    fetline.initAttributes(1)
                    fetline.setAttribute(0,pe)
                    vline.updateFeature(fetline)

                    featlines.append(fetline)
                    index_progress+=1

                    for pol_t in polygons_t:
                        poly_t = pol_t.geometry()
                        if poly_t.contains(point):
                            nometarget=pol_t.attributes()[idxtargetname]
                            messaggio="\nIl vettore drain"+str(nfeat)+" ha raggiunto l'area bersaglio denominata '"+str(nometarget)+"' con un volume pari a: "+str(round(pe,3))+" mm\n"
                            self.console.appendPlainText(messaggio)
                            currentdistance=length+1

                else:
                    pe=0
                    currentdistance=length+1
                self.progressBar.setValue(max_progress)
                # pyqtRemoveInputHook()
                # pdb.set_trace()

                old_x=x
                old_y=y

                currentdistance = currentdistance + intervallo


            prline.addFeatures(featlines)
            vline.updateFields()
            QgsProject.instance().addMapLayer(vline)




        tempoanalisi=time.time() - start_time
        tempostimato=time.strftime("%H:%M:%S", time.gmtime(tempoanalisi))
        messaggio="---------------------------------\n"
        messaggio+="Fine modellazione\n"
        messaggio+="\nTempo di analisi: "+tempostimato+"\n"
        messaggio+="---------------------------------\n\n"
        self.console.appendPlainText(messaggio)

        self.label_status.setText("In attesa di dati")
        self.label_status.setStyleSheet('color : green; font-weight:bold')
        self.progressBar.setValue(max_progress)
=#













using Plots
import ArchGDAL as agd


function flow( map::Matrix, sourceRow::Integer, sourceCol::Integer )
    source = map[sourceRow, sourceCol]
    
end


function flow!( map::AbstractArray, dem_band, row::Integer, col::Integer, noDataValue::Real, cicles::Integer )
    if cicles <= 0
        return nothing
    end
    # Limits of the raster
    maxR, maxC = size(dem_band)
    # Adjacent cells
    indexes = [
        ( dem_band[ row+i, col+j ], row+i, col+j )
        for i in -1:1
            for j in -1:1
                if ( i != 0 || j != 0 ) &&
                    ( row+i >= 1 && row+i <= maxR ) &&
                    ( col+j >= 1 && col+j <= maxC ) &&
                    dem_band[ row+i, col+j ] != noDataValue &&
                    ismissing( map[ row+i, col+j ] )
    ]
    # Adjacent cells' minimum height
    min_height =  minimum( t -> t[1], indexes )

    if height < min_height
        return nothing
    end

    for (h, r, c) in indexes
        if h == min_height
            map[r, c] = true
            flow!(map, dem_band, r, c, noDataValue, cicles-1)
        else
            map[r, c] = false
        end
    end

    return nothing
end


function flow!( stream_points::AbstractVector, visited_points::AbstractVecotr, dem_band, row::Integer, col::Integer, noDataValue::Real, cicles::Integer )
    if cicles <= 0
        return nothing
    end
    # Limits of the raster
    maxR, maxC = size(dem_band)
    # Adjacent cells
    indexes = [
        ( row+1, col   ),
        ( row+1, col+1 ),
        ( row,   col+1 ),
        ( row-1, col+1 ),
        ( row-1, col   ),
        ( row-1, col-1 ),
        ( row,   col-1 ),
        ( row+1, col-1 )
    ]
    # Valid adjacent cells
    condition(a, b) = ( a >= 1 && a <= maxR ) && ( b >= 1 && b <= maxC ) && dem_band[a, b] != noDataValue && ( (a, b) ∉ stream_points || (a, b) ∉ visited_points )
    visited = filter( (r, c) -> !condtion(r, c), indexes )
    push!(visited_points, visited...)
    filter!(condition, indexes)
    # Adjacent cells' minimum height
    min_height =  minimum( (r,c) -> dem_band[r,c], indexes )

    if height < min_height
        return nothing
    end

    for (r, c) in indexes
        if dem_band[r, c] == min_height
            push!( stream_points, (r, c) )
            flow!(stream_points, other_points, dem, r, c, noDataValue, cicles-1)
        else
            push!( visited_points, (r, c) )
        end
    end

    return nothing  
end



heatmap( 1:1000, 1:1000, test1, c=cgrad([:blue, :white, :yellow, :red]) )

mat = Array{Any}( missing, 1000, 1000 )
flow!( mat, test1, 300, 300, ndv, 10000 )
heatmap( 1:1000, 1:1000, mat, c =[:orange, :black, :blue] )



heatmap( 1:100, 1:100, test2, c=cgrad([:blue, :white, :yellow, :red]) )

mat = Array{Any}( missing, 100, 100 )
flow!( mat, test2, 60, 90, ndv, 1000 )
heatmap( 1:100, 1:100, mat, c =[:orange, :black, :blue] )

mat = Array{Any}( missing, 5, 5 )
flow!( mat, test3, 3, 3, ndv, 50 )


























# VEDERE @view, @inbound, @turbo, @fast, LoopedVectorization.jl e StableArrays.jl PER ULTERIORI OTTIMIZZAZIONI
function looseIn( tuple::Tuple{T1, T1}, tuples::Vector{ Tuple{T1, T1, T2} } ) where {T1 <: Number, T2 <: Number}
    for t in tuples
        if tuple[1] == t[1] && tuple[2] == t[2]
            return true
        end
    end
    return false
end

function connectivity_batch!( mat::AbstractArray{ Union{ Missing, Vector{Tuple{Int64, Int64, Float32} } } }, dem_band::AbstractArray{T}, noDataValue::Real ) where {T <: Number}
    rows, cols = size(dem_band)
    indexes = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]
    # For each cell of the dem's band
    for r in 1:rows, c in 1:cols
        if dem_band[r, c] != noDataValue
            # Indexes of adjacent cells
            for (i, j) in indexes
                if ( r+i >= 1 && r+i <= rows ) && ( c+j >= 1 && c+j <= cols ) && dem_band[r+i, c+j] != noDataValue
                    if !looseIn((i, j), mat[r, c])
                        push!(
                            mat[r, c],
                            ( i, j, dem_band[r+i, c+j] - dem_band[r, c] )
                        )
                    end
                    if !looseIn((-i, -j), mat[r+i, c+j])
                        push!(
                            mat[r+i, c+j],
                            ( -i, -j, dem_band[r+i, c+j] - dem_band[r, c] )
                        )
                    end
                end
            end
        end
    end
end

function connectivity( dem_band::Matrix{T}, batch_size::Integer, noDataValue::Real ) where {T <: Number}
    rows, cols = size(dem_band)
    n, m = ceil.( Int64, [rows, cols] ./ (batch_size - 1) )
    mat = Array{ Union{ Missing, Vector{Tuple{Int64, Int64, Float32} } } }(missing, rows, cols)
    for r in 1:rows, c in 1:cols
        if dem_band[r, c] != noDataValue
            mat[r, c] = Vector{Tuple{Int64, Int64, Float32}}()
        end
    end
    for i in 1:n, j in 1:m
        # Find the starting and ending indexes for the current slice of the matrix
         # if the batch is one of the ending ones its ending index will be the size of the matrix for that dimension
        rows_range::UnitRange{Int64} = ( (batch_size - 1) * (i - 1) + 1 ) : ( i != n ? (batch_size - 1) * i + 1 : rows )
        cols_range::UnitRange{Int64} = ( (batch_size - 1) * (j - 1) + 1 ) : ( j != m ? (batch_size - 1) * j + 1 : cols )
        # Use the indexes to run the function on a view of the matrix (passing also the corresponding view of the dem)
        connectivity_batch!(
 #= IL TIPO DEL VETTORE CHE CONTIENE mat E dem_band NON E' STABILE
            view.(
                [mat, dem_band],
                Ref( ( (batch_size - 1) * (i - 1) + 1 ) : ( i != n ? (batch_size - 1) * i + 1 : rows ) ),
                Ref( ( (batch_size - 1) * (j - 1) + 1 ) : ( j != m ? (batch_size - 1) * j + 1 : cols ) )
            )...,
 =#
            view(mat, rows_range, cols_range), view(dem_band, rows_range, cols_range), noDataValue )
    end
    return mat
end



function X_connectivity_batch!( mat::AbstractArray{ Union{ Missing, Vector{Tuple{Int64, Int64, Float32} } } }, dem_band::AbstractArray{T}, noDataValue::Real ) where {T <: Number}
    rows, cols = size(dem_band)
    indexes = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]
    # For each cell of the dem's band
    for r in 1:rows, c in 1:cols
        if dem_band[r, c] != noDataValue
            # Indexes of adjacent cells
            for (i, j) in indexes
                if ( r+i >= 1 && r+i <= rows ) && ( c+j >= 1 && c+j <= cols ) && dem_band[r+i, c+j] != noDataValue
                    if dem_band[r, c] > dem_band[r+i, c+j] && !looseIn( (i, j), mat[r, c] )
                        push!(
                            mat[r, c],
                            ( i, j, dem_band[r, c] - dem_band[r+i, c+j] )
                        )
                    end
                    if dem_band[r, c] < dem_band[r+i, c+j] && !looseIn( (-i, -j), mat[r+i, c+j] )
                        push!(
                            mat[r+i, c+j],
                            ( -i, -j, -(dem_band[r, c] - dem_band[r+i, c+j]) )
                        )
                    end
                end
            end
        end
    end
end

function X_connectivity( dem_band::Matrix{T}, batch_size::Integer, noDataValue::Real ) where {T <: Number}
    rows, cols = size(dem_band)
    n, m = ceil.( Int64, [rows, cols] ./ (batch_size - 1) )
    mat = Array{ Union{ Missing, Vector{Tuple{Int64, Int64, Float32} } } }(missing, rows, cols)
    for r in 1:rows, c in 1:cols
        if dem_band[r, c] != noDataValue
            mat[r, c] = Vector{Tuple{Int64, Int64, Float32}}()
        end
    end
    for i in 1:n, j in 1:m
        # Find the starting and ending indexes for the current slice of the matrix
         # if the batch is one of the ending ones its ending index will be the size of the matrix for that dimension
        rows_range::UnitRange{Int64} = ( (batch_size - 1) * (i - 1) + 1 ) : ( i != n ? (batch_size - 1) * i + 1 : rows )
        cols_range::UnitRange{Int64} = ( (batch_size - 1) * (j - 1) + 1 ) : ( j != m ? (batch_size - 1) * j + 1 : cols )
        # Use the indexes to run the function on a view of the matrix (passing also the corresponding view of the dem)
        X_connectivity_batch!( view(mat, rows_range, cols_range), view(dem_band, rows_range, cols_range), noDataValue )
    end
    return mat
end







import ArchGDAL as agd


dtm_file = split( @__DIR__ , "\\Porting\\")[1] * "\\Mappe\\DTM_wgs84.tiff"
dtm = agd.readraster(dtm_file)
band = agd.getband(dtm, 1)
band_mat = agd.read(band)
ndv = agd.getnodatavalue(band)
test1 = band[4001:5000, 6001:7000]
test1b = band[4001:6050, 6001:8050]
test2 = band[4501:4600, 6501:6600]
test3 = band[4001:4005, 6001:6005]
test4 = [ 10.0 10.0 10.0  3.0 10.0; 10.0 10.0  5.0 10.0 10.0; 10.0 10.0  8.0 10.0 10.0; 10.0 10.0  6.0  7.0 10.0; 10.0 10.0  4.0  6.0 10.0 ]
test5 = [ 1.0  ndv 10.0  3.0 10.0; 10.0  ndv  5.0 10.0 10.0; ndv 10.0  8.0 10.0 10.0;16.0 10.0  6.0  7.0 10.0;10.0  2.0  4.0  6.0 24.0 ]



mat = Array{ Union{ Missing, Vector{Tuple{Int64, Int64, Float32} } } }(missing, 1000, 1000)
for r in 1:5, c in 1:5
    if test1[r, c] != ndv
        mat[r, c] = Vector{Tuple{Int64, Int64, Float32}}()
    end
end



using BenchmarkTools

@code_warntype connectivity_batch!( mat, test5, ndv )
@code_warntype X_connectivity_batch!( mat, test5, ndv )

@code_warntype connectivity(test5, 5, ndv)
@code_warntype X_connectivity(test5, 5, ndv)


dim = 2048

@time mat = connectivity(test1b, dim, ndv)
@time mat = X_connectivity(test1b, dim, ndv)

@benchmark mat = connectivity(test1b, dim, ndv)
@benchmark mat = X_connectivity(test1b, dim, ndv)

mat = nothing
GC.gc()



@time mat = connectivity( band_mat, 1024, ndv )
@time mat = X_connectivity( band_mat, 1024, ndv )



#=
test1b
    256
        time
            1.415754 seconds (5.65 M allocations: 764.729 MiB, 29.85% gc time)

            0.831210 seconds (3.92 M allocations: 421.320 MiB, 42.32% gc time)

        benchmark
            BenchmarkTools.Trial: 3 samples with 1 evaluation.
            Range (min … max):  1.671 s …    2.135 s  ┊ GC (min … max): 42.34% … 42.65%
            Time  (median):     1.862 s               ┊ GC (median):    46.38%
            Time  (mean ± σ):   1.889 s ± 233.568 ms  ┊ GC (mean ± σ):  43.78% ±  2.25%
            Memory estimate: 764.73 MiB, allocs estimate: 5648918.

            BenchmarkTools.Trial: 7 samples with 1 evaluation.
            Range (min … max):  506.997 ms …    1.063 s  ┊ GC (min … max):  0.00% … 52.39%
            Time  (median):     730.994 ms               ┊ GC (median):    33.78%
            Time  (mean ± σ):   798.933 ms ± 192.962 ms  ┊ GC (mean ± σ):  36.31% ± 16.76%
            Memory estimate: 421.32 MiB, allocs estimate: 3917713.
    
    1024
        time
            1.511550 seconds (5.65 M allocations: 764.716 MiB, 34.69% gc time)

            0.891455 seconds (3.92 M allocations: 421.307 MiB, 43.98% gc time)
            
        benchmark
            BenchmarkTools.Trial: 3 samples with 1 evaluation.
            Range (min … max):  1.848 s …   2.010 s  ┊ GC (min … max): 42.80% … 42.39%
            Time  (median):     1.934 s              ┊ GC (median):    44.05%
            Time  (mean ± σ):   1.931 s ± 81.007 ms  ┊ GC (mean ± σ):  43.13% ±  0.95%
            Memory estimate: 764.72 MiB, allocs estimate: 5648846.

            BenchmarkTools.Trial: 7 samples with 1 evaluation.
            Range (min … max):  501.153 ms …    1.011 s  ┊ GC (min … max):  0.00% … 50.24%
            Time  (median):     748.714 ms               ┊ GC (median):    31.93%
            Time  (mean ± σ):   774.821 ms ± 168.350 ms  ┊ GC (mean ± σ):  35.65% ± 16.72%
            Memory estimate: 421.31 MiB, allocs estimate: 3917641.

    2048
        time
            1.565435 seconds (5.65 M allocations: 764.715 MiB, 34.72% gc time)

            0.907934 seconds (3.92 M allocations: 421.306 MiB, 42.38% gc time)

        benchmark
            BenchmarkTools.Trial: 3 samples with 1 evaluation.
            Range (min … max):  2.009 s …   2.126 s  ┊ GC (min … max): 46.20% … 39.91%
            Time  (median):     2.062 s              ┊ GC (median):    41.14%
            Time  (mean ± σ):   2.066 s ± 58.489 ms  ┊ GC (mean ± σ):  42.20% ±  3.43%
            Memory estimate: 764.72 MiB, allocs estimate: 5648841.

            BenchmarkTools.Trial: 5 samples with 1 evaluation.
            Range (min … max):  640.191 ms …    1.528 s  ┊ GC (min … max):  0.00% … 60.09%
            Time  (median):        1.044 s               ┊ GC (median):    41.98%
            Time  (mean ± σ):      1.060 s ± 335.365 ms  ┊ GC (mean ± σ):  44.70% ± 23.43%
            Memory estimate: 421.31 MiB, allocs estimate: 3917636.

full dtm
    256
        time
            2429.382277 seconds (59.62 M allocations: 7.986 GiB, 86.59% gc time, 0.02% compilation time)
=#



using DataFrames
using CSV





# Vecchio connectivity: 985.208187 seconds (118.07 M allocations: 13.669 GiB, 90.02% garbage collection time, 0.02% compilation time)
# Nuovo connectivity: 763.825265 seconds (88.54 M allocations: 11.471 GiB, 86.02% gc time, 0.01% compilation time)
@time mat = connectivity(band_mat, ndv)

# 70.208267 seconds (71.14 M allocations: 7.088 GiB, 60.86% gc time, 0.11% compilation time)
@time dmat = direct_connectivity(band_mat, ndv)


df = DataFrame(mat, :auto)
CSV.write("D:\\Connectivity Matrix\\direct_connectivity.csv", df)
CSV.write("C:\\Users\\DAVIDE-FAVARO\\Desktop\\connectivity.csv", df)









import ArchGDAL as agd

include("../Library/Functions.jl")


mat = connectivity(band_mat, 2048, ndv)

rows, cols = size(mat)
dtm2 = agd.read(dtm_file)
refsys = agd.getproj(dtm2)
geotransform = agd.getgeotransform(dtm2)
#=
#   minX, maxY = Functions.getSidesDistances(dtm2)
minX, maxY = Functions.toCoords(dtm2, 4001, 7000)
a, b = Functions.getCellDims(dtm2)
=#
gtiff_driver = agd.getdriver("GTiff")

target_ds = agd.create( "C:\\Users\\DAVIDE-FAVARO\\Desktop\\connectivity2.tiff", driver=gtiff_driver, width=rows, height=cols, nbands=8, dtype=Float32 )
#   agd.setgeotransform!( target_ds, [ minX, a, 0.0, maxY, 0.0, b ] )
agd.setgeotransform!( target_ds, geotransform )
agd.setproj!( target_ds, refsys )
valNoData = -9999.0
bands = [ agd.getband( target_ds, i ) for i in 1:8 ]
agd.setnodatavalue!.(bands, valNoData)
agd.fillraster!.(bands, valNoData)
band_mats = agd.read.(bands)

index_dict = Dict( 
    (-1, -1) => 1,
    (-1, 0)  => 2,
    (-1, 1)  => 3,
    (0, -1)  => 4,
    (0, 1)   => 5,
    (1, -1)  => 6,
    (1, 0)   => 7,
    (1, 1)   => 8
)

for r in 1:rows, c in 1:cols
    if !ismissing(mat[r, c]) && !isempty(mat[r, c])
        for (i, j, val) in mat[r, c]
            band_mats[ index_dict[(i, j)] ][r, c] = val
        end
    end
end

for i in 1:8
    agd.write!( target_ds, band_mats[i], i )
end

agd.destroy(target_ds)






function createRasterizedConnectivity( file::AbstractString, dtm_path::AbstractString )
    mat = connectivity(band_mat, 2048, ndv)
    dtm = agd.read(dtm_path)
    rows, cols = size(mat)
    res = agd.create( file, driver=agd.getdriver("GTiff"), width=rows, height=cols, nbands=8, dtype=Float32 )
    agd.setgeotransform!( res, agd.getgeotransform(dtm) )
    agd.setproj!( res, agd.getproj(dtm) )
    valNoData = agd.getnodatavalue(dtm)
    bands = [ agd.getband( res, i ) for i in 1:8 ]
    agd.setnodatavalue!.(bands, valNoData)
    agd.fillraster!.(bands, valNoData)
    band_mats = agd.read.(bands)

    index_dict = Dict( 
        (-1, -1) => 1,
        (-1, 0)  => 2,
        (-1, 1)  => 3,
        (0, -1)  => 4,
        (0, 1)   => 5,
        (1, -1)  => 6,
        (1, 0)   => 7,
        (1, 1)   => 8
    )

    for r in 1:rows, c in 1:cols
        if !ismissing(mat[r, c]) && !isempty(mat[r, c])
            for (i, j, val) in mat[r, c]
                band_mats[ index_dict[(i, j)] ][r, c] = val
            end
        end
    end

    for i in 1:8
        agd.write!( target_ds, band_mats[i], i )
    end

    agd.destroy(target_ds)
end




#= TENTATIVO DI SISTEMARE IL TIFF
    import ArchGDAL as agd

    dtm_file = split( @__DIR__ , "\\Porting\\")[1] * "\\Mappe\\DTM_32.tiff"
    dtm = agd.read(dtm_file)
    res = agd.read("C:\\Users\\DAVIDE-FAVARO\\Desktop\\connectivity.tiff")

    dtm_band = agd.getband(dtm, 1)
    dtm_band_mat = agd.read(dtm_band)

    res_bands = [ agd.getband(res, i) for i in 1:8 ]
    res_band_mats = agd.read.(res_bands) 

    rows, cols = size(dtm_band)
    refsys = agd.getproj(dtm)
    geotransform = agd.getgeotransform(dtm)
    gtiff_driver = agd.getdriver("GTiff")

    target_ds = agd.create( "C:\\Users\\DAVIDE-FAVARO\\Desktop\\connectivity2.tiff", driver=gtiff_driver, width=rows, height=cols, nbands=8, dtype=Float32 ) 
    agd.setgeotransform!( target_ds, geotransform )
    agd.setproj!( target_ds, refsys )
    valNoData = -9999.0
    bands = [ agd.getband( target_ds, i ) for i in 1:8 ]
    agd.setnodatavalue!.(bands, valNoData)
    agd.fillraster!.(bands, valNoData)
    band_mats = agd.read.(bands)

    #=
    for i in 1:8
        for r in 1:rows, c in 1:cols
            if res_band_mats != 0.0
                band_mats[i][r, c] = res_band_mats[i][r, c]
            end
        end
    end
    =#
    for i in 1:8
        agd.write!( target_ds, res_band_mats[i], i )
    end

    agd.destroy(target_ds)
=#























function run_runoff( dem, source, target, landcover, soil_text::AbstractString, resolution::Integer, folder::AbstractString=".\\" )

 """ NON SO QUALE SIA L'EQUIVALENTE
    if not self.dem.isValid():
        QMessageBox.warning(self,"Warning", "The dem file is not valid" )
        return
 """

    if agd.geomdim(source) != 0
        throw(DomainError(source, "`source` must be a point"))
    end

    target_layer, landcover_layer, soil_layer = agd.getlayer([target, landcover, soil], 0) 

    if agd.geomdim(target_layer) != 2
        throw(DomainError(source, "`target` must be a polygon"))
    end

    if agd.geomdim(landcover_layer) != 2
        throw(DomainError(source, "Not a valid `landcover` geometry"))
    end

    refsys = agd.getspatialref(source)

    if agd.getspatialref(target_layer) != agd.getspatialref(source) || agd.getspatialref(landcover) != agd.getspatialref(source) || agd.getspatialref(dem) != agd.getspatialref(source)
        throw(DomainError("The reference systems are not uniform. Aborting analysis."))
    end

    path_temp_landcover = folder * "\\temp_lc.tiff"
    path_temp_soil = folder * "\\temp_soil.tiff"


 """ NON SO COSA SIANO QUESTE VARIABILI
    self.text_p = str(self.combo_fieldp.currentText())
 """

    clc_list = Functions.cn_list_extract()

    soil_control = 0

 """ PRINT DI COSE
    messaggio+='ALGORITMO UTILIZZATO: calcolo della separazione delle componenti infiltrazione e ruscellamento tramite metodo SCS-CN; US Department of Agriculture Soil Conservation Service, 1972. National Engineering Handbook, Section 4, Hydrology. US Government Printing Office, Washington, DC, 544pp.\n\n'
 """

    output_path = folder * "\\runoff.tiff"
    intervallo = max(getCellDims(dem))




    # NON SO SE FUNZIONA
    bbox_src = agd.boundingbox(source)
    bbox_trgt = agd.boundingbox(target)
    area = agd.union( bbox_src, bbox_trgt )
    # Geometry containing both the source and the target
     # Se union non ritorna una geometria unica ma un'unica geometria composta di due elementi disgiunti
     #   area = agd.boundingbox(agd.union( bbox_src, bbox_trgt ))
    area = agd.union( bbox_src, bbox_trgt )

 """ PRENDE LA PORZIONE DI `landcover` RAPPRESENTATA DA `area`
    lc_clip_proc = processing.run('qgis:clip', {'INPUT':self.lc, 'OVERLAY':self.areastudio, 'OUTPUT':self.path_working+'/clip.gpkg'})
    lc_clip=QgsVectorLayer(lc_clip_proc['OUTPUT'], 'lc_clip', 'ogr')
    lc_clip.setCrs(self.source.crs())

    #path__layer_lc=lc_clip['OUTPUT'].dataProvider().dataSourceUri()
    path__layer_lc=lc_clip.dataProvider().dataSourceUri()
    path_lc=path__layer_lc.split("|")
    source_ds_lc = ogr.Open(path_lc[0])
    lc_layer = source_ds_lc.GetLayer()
 """
    # NON SONO CERTO PRESERVI LE INFORMAZIONI DI landcover
    landcover_clip = agd.intersects(landcover, area)
    landcover_layer = agd.getlayer(landcover_clip, 0)   

    rows = agd.height(landcover_layer)
    cols = agd.width(landcover_layer)

    # NON SO COME OTTENERE minX E maxY
    minX, maxY = ?

    landcover_ds = agd.create( path_temp_landcover, gtiff_driver, rows, cols, 1, agd.GDAL.GDT_Float32 )
    agd.setgeotransform!(landcover_ds, [ minX, 25.0, 0.0, maxY, 0.0, -25.0 ])
    agd.setproj!(landcover_ds, refsys)
    band_lc = agd.getband(landcover_ds, 1)
    agd.setnodatavalue!( band_lc, Float64(valNoData) )
    bandlc = agd.read(band_lc)
    agd.fillraster!(bandlc, valNoData)

    
 """ TRASFORMA IN RASTER LA PORZIONE DI `landcover` PRESA PRIMA """
    gdal.RasterizeLayer(lc_ds, [1], lc_layer,options=["ATTRIBUTE="+self.text_lcfield])
    lc_ds=None

    lc_layer=QgsRasterLayer(self.path_temp_lc,"lc_layer")
 """"""
    landcover_layer = agd.gdalrasterize( x -> x, landcover_ds )



    if soil_text == "Valore campo"
        soil_control = 1

     # NON SO SE SIA EQUIVALENTE
        #   source_ds_soil = ogr.Open(path_lc[0])
        #   soil_layer = source_ds_soil.GetLayer()
        agd.getlayer(path_lc, 0)

        soil_ds = agd.create( path_temp_soil, gtiff_driver, round(Int64, x_res), round(Int64, y_res), 1, agd.GDAL.GDT_Float32 )
        agd.setgeotransform!(soil_ds, [ x_min, 25.0, 0.0, y_max, 0.0, -25.0 ])
        agd.setproj!(soil_ds, refsys)
        band_sl = agd.getband(soil_ds, 1)
        agd.setnodatavalue!( band_sl, Float64(valNoData) )
        bandsl = agd.read(band_sl)
        agd.fillraster!(bandsl, valNoData)

     """ RASTERIZZA IL  VETTORIALE DEL TIPO DI SUOLO """
        gdal.RasterizeLayer(lc_ds, [1], soil_layer,options=["ATTRIBUTE="+self.text_soilfield])
        soil_ds=None

        soil_layer=QgsRasterLayer(self.path_temp_soil,"soil_layer")
     """"""
        soil_layer = agd.gdalrasterize( x -> x, landcover_ds )
    end

    
 """ CALCOLA UN VETTORIALE TRAMITE `r.drain` SI PUO' FARE CON Omniscape.jl O SIMILI"""
  # FORSE C'E' ANCHE Anasol E IL FRAMEWORK MADS

    grass_area=str(x_min)+','+str(x_max)+','+str(y_min)+','+str(y_max)+' ['+str(self.areastudio.crs().authid())+']'
    # grass_coord=str(x_source)+','+str(y_source)+' ['+str(self.source.crs().authid())+']'

    namewatershed=self.path_working+'/watershed'
    namedrain=self.path_working+'/wshed.shp'

    params = { 'GRASS_RASTER_FORMAT_OPT' : '','GRASS_REGION_CELLSIZE_PARAMETER' : 0, 'GRASS_REGION_PARAMETER' :grass_area,
               'GRASS_VECTOR_EXPORT_NOCAT' : False, '-a' : False, 'start_coordinates' : None,  '-n' : False,
               'input' : self.dem.dataProvider().dataSourceUri(),'-c' : True, 'drain' : namedrain, 'GRASS_MIN_AREA_PARAMETER' : 0.0001,
               'start_points' : self.source.dataProvider().dataSourceUri(), 'output' : namewatershed }

    waterwshed_proc = processing.run('grass7:r.drain', params)


    #aggiungo per controllo la viewshed alla toc
    #iface.addVectorLayer(namedrain,'watershed','ogr')
    #watershed=QgsProject.instance().mapLayersByName('watershed.shp')


    vdrain = QgsVectorLayer(namedrain, 'vdrain', 'ogr')

    idxlevel = self.source.fields().indexFromName(self.text_p)
    idxcat = vdrain.fields().indexFromName('cat')

    idxtargetname = self.target.fields().indexFromName(self.text_targetfield)
 """"""




    features = agd.getgeom.(collect(agd.getlayer(vdrain)))
    nfeat = 0
    polygons_t = collect(target_layer) 

    #   start_time = time.time()

    for f in features
        length = agd.geomlength(f)
        currentdistance = intervallo
        nfeat += 1
        featlines = []



        #   firstpoint=geom.interpolate(0)
        firstpoint = agd.getpoint(geom, 0)
        old_x = agd.getx(firstpoint, 0)
        old_y = agd.gety(firstpoint, 0)

        fileoutput = folder * "\\drain$nfeat.shp"



     """ CREA UN LAYER
        vline = QgsVectorLayer("LineString?crs=EPSG:"+self.refsys, "drain"+str(nfeat), "memory")

        prline = vline.dataProvider()
        prlfield=prline.addAttributes( [ QgsField("concentrazione", QVariant.Double) ] )

        idf=f.attributes()[idxcat]
        feat_drain = next(self.source.getFeatures(QgsFeatureRequest().setFilterFid(idf-1)))
        p00=feat_drain.attributes()[idxlevel]
     """
        vline = agd.createlayer( x -> x,  name="drain$nfeat", geom=agd.wkbLineString, spatialref=refsys )
        agd.addfielddefn!(vline, "concentrazione", agd.OFTReal)
     # idxcat SI OTTIENE DA vdrain LA CUI GENERAZIONE NON E' ANCORA IMPLEMENTATA
        idf = agd.getfield(f, idxcat)
     # NON SONO SICURO DEL idf-1 E NEMMENO DELL'EFFETTO DI next
        feat_drain = getFeatureByFid( collect(agd.getlayer(source)), idf )
      # idxlevel SI OTTIENE DA vdrain LA CUI GENERAZIONE NON E' ANCORA IMPLEMENTATA  
        p00 = agd.getfield(feat_drain, idxlevel)
        
        index_progress = 0
        while currentdistance < length
            if index_progress == 0
                p0 = p00
            else
                p0 = pe
            end

            #   point = geom.interpolate(currentdistance)
            point = agd.getpoint(geom, currentdistance)
            x = agd.getx(point, 0)
            y = agd.gety(point, 0)
            r, c = toIndexes(landcover_layer, x, y)
            clc = landcover_layer[r, c]
            if soil_control == 1
                r, c = toIndexes(soil_layer, x, y)
                soil = soil_layer[r, c]
            else
                soil = text_soil
            end
            try
             # MANCA "classisoil"
                cn = listaclc[clc][ classisoil[soil] ]
                S = 254.0((100 / cn) - 1)
            catch
                S = 0
            end

            pcheck = (p0 - 0.2S)^2 / (p0 - 0.2S + S)
            if pcheck > 0.2S
                pe = pcheck


             """ CREA UNA FEATURE
               fetline = QgsFeature()
               fetline.setGeometry( QgsGeometry.fromPolyline( [QgsPoint(old_x,old_y),QgsPoint(x,y)] ))
               fetline.initAttributes(1)
               fetline.setAttribute(0,pe)
               vline.updateFeature(fetline)

               featlines.append(fetline)
             """
                fetline = agd.createfeature( x -> x, agd.getfeaturedefn(vline) )
                line = agd.createlinestring( Flloat64.([old_x, old_y]), Flloat64.([x, y]) )
                agd.setgeom!(fetline, line)
                agd.fillunsetwithdefault!(fetline)
                agd.setfield!(fetline, 0, pe)
                push!(featlines, fetline)

                index_progress += 1

                for polygon in polygons_t
                    p_geom = agd.getgeom(polygon)
                    if agd.within(point, p_geom)
                        #   nometarget = pol_t.attributes()[idxtargetname]
                        nometarget = agd.getfield(polygon, idxtargetname)
                        println("\nIl vettore drain$nfeat ha raggiunto l'area bersaglio denominata $nometarget con un volume pari a: $(round(pe,3))mm\n")
                        currentdistance = length + 1
                    end
                end
            else
                pe = 0
                currentdistance = length + 1
            end
                old_x = x
                old_y = y
                currentdistance += intervallo
            end

         """ AGGIUNGE LE FEATURES CREATE
            prline.addFeatures(featlines)
            vline.updateFields()
            QgsProject.instance().addMapLayer(vline)
         """
            agd.addfeature!.(Ref(vline), featlines)
    end
end



end # module








































































function run_runoff( dem, source, target, landcover, soil_text::AbstractString, resolution::Integer, folder::AbstractString=".\\" )

    """ NON SO QUALE SIA L'EQUIVALENTE
       if not self.dem.isValid():
           QMessageBox.warning(self,"Warning", "The dem file is not valid" )
           return
    """
   
       if agd.geomdim(source) != 0
           throw(DomainError(source, "`source` must be a point"))
       end
   
       if agd.geomdim(targer) != 2
           throw(DomainError(source, "`target` must be a polygon"))
       end
   
       if agd.geomdim(landcover) != 2
           throw(DomainError(source, "Not a valid `landcover` geometry"))
       end
   
       if agd.getspatialref(target) != agd.getspatialref(source) || agd.getspatialref(landcover) != agd.getspatialref(source) || agd.getspatialref(dem) != agd.getspatialref(source)
           throw(DomainError("The reference systems are not uniform. Aborting analysis."))
       end
   
       refsys = agd.importEPSG(agd.fromWKT(agd.getspatialref(source)))
   
       path_temp_landcover = folder * "\\temp_lc.tiff"
       path_temp_soil = folder * "\\temp_soil.tiff"
   
   
    """ NON SO COSA SIANO QUESTE VARIABILI
       self.text_p = str(self.combo_fieldp.currentText())
    """
   
       clc_list = Functions.cn_list_extract()
   
       soil_control = 0
   
    """ PRINT DI COSE
       messaggio+='ALGORITMO UTILIZZATO: calcolo della separazione delle componenti infiltrazione e ruscellamento tramite metodo SCS-CN; US Department of Agriculture Soil Conservation Service, 1972. National Engineering Handbook, Section 4, Hydrology. US Government Printing Office, Washington, DC, 544pp.\n\n'
    """
   
       output_path = folder * "\\runoff.tiff"
   
   
   
   
   
   
   
   
   
   
   
   
   
   
       valNoData = -9999
   
       gtiff_driver = agd.getdriver("GTiff")
       target_ds = agd.create( output_path, gtiff_driver, round(Int64, ), round(Int64, ), 1, agd.GDAL.GDT_Float32 )
       agd.setgeotransform!(target_ds, [ , resolution, 0.0, , 0.0, -resolution ])
       agd.setproj!(target_ds, refsys)
    """ NON SO QUALE SIA IL COMANDO PER SETTARE I METADATI CON `ArchGDAL`
       target_ds.SetMetadata(
           Dict(
               "credits" => "Envifate - Francesco Geri, Oscar Cainelli, Paolo Zatelli, Gianluca Salogni, Marco Ciolli - DICAM Università degli Studi di Trento - Regione Veneto",
               "modulo" => "Analisi ruscellamento",
               "descrizione" => "Analisi di ruscellamento di un inquinante attraverso il metodo della separazione delle componenti",
               "srs" => refsys,
               "data" => today()
           )
       )
    """
       band1 = agd.getband(target_ds, 1)
       agd.setnodatavalue!( band1, Float64(valNoData) )
       band = agd.read(band1)
       agd.fillraster!(band, valNoData)
       xsize = agd.width(band)
       ysize = agd.height(band)
   
   
   
   
   
       intervallo = getCellDims(dem)
   
   
   
   
   
   
   
    """ PRENDE LA PORZIONE DI `landcover` RAPPRESENTATA DA `area` """
       lc_clip_proc = processing.run('qgis:clip', {'INPUT':self.lc, 'OVERLAY':self.areastudio, 'OUTPUT':self.path_working+'/clip.gpkg'})
       lc_clip=QgsVectorLayer(lc_clip_proc['OUTPUT'], 'lc_clip', 'ogr')
       lc_clip.setCrs(self.source.crs())
   
       #path__layer_lc=lc_clip['OUTPUT'].dataProvider().dataSourceUri()
       path__layer_lc=lc_clip.dataProvider().dataSourceUri()
       path_lc=path__layer_lc.split("|")
       source_ds_lc = ogr.Open(path_lc[0])
       lc_layer = source_ds_lc.GetLayer()
    """"""
       
   
       landcover_ds = agd.create( path_temp_landcover, gtiff_driver, round(Int64, x_res), round(Int64, y_res), 1, agd.GDAL.GDT_Float32 )
       agd.setgeotransform!(landcover_ds, [ x_min, 25.0, 0.0, y_max, 0.0, -25.0 ])
       agd.setproj!(landcover_ds, refsys)
       band_lc = agd.getband(landcover_ds, 1)
       agd.setnodatavalue!( band_lc, Float64(valNoData) )
       bandlc = agd.read(band_lc)
       agd.fillraster!(bandlc, valNoData)
       xsize = agd.width(band)
       ysize = agd.height(band)
   
    """ TRASFORMA IN RASTER LA PORZIONE DI `landcover` PRESA PRIMA """
       gdal.RasterizeLayer(lc_ds, [1], lc_layer,options=["ATTRIBUTE="+self.text_lcfield])
       lc_ds=None
   
       lc_layer=QgsRasterLayer(self.path_temp_lc,"lc_layer")
    """"""
   
   
   
   
       if soil_text == "Valore campo"
           soil_control = 1
   
        # NON SO SE SIA EQUIVALENTE
           #   source_ds_soil = ogr.Open(path_lc[0])
           #   soil_layer = source_ds_soil.GetLayer()
           agd.getlayer(path_lc, 0)
   
           soil_ds = agd.create( path_temp_soil, gtiff_driver, round(Int64, x_res), round(Int64, y_res), 1, agd.GDAL.GDT_Float32 )
           agd.setgeotransform!(soil_ds, [ x_min, 25.0, 0.0, y_max, 0.0, -25.0 ])
           agd.setproj!(soil_ds, refsys)
           band_sl = agd.getband(soil_ds, 1)
           agd.setnodatavalue!( band_sl, Float64(valNoData) )
           bandsl = agd.read(band_sl)
           agd.fillraster!(bandsl, valNoData)
   
        """ RASTERIZZA IL  VETTORIALE DEL TIPO DI SUOLO """
           gdal.RasterizeLayer(lc_ds, [1], soil_layer,options=["ATTRIBUTE="+self.text_soilfield])
           soil_ds=None
   
           soil_layer=QgsRasterLayer(self.path_temp_soil,"soil_layer")
        """"""
       end
   
       
    """ CALCOLA UN VETTORILE TRAMITE `r.drain` """
       grass_area=str(x_min)+','+str(x_max)+','+str(y_min)+','+str(y_max)+' ['+str(self.areastudio.crs().authid())+']'
       # grass_coord=str(x_source)+','+str(y_source)+' ['+str(self.source.crs().authid())+']'
   
       namewatershed=self.path_working+'/watershed'
       namedrain=self.path_working+'/wshed.shp'
   
       params = { 'GRASS_RASTER_FORMAT_OPT' : '','GRASS_REGION_CELLSIZE_PARAMETER' : 0, 'GRASS_REGION_PARAMETER' :grass_area,
                  'GRASS_VECTOR_EXPORT_NOCAT' : False, '-a' : False, 'start_coordinates' : None,  '-n' : False,
                  'input' : self.dem.dataProvider().dataSourceUri(),'-c' : True, 'drain' : namedrain, 'GRASS_MIN_AREA_PARAMETER' : 0.0001,
                  'start_points' : self.source.dataProvider().dataSourceUri(), 'output' : namewatershed }
   
       waterwshed_proc = processing.run('grass7:r.drain', params)
   
   
       #aggiungo per controllo la viewshed alla toc
       #iface.addVectorLayer(namedrain,'watershed','ogr')
       #watershed=QgsProject.instance().mapLayersByName('watershed.shp')
   
   
       vdrain = QgsVectorLayer(namedrain, 'vdrain', 'ogr')
   
       idxlevel = self.source.fields().indexFromName(self.text_p)
       idxcat = vdrain.fields().indexFromName('cat')
   
       idxtargetname = self.target.fields().indexFromName(self.text_targetfield)
    """"""
   
   
   
   
       features = agd.getgeom.(collect(agd.features(vdrain)))
       nfeat = 0
       polygons_t = collect(agd.getfeature(target)) 
   
       #   start_time = time.time()
   
       for f in features
           length = agd.geomlength(f)
           currentdistance = intervallo
           nfeat += 1
           featlines = []
   
   
   
           #   firstpoint=geom.interpolate(0)
           firstpoint = agd.getpoint(geom, 0)
           old_x = agd.getx(firstpoint, 0)
           old_y = agd.gety(firstpoint, 0)
   
           fileoutput = folder * "\\drain$nfeat.shp"
   
   
   
        """ NON SO COSA FACCIA STA ROBA """
           vline = QgsVectorLayer("LineString?crs=EPSG:"+self.refsys, "drain"+str(nfeat), "memory")
   
           prline = vline.dataProvider()
           prlfield=prline.addAttributes( [ QgsField("concentrazione", QVariant.Double) ] )
   
           idf=f.attributes()[idxcat]
           feat_drain = next(self.source.getFeatures(QgsFeatureRequest().setFilterFid(idf-1)))
           p00=feat_drain.attributes()[idxlevel]
        """"""
   
           
           index_progress = 0
           while currentdistance < length
               if index_progress == 0
                   p0 = p00
               else
                   p0 = pe
               end
   
               #   point = geom.interpolate(currentdistance)
               point = agd.getpoint(geom, currentdistance)
               x = agd.getx(point, 0)
               y = agd.gety(point, 1)
               clc = extract_values(lc_layer, x, y)
               if soil_control == 1
                   r, c = toIndexes(soil_layer, x, y)
                   soil = soil_layer[r, c]
               else
                   soil = text_soil
               end
               try
                # MANCA "classisoil"
                   cn = listaclc[clc][classisoil[soil]]
                   S = 254.0((100 / cn) - 1)
               catch
                   S = 0
               end
   
               pcheck = (p0 - 0.2S)^2 / (p0 - 0.2S + S)
               if pcheck > 0.2S
                   pe = pcheck
   
   
                """ NON SO COSA FACCIA STA ROBA """
                  fetline = QgsFeature()
                  fetline.setGeometry( QgsGeometry.fromPolyline( [QgsPoint(old_x,old_y),QgsPoint(x,y)] ))
                  fetline.initAttributes(1)
                  fetline.setAttribute(0,pe)
                  vline.updateFeature(fetline)
   
                  featlines.append(fetline)
                """"""
                
   
                   index_progress += 1
   
                   for polygon in polygons_t
                       p_geom = agd.getgeom(polygon)
                       if agd.within(point, p_geom)
                           nometarget = pol_t.attributes()[idxtargetname]
                        """ MESSAGIO
                           messaggio = "\nIl vettore drain$nfeat ha raggiunto l'area bersaglio denominata $nometarget con un volume pari a: $(round(pe,3))mm\n"
                           self.console.appendPlainText(messaggio)
                        """
                           currentdistance = length + 1
                       end
                   end
               else
                   pe = 0
                   currentdistance = length + 1
               end
                   old_x = x
                   old_y = y
                   currentdistance += intervallo
               end
   
   
            """ NON SO COSA FACCIA STA ROBA """
               prline.addFeatures(featlines)
               vline.updateFields()
               QgsProject.instance().addMapLayer(vline)
            """"""
   
       end
   end
   