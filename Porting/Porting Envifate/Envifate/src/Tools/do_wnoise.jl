module WaterNoise

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
            os.system("start "+os.path.dirname(__file__)+"/../tutorial/manuale_envifate_rumore_in_acqua.pdf")
        if platform.uname()[0]=="Linux":
            os.system("xdg-open "+os.path.dirname(__file__)+"/../tutorial/manuale_envifate_rumore_in_acqua.pdf")
        else:
            os.system("open "+os.path.dirname(__file__)+"/../tutorial/manuale_envifate_rumore_in_acqua.pdf")
=#


import ArchGDAL as agd

include("..\\Library\\Functions.jl")


function run_f1(S,T)
    return √(0.78(S / 35)) * ℯ^(T/26)
end



function run_f2(T)
    return 42ℯ^(T/17)
end



function run_spread( dem, source, depth::Real, salinity::Real, pH::Real, temperature::Real, frequencys::Real, resolution::Integer, folder::AbstractString=".\\" )

    if pH < 0 || pH > 14
        throw(DomainError(pH, "`pH` must be a number between 0 and 14"))
    end

    if agd.geomdim(source) != 0
        throw(DomainError(source, "`source` must be a point"))
    end

    if agd.geomdim(area) != 2
        throw(DomainError(source, "`area` must be a polygon"))
    end

    refsys = agd.getspatialref(source)
    measure_dist = 15.0

    path_temp_lc = folder*"\\temp_lc.tiff"




    for freq in frequencys
        output_path = folder*"\\sound_level_$freq.tiff"

        # messaggio+='ALGORITMO UTILIZZATO: Ainslie, M. A., & McColm, J. G. (1998). A simplified formula for viscous and chemical absorption in sea water. The Journal of the Acoustical Society of America, 103(3), 1671-1672.)\n\n'

        nfeature=0
        features = collect(agd.getlayer(source))
        for feature in features
            geom = agd.getgeom(feature)
            x_source = agd.getx(geom,0)
            y_source = agd.gety(geom,0)
 
            soundlevel = agd.getfield(feature, :level)
            nfeature+=1

            #   start_time = time.time()
            
            for row in 1:rows
                for col in 1:cols
                    x, y = (col, row) * resolution + (x_min, y_min) + (resolution / 2)

                    # calcolo distanza
                    Δx = x - x_source
                    Δy = y - y_source
                    dist = √( Δy^2 + Δx^2 )

                    f1 = √(0.78(salinity / 35))ℯ^(temperature / 26)
                    f2 = 42ℯ^(temperature / 17)
                    α1 = 0.106( (f1 * freq_f^2) / (freq_f^2 + f1^2) )ℯ^( (pH - 8) / 0.56 )
                    α2 = 0.52( 1 + (temperature / 43) ) * (salinity / 35) * ( (f2 * freq_f^2) / (freq_f^2 + f1^2) )ℯ^(-depth / 6)
                    α3 = 0.00049freq_f^2 * ℯ^( -((temperature / 27) + (depth / 17)) )
                    α = α1 + α2 + α3
                    tl = 20log(dist) + α
                    total_loss = soundlevel - tl

                    if nfeature == 1
                        if total_loss > 0
                            outData[row,col]=total_loss
                        else
                            outData[row,col]=0
                        end
                    else
                        if total_loss > 0
                            outData[row,col]=outData[row,col]+total_loss
                        else
                            outData[row,col]=outData[row,col]
                        end
                    end

                    ##### inizio scrittura file temporanei
                    # outData_eucdist[row,col]=self.soundlevel
                    # outData_aal[row,col]=ssl_loss
                    # outData_mvl[row,col]=sslaal
                    # outData_bar[row,col]=sslaal1
                    # outData_wind[row,col]=wind_loss
                    ##### fine scrittura file temporanei
                end
            end
        end







        gtiff_driver = agd.getdriver("GTiff")
        target_ds = agd.create( path, gtiff_driver, rows, cols, 1, agd.GDAL.GDT_Float32 )
     # NON SONO CERTO CHE IL GEOTRASFORM VADA BENE
        agd.setgeotransform!( target_ds, [ minX, resolution, 0.0, maxY, 0.0, -resolution ] )
        agd.setproj!( target_ds, refsys )
     """ NON SO QUALE SIA IL COMANDO PER SETTARE I METADATI CON `ArchGDAL`
        target_ds.SetMetadata(
            Dict(
                "credits" => "Envifate - Francesco Geri, Oscar Cainelli, Paolo Zatelli, Gianluca Salogni, Marco Ciolli - DICAM Università degli Studi di Trento - Regione Veneto",
                "modulo" => "Dispersione rumore in acqua",
                "descrizione" => "Analisi della dispersione acustica in acqua",
                "srs" => refsys,
                "data" => today()
            )
        )
     """
        valNoData = -9999.0
        band1 = agd.getband( target_ds, 1 )
        agd.setnodatavalue!(band1, valNoData)
        agd.fillraster!(band1, valNoData)
        band = agd.read(band1)





 """
        outData_raster=outData[::-1]
        band.WriteArray(outData_raster)

        band= None
        target_ds = None

        base_raster_name=os.path.basename(self.path_output)
        raster_name=os.path.splitext(base_raster_name)[0]
        self.outputlayer=self.iface.addRasterLayer(self.path_output, raster_name)

        layer=None
        for lyr in list(QgsProject.instance().mapLayers().values())
            if lyr.name() == raster_name
                layer = lyr
            end
        end
"""
end

end # module
