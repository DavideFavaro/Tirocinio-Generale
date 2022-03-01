module Runoffs



import ArchGDAL as agd

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



# Given a tuple with both values in "-1:1" (except "(0, 0)") return a value in "1:8"
function hash_adjacent( t::Tuple{Int64, Int64} )
    # 3x3 matrix linearization
     # Its necessary to sum 2 to both values of t as they reppresent index displacements from the center of the matrix
      # the linearization process is based on starting index 0, hence the need to sum 2 instead of 3
      # for the same reason 1 needs to be substracted from the first index, ending up with "t[1] + 1"  
    res = 3(t[1] + 1) + (t[2] + 2)
    # The cell of index 5 is the center and it's not needed
    return res >= 5 ? res - 1 : res 
end

# Given two values in "-1:1" (except "0" and "0") return a value in "1:8"
function hash_adjacent( r::Int64, c::Int64 )
    # 3x3 matrix linearization
     # Its necessary to sum 2 to both values of t as they reppresent index displacements from the center of the matrix
      # the linearization process is based on starting index 0, hence the need to sum 2 instead of 3
      # for the same reason 1 needs to be substracted from the first index, ending up with "t[1] + 1"  
    res = 3(r + 1) + (c + 2)
    # The cell of index 5 is the center and it's not needed
    return res >= 5 ? res - 1 : res 
end

# Given a value in "1:8" return the corresponding tuple with both values in "-1:1", except "(0, 0)"
function hash_adjacent( i::Int64 )
    # To account for the fact that "i" = 5 for the tuple "(0, 0)" 
    if i > 4
        i += 1
    end
    r = ceil(Int64, i / 3)
    c = ((i - 1) % 3) + 1
    # Substract 2 to obtain the index displacement from the center
    return r - 2, c - 2 
end



# VEDERE @view, @inbound, @turbo, @fast, LoopedVectorization.jl e StableArrays.jl PER ULTERIORI OTTIMIZZAZIONI
function connectivity_batch!( mat, dem_band::AbstractArray{T}, noDataValue::Real ) where {T <: Number}
    rows, cols = size(dem_band)
    # For each cell of the dem's band
    @inbounds for r in 1:rows, c in 1:cols
        if dem_band[r, c] != noDataValue
            # Indexes of adjacent cells
            for (i, j) in keys(indexes)
                if ( r+i >= 1 && r+i <= rows ) && ( c+j >= 1 && c+j <= cols ) && dem_band[r+i, c+j] != noDataValue
                    if mat[ r, c, hash_adjacent(i, j) ] == noDataValue
                        mat[ r, c, hash_adjacent(i, j) ] = dem_band[r, c] - dem_band[r+i, c+j]
                    end
                    if mat[ r+i, c+j, hash_adjacent(-i, -j) ] == noDataValue
                        mat[ r+i, c+j, hash_adjacent(-i, -j) ] = dem_band[r+i, c+j] - dem_band[r, c]
                    end
                end
            end
        end
    end
end


function connectivity( dem_band::Matrix{T}, batch_size::Integer, noDataValue::Real ) where {T <: Number}
    rows, cols = size(dem_band)
    mat = fill( convert(Float32, noDataValue) , rows, cols, 8 )
    if batch_size == 1
        connectivity_batch!(mat, dem_band, noDataValue)
        return mat
    end
    n, m = ceil.( Int64, [rows, cols] ./ (batch_size - 1) )
    @inbounds for i in 1:n, j in 1:m
        # Find the starting and ending indexes for the current slice of the matrix
         # if the batch is one of the ending ones its ending index will be the size of the matrix for that dimension
        rows_range::UnitRange{Int64} = ( (batch_size - 1) * (i - 1) + 1 ) : ( i != n ? (batch_size - 1) * i + 1 : rows )
        cols_range::UnitRange{Int64} = ( (batch_size - 1) * (j - 1) + 1 ) : ( j != m ? (batch_size - 1) * j + 1 : cols )
        # Use the indexes to run the function on a view of the matrix (passing also the corresponding view of the dem)
        connectivity_batch!( view(mat, rows_range, cols_range, :), view(dem_band, rows_range, cols_range), noDataValue )
    end
    return mat
end


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

    for r in 1:rows, c in 1:cols
        if !ismissing(mat[r, c]) && !isempty(mat[r, c])
            for (i, j, val) in mat[r, c]
                band_mats[hash_adjacent(i, j)][r, c] = val
            end
        end
    end

    for i in 1:8
        agd.write!( target_ds, band_mats[i], i )
    end

    agd.destroy(target_ds)
end



























using Plots

import ArchGDAL as agd

include("../Library/Functions.jl")



# Find the direction with the greatest difference in height from (r, c)
function mindir( mat, r::Int64, c::Int64 )
    dir = 0
    min = Inf
    @inbounds for i in 1:8
        if mat[r, c, i] != -9999.0 && mat[r, c, i] < min
            dir = i
            min = mat[i][r, c]
        end
    end
    return (dir, min)
end
 #= MINDIR TEST VALUES 
    arr = [
        [ 1  2  -9999.0; 1  2  3 ],
        [ 3  4  -9999.0; 4  5  6 ],
        [ 5  6  -9999.0; 7  8  9 ],
        [ 7  8  -9999.0; 10 11 12 ],
        [ 9  10 -9999.0; 13 14 15 ],
        [ 11 12 -9999.0; 16 17 18 ],
        [ 13 14 -9999.0; 19 20 21 ],
        [ 17 16 -9999.0; 22 23 24 ]
    ]
 =#

# Return the amount of water flowing along the path formed by the points in "flowpoints" at "instant" and update the value corresponding to the volume of water at each point
function moving_water!( flowpoints::Vector, instant::Int64, rain::Matrix, permeability::Matrix )
    v_moving = 0 # Water flowing from already visited cells following the path
    for (rₚ, cₚ, Δhminₚ, vₚ) in flowpoints
        vₚ = (vₚ + rain[instant][rₚ, cₚ] + v_moving) * (1 - permeability[rₚ, cₚ])
        if Δhmin <= 0 # If the following cell in the path is below the current one
            # All the water on the cell flows to the next one
            v_moving = vₚ
            vₚ = 0
        else # If the cell is below all the surrounding cells
            # The water accumulates on the cell, if the volume of the water added to the height of the cell exceeds the height of the lowest adjacent cell
             # the water overflows 
            if vₚ > Δhminₚ
                v_moving = vₚ - Δhminₚ
                vₚ -= v_moving
            else
                v_moving = 0
            end
        end
    end
    return v_moving
end

function flow( flowpoints::Vector, r::Int64, c::Int64, rows::Int64, cols::Int64, instants::Int64, band::Matrix{Float32}, rev_rain_band::Matrix, permeability_band::Matrix )
    while 1 < r < rows && 1 < c < cols
        # While its still raining
        v_moving = 0 # Volume of water flowing FROM the current cell to the next in the current instant
         # "moving_water!" computes the water flowing TO the current cell and so is an instant behind accounting for "v_moving" 
        while instants > 0
            # Ad ogni istante ricalcoliamo il volume d'acqua che sta scorrendo lungo il percorso del flusso, in questo modo il volume d'acqua
             # di ogni nuova cella tiene conto del fatto che la pioggia influisce su tutte le celle appertenenti al percorso
            # Volume of water on the current cell
            v = ( sum(rain_band[instants:end][r, c]) + v_moving + moving_water!(flowpoints, instant, rev_rain_band, permeability_band) ) * (1 - permeability_band[r, c])
            dir, Δhmin = mindir(band, r, c)
            if Δhmin <= 0 # If there is a lower adjacent cell all the water flows there
                v_moving = v
                v = 0
                instants -= 1
            else # If the cell is a local minimum the water keeps accumulating untill it overflows
                while instants > 0 && Δhmin > v
                    v = ( v + rain_band[instants][r, c] + moving_water!(flowpoints, instant, rev_rain_band, permeability_band) ) * (1 - permeability[r, c])
                    instants -= 1
                end
                # If the rain stops before the water overflows the stream alts on the cell
                v_moving = instants == 0 ? 0 : v - Δhmin
                v -= v_moving
            end

            push!(flowpoints, (r, c, Δhmin, v))
            
            # If there is no water flowing the flow stops
            if v_moving == 0 && instants == 0
                break
            end

            # Flow to the next cell
            Δr, Δc = hash_adjacent(dir)
            r += Δr
            c += Δc
        end
        # After the rain stops but there is still water flowing
        while v_moving > 0
        end
    end
end




function mindir( mat, r::Int64, c::Int64, noDataValue::Real )
    dir = [0]
    min = Inf
    @inbounds for i in 1:8
        #   println( "$i) $(mat[i][r, c])" )
        if mat[r, c, i] != noDataValue
            if mat[r, c, i] < min
                dir = [i]
                min = mat[r, c, i]
                continue
            end
            if mat[r, c, i] == min
                push!(dir, i)
            end
        end
    end
    return (dir, min)
end

function flow( band::Matrix{Float32}, rain, permeability, noDataValue::Real )
    rows, cols = size(band)
    water = zeros(Float64, rows, cols, 2) # water[:,:,1]: water volume, water[:,:,2]: volume of incoming water
    flow = Array{ Union{ Vector{Float64}, Float64 } }(noDataValue, rows, cols)

    for rain_epoch in rain
        # Compute the amount of water for every cell
        water[:,:,1] = (water[:,:,1] + rain_epoch + water[:,:,2]) .* (1 .- permeability)

        # For every cell compute the water flowing from that cell to another and update the volume of water in the cell
        @inbounds for r in 1:rows, c in 1:cols
            if band[r,c] == noDataValue
                water[r, c, 1] = water[r, c, 2] = noDataValue
                flow[r, c] = noDataValue
                continue
            end 
            # Get the minimum difference in height with the neighbouring cells and the directions to the cells at that height
            dir, Δhmin = mindir(band, r, c)
            # Get the displacements to reach said cells from the current one
            flow[r, c] = hash_adjacent.(dir)

            # If the volume of water is greater than the difference in height between the current cell and the lowest adjacent, the water flows
            if water[r, c, 1] > Δhmin
                ld = length(flow[r,c])
                # The volume of flowing water depends on the difference of height between the cell and the adjacent one
                v = Δhmin <= 0 ? water[r, c, 1] : water[r, c, 1] - Δhmin
                # To the incoming water of cell "(r+Δr, c+Δc)" is added a portion of the water of cell "(r, c)" depending on the number of paths
                for delta in flow[r,c]
                    water[r+delta[1], c+delta[2], 2] += v / ld
                end
                # The remaining water is equal to the original volume less the amount that passed to the adjacent cell
                water[r, c, 1] -= v 
            end
        end
    end

    return flow, water
end


#=
Problema:
Trovare il percorso seguito dall'acqua tenendo conto che se il percorso trova un'avvallamento e si ferma li, l'acqua si accumula e con il
    passare del tempo potrebbe riuscire a fuoriuscire e proseguire la corsa.

Idea:
Trovare il percorso fino al primo stop (primo avvallamento) calcolando il cammino minimo discendente (il cammino in cui da ongi punto al successivo si scende),
    trovato il primo stop, verificare se l'acqua ad un certo punto (temporale) può straripare (se l'avvallamento è profondo come il fosso delle marianne non succederà, o,
    qunatomeno, non in tempi abbastanza brevi da avere senso per lo studio), se ciò può avvenire (perchè la profondità dell'avvallamento lo consente), riapplicare il calcolo
    del cammino minimo da questo punto fino al prossimo, ripetere finchè non si trova uno stop forte (= buco molto profondo).
Calcolato questo percorso ideale utlizzando i dati della piovosità ecc calcolare l'effettivo flusso d'acqua (il flusso seguirà per forza questo percorso essendo il
    minimo, quindi la problematica principale sarà individuare il punto di stop del flusso).
=#




function flowpath!( mat::Matrix{Float32}, band::Matrix{Float32}, noDataValue::Real, row::Int64, col::Int64, )
    if band[row, col] != noDataValue
        r, c = row, col
        mat[r, c] = 1.f0
        directions, Δhmin = mindir(band, row, col)
        while Δhmin < 0
            if lenght(directions) > 1
                for dir in directions[2:end]
                    Δr, Δc = hash_adjacent(dir)
                    flowpath!( mat, band, noDataValue, r + Δr, c + Δc )
                end
            end
            r, c .+= hash_adjacent(directions[1])
            directions, Δhmin = mindir(band, r, c)
        end
    end
end














GC.gc()
# Per porre a 1 tutte le celle di un raster che stanno dentro and un poligono delimitato da uno shapefile
using Rasters
using Shapefile

dtm_file = split( @__DIR__ , "\\Porting\\")[1] * "\\Mappe\\DTM_wgs84.tiff"
# sat_file = split( @__DIR__ , "\\Porting\\")[1] * "\\Mappe\\sat WGS84\\sette_sorelle.shp"
csoil_file = "D:\\Z_Tirocinio_Dati\\Classi suolo WGS84\\Classi suolo.shp"
perm_file = "D:\\Z_Tirocinio_Dati\\Permeabilità suolo WGS84\\Permeabilità suolo.shp"
ccs_file = "D:\\Z_Tirocinio_Dati\\ccs WGS84\\ccs.shp"
ccs_file = "C:\\Users\\DAVIDE-FAVARO\\Desktop\\Dati\\ccs WGS84\\ccs.shp"


dtm = Raster(dtm_file)
# sat_shp = Shapefile.Handle(sat_file).shapes[1]
csoil_shp = Shapefile.Table(csoil_file)
perm_shp = Shapefile.Table(perm_file)
ccs_shp = Shapefile.Table(ccs_file)

# Creation of the final raster
originX = dtm.dims[1][1]
originY = dtm.dims[2][1]
stepX = dtm.dims[1][2] - originX
stepY = dtm.dims[2][2] - originY
dims =  X( Projected( originX:stepX:dtm.dims[1][end]+stepX; order=Rasters.ForwardOrdered(), span=Rasters.Regular(stepX), sampling=Rasters.Intervals(Rasters.Start()), crs=EPSG(4326) ) ),
        Y( Projected( originY:stepY:dtm.dims[2][end]+stepY; order=Rasters.ReverseOrdered(), span=Rasters.Regular(stepY), sampling=Rasters.Intervals(Rasters.Start()), crs=EPSG(4326) ) )

#   csoil = Raster( fill(dtm.missingval, dims); missingval=dtm.missingval )
#   perm = Raster( fill(dtm.missingval, dims); missingval=dtm.missingval )
ccs = Raster( fill(dtm.missingval, dims); missingval=dtm.missingval )

# CAMBIARE IL VALORE ASSOCIATO A "X" IN missingval O SIMILI, METTERE VALORI SENSATI E SISTEMARE I VALORI DELLE COPPIE
#   cdict = Dict( val => Float32(i) for (i, val) in enumerate(unique(csoil_shp.gr_idrolog)) )
#   for (shape, value) in zip(csoil_shp.geometry, csoil_shp.gr_idrolog)
#       rasterize!( csoil, shape; fill=cdict[value], order=(X, Y) )
#   end
#   
#   Rasters.write( "D:\\Z_Tirocinio_Dati\\Classi suolo WGS84\\Classi suolo.tiff", csoil )


#   pdict = Dict( val => isnothing(tryparse(Float32, val)) ? 0.f0 : tryparse(Float32, val) for (i, val) in enumerate(unique(perm_shp.perm_clas)) )
#   for (shape, value) in zip(perm_shp.geometry, perm_shp.perm_clas)
#       rasterize!( perm, shape; fill=pdict[value], order=(X, Y) )
#   end
#   
#   Rasters.write( "D:\\Z_Tirocinio_Dati\\Permeabilità suolo WGS84\\Permeabilità suolo.tiff", perm )


ccdict = Dict( desc => Float32(i) for (i, desc) in enumerate(unique(ccs_shp.legenda)) )
for (shape, value) in zip(ccs_shp.geometry[1:2], ccs_shp.legenda[1:2])
    rasterize!( ccs, shape; fill=ccdict[value], order=(X, Y) )
end

rasterize!( ccs, ccs_shp.geometry[1]; fill=1.f0, order=(X, Y) )
rasterize!( ccs, ccs_shp.geometry[394066]; fill=3.f0, order=(X, Y) )
rasterize!( ccs, ccs_shp.geometry[2]; fill=2.f0, order=(X, Y) )

Rasters.write( "D:\\Z_Tirocinio_Dati\\test.tiff", ccs )

# I poligoni 1 e 2 sono molto vicini ma non sono esattamente adiacenti (i loro mbr si intersecano) in mezzo a loro due c'è il poligono 394066.
#   Dopo la rasterizzazione i due poligoni si toccano coprendo il poligono che dovrebbe stare tra i due.
#   Rasterizzado tutti e tre non sembrano esserci sovrapposizioni nelle celle, ma i poligoni sono di piccole dimensioni e non sono rasterizzati con precisione.


using ArchGDAL
const agd = ArchGDAL

# DA "unique(ccs)" VENGONO 58 VALORI DIVERSI, AL POSTO DI 86 (un valore per ogni poligono più uno per il "-9999.0") QUINDI DEI POLIGONI VENGONO COPERTI DA ALTRI
res_ids = tryparse.( Int64, split( read("C:\\Users\\DAVIDE-FAVARO\\Desktop\\Dati\\res_ids.txt", String), "\n", keepempty=false ) )
#   ccs_shp = agd.read("D:\\Z_Tirocinio_Dati\\ccs WGS84\\ccs.shp")
#= IL SORTING VIENE EFFETUATO PRIMA DI SALVARE GLI ID
features =  collect( agd.getlayer( agd.read("C:\\Users\\DAVIDE-FAVARO\\Desktop\\Dati\\ccs WGS84\\ccs.shp"), 0 ) )
sort!(res_ids, lt=(x, y) -> isless( agd.geomarea.( agd.getgeom.(features[[x, y]]) )... ) )
=#
for (i, id) in enumerate(resids)
    rasterize!( ccs, ccs_shp.geometry[id], fill=Float32(i), order=(X,Y) )
end
Rasters.write( "C:\\Users\\DAVIDE-FAVARO\\Desktop\\Dati\\sat_polys.tiff", ccs )
Rasters.write( "D:\\Z_Tirocinio_Dati\\sat_polys.tiff", ccs )











using Rasters
using SpatialGraphs

dtm = Raster(split( @__DIR__ , "\\Porting\\")[1] * "\\Mappe\\DTM_wgs84.tiff")
csoil = Raster("C:\\Users\\DAVIDE-FAVARO\\Desktop\\Dati\\Classi suolo WGS84\\Classi suolo.tiff")
perm = Raster("C:\\Users\\DAVIDE-FAVARO\\Desktop\\Dati\\Permeabilità suolo WGS84\\Permeabilità suolo.tiff")
perm = Raster("D:\\Z_Tirocinio_Dati\\Permeabilità suolo WGS84\\Permeabilità suolo.tiff")

perm2 = Raster(
    replace(
        perm.data,
        0.f0 => -9999.f0,
        2.f0 => 0.036f0,
        3.f0 => 0.36f0,
        4.f0 => 3.6f0,
        5.f0 => 36.f0,
        6.f0 => 360.f0
    ),
    perm.dims,
    missingval = perm.missingval
)

# "weights" sarà il raster delle resistenze di ongi cella potrebbe essere semplicemente il raster della permeabilità, oppure una combinazione di permeabilità e altezza
wrg = weightedrastergraph(
    perm2,
    directed = true,
    condition_raster = dtm,
    condition = ( hs, hd ) -> hs != dtm.missingval && hd != dtm.missingval && hs >= hd
)










function run_runoff( dem, ccs, source, target, resolution::Integer, folder::AbstractString=".\\" )

    if agd.geomdim(source) != 0
        throw(DomainError(source, "`source` must be a point"))
    end

    target_layer, landcover_layer = agd.getlayer([target, ccs], 0) 

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
    # Geometry containing both the source and the target
    # area = agd.union( bbox_src, bbox_trgt )
    # Se union non ritorna una geometria unica ma un'unica geometria composta di due elementi disgiunti
    area = agd.boundingbox(agd.union( bbox_src, bbox_trgt ))
    

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
    landcover_clip = agd.intersection(landcover, area)
    landcover_layer = agd.getlayer(landcover_clip, 0)   

    rows = agd.height(landcover_layer)
    cols = agd.width(landcover_layer)
    valNoData = -9999.0
    gtiff_driver = agd.getdriver("GTiff")


    landcover_ds = agd.create( path_temp_landcover, driver=gtiff_driver, width=rows, height=cols, nbands=1, dtype=Float32 )
    agd.setgeotransform!( landcover_ds, agd.getgeotransform(dem) )
    agd.setproj!(landcover_ds, refsys)
    band_lc = agd.getband(landcover_ds, 1)
    agd.setnodatavalue!(band_lc, valNoData)
    agd.fillraster!(band_lc, valNoData)
    bandlc = agd.read(band_lc)
    
 """ TRASFORMA IN RASTER LA PORZIONE DI `landcover` PRESA PRIMA """
    gdal.RasterizeLayer(lc_ds, [1], lc_layer,options=["ATTRIBUTE="+self.text_lcfield])
    lc_ds=None

    lc_layer=QgsRasterLayer(self.path_temp_lc,"lc_layer")
 """"""
    landcover_layer = agd.gdalrasterize( x -> x, landcover_ds )



    if soil_text == "Valore campo"
        soil_control = 1

        agd.getlayer(path_lc, 0)

        soil_ds = agd.create( path_temp_soil, driver=gtiff_driver, width=round(Int64, x_res), height=round(Int64, y_res), nbands=1, dtype=Float32 )
        agd.setgeotransform!( soil_ds, agd.getgeotransform(dem) )
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





 # VANNO PASSATI A flow I RANGE DELL'AREA IN CUI PIOVE 

 # ASSUMENDO CHE target VENGA PASSATO COME FILE VETTORIALE
    target_layer = collect(agd.getlayer(target, 0))
    target_geom = agd.getgeom(agd.getgeom(target_layer[1], 0), 0)
 # SE target VIENE PASSATO COME GEOMETRIA ELIMINARE LE DUE RIGHE SOPRA
    length = agd.ngeom(target_geom)
    target_points = [ agd.getpoint(target_geom, i)[1:2] for i in 1:length ]

    x_source = agd.getx(source, 0)
    y_source = agd.gety(source, 0)
    mat = agd.read(".\\connectivity.tiff")
    res = flow(mat, x_source, y_source)

    volumes = []
    for point in res
        lc = landcover_layer[point...]
        if soil_control == 1
            soil = soil_layer[point...]
        else
            soil = text_soil
        end
        # IL BLOCCO TRY CATCH NON LO VEDO PARTICOLARMENTE NECESSARIO
        try
         # DA MODIFICARE TENENDO CONTO DEI VETTORIALI DEL SUOLO
         # PER LA CLASSE DEL SUOLO SI PUO' FAR RIFERIMENTO AL DIZIONARIO DEL FILE
         # CREDO STIA PRENDENDO LA RESISTENZA DEL SUOLO NEL PUNTO CORRENTE
            # Se passiamo come parameto un raster delle resisteze questa riga va sostituita
            cn = listaclc[lc][ classisoil[soil] ] # Classisoi ottenuto da "substance.db"
            S = 254.0((100 / cn) - 1)
        catch
            S = 0
        end
        pcheck = (p0 - 0.2S)^2 / (p0 - 0.2S + S)
        push!(volumes, pcheck > 0.2S ? pcheck : 0 )
    end



end



end # module



#=
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
               clc = lc_layer[x, y]
               if soil_control == 1
                   r, c = toIndexes(soil_layer, x, y)
                   soil = soil_layer[r, c]
               else
                   soil = text_soil
               end
               try
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
=#