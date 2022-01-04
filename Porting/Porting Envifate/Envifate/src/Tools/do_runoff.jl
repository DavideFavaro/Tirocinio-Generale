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



import ArchGDAL as agd

# VEDERE @view, @inbound, @turbo, @fast, LoopedVectorization.jl e StableArrays.jl PER ULTERIORI OTTIMIZZAZIONI


# 3 DIMENSIONAL MATRIX
function connectivity_batch!( mat, dem_band::AbstractArray{T}, noDataValue::Real ) where {T <: Number}
    rows, cols = size(dem_band)
    indexes = Dict(
        (-1, -1) => 1,
        (-1, 0)  => 2,
        (-1, 1)  => 3,
        (0, -1)  => 4,
        (0, 1)   => 5,
        (1, -1)  => 6,
        (1, 0)   => 7,
        (1, 1)   => 8
    )
    # For each cell of the dem's band
    @inbounds for r in 1:rows, c in 1:cols
        if dem_band[r, c] != noDataValue
            # Indexes of adjacent cells
            for (i, j) in keys(indexes)
                if ( r+i >= 1 && r+i <= rows ) && ( c+j >= 1 && c+j <= cols ) && dem_band[r+i, c+j] != noDataValue
                    if mat[ r, c, indexes[(i,j)] ] == noDataValue
                        mat[ r, c, indexes[(i, j)] ] = dem_band[r, c] - dem_band[r+i, c+j]
                    end
                    if mat[ r+i, c+j, indexes[(-i, -j)] ] == noDataValue
                        mat[ r+i, c+j, indexes[(-i, -j)] ] = dem_band[r+i, c+j] - dem_band[r, c]
                    end
                end
            end
        end
    end
end

function connectivity( dem_band::Matrix{T}, batch_size::Integer, noDataValue::Real ) where {T <: Number}
    rows, cols = size(dem_band)
    n, m = ceil.( Int64, [rows, cols] ./ (batch_size - 1) )
    mat = fill( convert(Float32, noDataValue) , rows, cols, 8 )
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

@time mat = connectivity( band_mat, 1024, ndv )
dims = size(mat)
flattened = reshape(mat, dims[1], :, 1)[:, :, 1]
io = open("D:\\Connectivity Data\\connectivity.txt", "w")
write(io, flattened)
close(io)


using BenchmarkTools

@benchmark mat = connectivity( test1b, 256, ndv )
@benchmark mat = connectivity( test1b, 1024, ndv )
@benchmark mat = connectivity( test1b, 2048, ndv )


#=
function direct_connectivity_batch!( mat, dem_band::AbstractArray{T}, noDataValue::Real ) where {T <: Number}
    rows, cols = size(dem_band)
    indexes = Dict(
        (-1, -1) => 1,
        (-1, 0)  => 2,
        (-1, 1)  => 3,
        (0, -1)  => 4,
        (0, 1)   => 5,
        (1, -1)  => 6,
        (1, 0)   => 7,
        (1, 1)   => 8
    )
    # For each cell of the dem's band
    @inbounds for r in 1:rows, c in 1:cols
        if dem_band[r, c] != noDataValue
            # Indexes of adjacent cells
            for (i, j) in keys(indexes)
                if ( r+i >= 1 && r+i <= rows ) && ( c+j >= 1 && c+j <= cols ) && dem_band[r+i, c+j] != noDataValue
                    if dem_band[r, c] > dem_band[r+i, c+j] && mat[ r, c, indexes[(i,j)] ] == noDataValue
                        mat[ r, c, indexes[(i, j)] ] = dem_band[r, c] - dem_band[r+i, c+j]
                    end
                    if dem_band[r, c] < dem_band[r+i, c+j] && mat[ r+i, c+j, indexes[(-i, -j)] ] == noDataValue
                        mat[ r+i, c+j, indexes[(-i, -j)] ] = dem_band[r+i, c+j] - dem_band[r, c]
                    end
                end
            end
        end
    end
end

function direct_connectivity( dem_band::Matrix{T}, batch_size::Integer, noDataValue::Real ) where {T <: Number}
    rows, cols = size(dem_band)
    n, m = ceil.( Int64, [rows, cols] ./ (batch_size - 1) )
    mat = fill( convert(Float32, noDataValue) , rows, cols, 8 )
    for i in 1:n, j in 1:m
        # Find the starting and ending indexes for the current slice of the matrix
         # if the batch is one of the ending ones its ending index will be the size of the matrix for that dimension
        rows_range::UnitRange{Int64} = ( (batch_size - 1) * (i - 1) + 1 ) : ( i != n ? (batch_size - 1) * i + 1 : rows )
        cols_range::UnitRange{Int64} = ( (batch_size - 1) * (j - 1) + 1 ) : ( j != m ? (batch_size - 1) * j + 1 : cols )
        # Use the indexes to run the function on a view of the matrix (passing also the corresponding view of the dem)
        direct_connectivity_batch!( view(mat, rows_range, cols_range, :), view(dem_band, rows_range, cols_range), noDataValue )
    end
    return mat
end

@time matd = direct_connectivity( band_mat, 1024, ndv )
=#








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


@time mat = connectivity(band_mat, 2048, ndv)
df = DataFrame(mat, :auto)
CSV.write("D:\\Connectivity Matrix\\direct_connectivity.csv", df)
CSV.write("C:\\Users\\DAVIDE-FAVARO\\Desktop\\connectivity.csv", df)








import ArchGDAL as agd


include("../Library/Functions.jl")



mat = connectivity(band_mat, 2048, ndv)

rows, cols = size(mat)
dtm2 = agd.read(dtm_file)
target_ds = agd.create( "C:\\Users\\DAVIDE-FAVARO\\Desktop\\connectivity2.tiff", driver=agd.getdriver("GTiff"), width=rows, height=cols, nbands=8, dtype=Float32 )
#   agd.setgeotransform!( target_ds, [ minX, a, 0.0, maxY, 0.0, b ] )
agd.setgeotransform!( target_ds, agd.getgeotransform(dtm2) )
agd.setproj!( target_ds, agd.getproj(dtm2) )
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
                band_mats[ index_dict[(i, j)] ][r, c] = val
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








# "flow" E "mindir" CHE TENGONO CONTO DELLA POSSIBILITA' CHE CI SIANO PIU' CAMMINI MINIMI
function mindir( mat, r::Int64, c::Int64 )
    dir = [0]
    min = Inf

    @inbounds for i in 1:8
        #   println( "$i) $(mat[i][r, c])" )
        if mat[r, c, i] != -9999.0
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

function multi_flow( source_volumes, x::Real, y::Real, x_interval, y_interval, file::AbstractString )
    delta_dict = Dict(
        1 => (-1, -1),
        2 => (-1, 0),
        3 => (-1, 1),
        4 => (0, -1),
        5 => (0, 1),
        6 => (1, -1),
        7 => (1, 0),
        8 => (1, 1)
    )

    mat = agd.read(file)
    bands = [ agd.getband(mat, i) for i in 1:8 ]
    rows, cols = size(bands[1])
    flowpoints = Vector{Tuple{Int64, Int64}}()
    r, c = Functions.toIndexes(mat, x, y)
    r_min, c_min = Funtions.toIndexes( x_interval[1], y_interval[1] )
    r_max, c_max = Funtions.toIndexes( x_interval[2], y_interval[2] )

    x₀ = source_volumes[1]
    push!(flowpoints, (r, c, x₀))
    val = -Inf
    for rain_i in source_volumes[2:end]
        if r < r_min || r > r_max || c < c_min || c > c_max || # If outside the rain zone
           r < 1 || r > rows || c < 1 || c > cols ||           # If outside the raster
           val > 0 || x₀ <= 0 ||                               # If the flow stops
            break
        end 
        # flow from (r, c) to (r+Δr, c+Δc)
        dir, val = mindir(bamds, r, c)
        for d in dir
            Δr, Δc = delta_dict[d]
            # NON SO COME SI OTTENGONO I VALORI PER L'ACQUA CHE SCENDE ("desc") E QUELLA CHE ESCE ("out")
            x₁ = x₀ - desc - out + rain_i
            x₀ = x₁
            push!(flowpoints, (r+Δr, c+Δc, x₀))
        end
    end
    # If the cicle ends because the flow exits the rain zone or the rainfall stops
    while 1 < r < rows && 1 < c < cols && val < 0 && x₀ > 0
        # Flow from (r, c) to (r+Δr, c+Δc)
        dir, val = mindir(bands, r, c)
        for d in dir
            Δr, Δc = delta_dict[d]
            # NON SO COME SI OTTENGONO I VALORI PER L'ACQUA CHE SCENDE ("desc") E QUELLA CHE ESCE ("out")
            x₁ = x₀ - desc - out
            x₀ = x₁
            push!(flowpoints, (r+Δr, c+Δc, x₀))
        end
    end
    return flowpoints
end
















# Find the direction with the greatest difference in height from (r, c)
function mindir( mat, r::Int64, c::Int64 )
    dir = 0
    min = Inf
    @inbounds for i in 1:8
        #   println( "$i) $(mat[i][r, c])" )
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

# "source_volumes" sarà un vettore contenente la quantità di pioggia caduta all'istante "i", si assume che la pioggia cada in modo uniforme in tutta l'area
 # interessata dal fenomeno
# "x_interval" contiene i valori minimo e massimo della dimensione "x" dell'area interessata dalla pioggia
# Similmente "y_interval" contiene i valori su "y"
function flow( source_volumes, x::Real, y::Real, x_interval, y_interval, file::AbstractString )
    delta_dict = Dict(
        1 => (-1, -1),
        2 => (-1, 0),
        3 => (-1, 1),
        4 => (0, -1),
        5 => (0, 1),
        6 => (1, -1),
        7 => (1, 0),
        8 => (1, 1)
    )

    mat = agd.read(file)
    bands = [ agd.getband(mat, i) for i in 1:8 ]
    rows, cols = size(bands[1])
    flowpoints = Vector{Tuple{Int64, Int64}}()
    r, c = Functions.toIndexes(mat, x, y)
    r_min, c_min = Funtions.toIndexes( x_interval[1], y_interval[1] )
    r_max, c_max = Funtions.toIndexes( x_interval[2], y_interval[2] )

    x₀ = source_volumes[1]
    push!(flowpoints, (r, c, x₀))
    val = -Inf
    for rain_i in source_volumes[2:end]
        if r < r_min || r > r_max || c < c_min || c > c_max || # If outside the rain zone
           r < 1 || r > rows || c < 1 || c > cols ||           # If outside the raster
           val > 0 || x₀ <= 0 ||                               # If the flow stops
            break
        end 
        # flow from (r, c) to (r+Δr, c+Δc)
        dir, val = mindir(bamds, r, c)
        Δr, Δc = delta_dict[dir]
        r += Δr
        c += Δc

        # NON SO COME SI OTTENGONO I VALORI PER L'ACQUA CHE SCENDE ("desc") E QUELLA CHE ESCE ("out")
        x₁ = x₀ - desc - out + rain_i
        x₀ = x₁
        push!(flowpoints, (r, c, x₀))
    end
    # If the cicle ends because the flow exits the rain zone or the rainfall stops
    while 1 < r < rows && 1 < c < cols && val < 0 && x₀ > 0
        dir, val = mindir(bands, r, c)
        #   println()
        #   println(val)
        #   println()
        #   println()
        
        # Flow from (r, c) to (r+Δr, c+Δc)
        Δr, Δc = delta_dict[dir]
        r += Δr
        c += Δc
        # NON SO COME SI OTTENGONO I VALORI PER L'ACQUA CHE SCENDE ("desc") E QUELLA CHE ESCE ("out")
        x₁ = x₀ - desc - out
        x₀ = x₁
        push!(flowpoints, (r, c, x₀))
    end
    return flowpoints
end

# "flow" to apply whene there is a singular source
# NON SO SE SI DEBBA SOMMARE AD OGNI PASSO IL VALORE INZIALE O SE SI DEBBA CONSIDERARE COMUNQUE IL FATTORE TEMPORALE (SE NELLA FONTE ARRIVA ACQUA PER UN CERTO TEMPO)
function flow( source_volume::Real, x::Real, y::Real, file::AbstractString )
    delta_dict = Dict(
        1 => (-1, -1),
        2 => (-1, 0),
        3 => (-1, 1),
        4 => (0, -1),
        5 => (0, 1),
        6 => (1, -1),
        7 => (1, 0),
        8 => (1, 1)
    )

    mat = agd.read(file)
    bands = [ agd.getband(mat, i) for i in 1:8 ]
    rows, cols = size(bands[1])
    flowpoints = Vector{Tuple{Int64, Int64}}()
    r, c = Functions.toIndexes(mat, x, y)
    x₀ = source_volume
    push!(flowpoints, (r, c, x₀))
    val = -Inf
    while 1 < r < rows && 1 < c < cols && val < 0 && x₀ > 0
        dir, val = mindir(bands, r, c)
        #   println()
        #   println(val)
        #   println()
        #   println()
        
        # Flow from (r, c) to (r+Δr, c+Δc)
        Δr, Δc = delta_dict[dir]
        r += Δr
        c += Δc
        # NON SO COME SI OTTENGONO I VALORI PER L'ACQUA CHE SCENDE ("desc") E QUELLA CHE ESCE ("out")
        x₁ = x₀ - desc - out
        x₀ = x₁
        push!(flowpoints, (r, c, x₀))
    end
    return flowpoints
end


#   x = 726467.4299990014
#   y = 5.025981399455068e6
#   path = "C:\\Users\\DAVIDE-FAVARO\\Desktop\\Connectivity Data\\connectivity.tiff"
#
#   res = flow( x, y, path )











function flow( x::Real, y::Real, band::Matrix{Float32}, rain_band::Matrix{Float32}, permeability_band::Matrix{Float32} )
    delta_dict = Dict(
        1 => (-1, -1),
        2 => (-1, 0),
        3 => (-1, 1),
        4 => (0, -1),
        5 => (0, 1),
        6 => (1, -1),
        7 => (1, 0),
        8 => (1, 1)
    )

    rows, cols = size(band)
    flowpoints = Vector{Tuple{Int64, Int64}}()
    r, c = Functions.toIndexes(mat, x, y)
    x = rain_band[r, c]
    val = -Inf
    while 1 < r < rows && 1 < c < cols && val < 0
        dir, val = mindir(bands, r, c)
        #   println()
        #   println(val)
        #   println()
        #   println()
        
        # Flow from (r, c) to (r+Δr, c+Δc)
        Δr, Δc = delta_dict[dir]
        r += Δr
        c += Δc
        push!(flowpoints, (r, c, x))
        x = x + rain_band[r, c] - ( x * permeability_band[r, c] ) - ( val ) # La quantità d'acqua in uscita dovrebbe dipendere dalla differenza in altezza delle celle  
    end
    return flowpoints
end


















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