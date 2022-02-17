module Lights



import ArchGDAL as agd

include("..\\Library\\Functions.jl")



function run_light( dem, source, resolution::Integer, intensity::Real, source_height::Real=0.0, observer_height::Real=1.75, rarefraction::Real=0.14286;
                    output_path::AbstractString=".\\light_intensity.tiff"  )
 # *VERSIONE CON INTENSITA'/FONTI MULTIPLE
 #= *
    if any( i -> i < 0, intensity )
        throw(DomainError(intenisty, "`intenisty` must be greater than 0."))
    end
 =#
    if intensity < 0
        throw(DomainError(intenisty, "`intenisty` must be greater than 0."))
    end

    src_layer = agd.getlayer(source, 0)
 # * src_geoms = agd.getgeom.(collect(src_layer))
    src_geom = agd.getgeom(collect(src_layer)[1])

 #= *
    if any( geom -> agd.geomdim(geom) != 0, src_geoms )
        throw(DomainError(source, "`source` must be a point."))
    end
 =#
    if agd.geomdim(src_geom) != 0
        throw(DomainError(source, "`source` must be a point."))
    end

    refsys = agd.getspatialref(src_layer)

    if agd.importWKT(agd.getproj(dem)) != refsys
        throw(DomainError("The reference systems are not uniform. Aborting analysis."))
    end

 #= *
    nfeature = 0
    features = agd.getgeom.( agd.getlayer(source, 0) )
    for feature in features
        x_source = agd.getx(feature, 0)
        y_source = agd.gety(feature, 0)
        nfeature += 1
        # The assumption is that the lightsource can only project in a downward emisphere (like a streetlight) and so all terrain that is higher than the source will
        #   block the light for the terrain beyond.
        vis_mat = Viewshed.viewshed( dem, x_source, y_source, source_height )
    end
 =#
    x_source = agd.getx(src_geom, 0)
    y_source = agd.gety(src_geom, 0)
    # The assumption is that the lightsource can only project in a downward emisphere (like a streetlight) and so all terrain that is higher than the source will
    #   block the light for the terrain beyond.
    vis_mat = Viewshed.viewshed( dem, x_source, y_source, source_height )
    rows, cols = size(vis_mat)
    noDataValue = -9999.f0
    data = fill(noDataValue, rows, cols)
    @inbounds for r in 1:rows, c in 1:cols
        if dem[r, c] != noDataValue
            if !vis_mat[r, c]
                data[r, c] = 0.f0
            else
                # I₂ = ( d₁^2 / d₂^2 ) * I₁
                #  d₁ = 1.0 -> I₂ = I₁ * ( 1 / d₂^2 ) = I₁ / d₂^2
                data[r, c] = intensity / Viewshed.distance( toCoords(dem, r, c), (x_source, y_source) )^2
            end    
        end
    end
    Functions.writeRaster( data, agd.getdriver("GTiff"), agd.getgeotransform(dem), resolution, refsys, noDataValue, output_path, false )
end



end # module




# ------------------------------------------------------------------------------- TESTING ---------------------------------------------------------------------------------------------

#= CON "GeoArrays.jl"
    using GeoArrays
    const ga = GeoArrays

    rotate_x( xp::In64, yp::In64, xc::In64, yc::In64, θ::In64 ) = round( Int64, (xp - xc)cos(deg2rad(θ)) - (yp - yc)sin(deg2rad(θ)) + xc )
    rotate_y( xp::In64, yp::In64, xc::In64, yc::In64, θ::In64 ) = round( Int64, (xp - xc)sin(deg2rad(θ)) + (yp - yc)cos(deg2rad(θ)) + yc )
    rotate_point( xp::In64, yp::In64, xc::In64, yc::In64, θ::In64 )::Tuple{Int64, Int64} = θ == 0 ? (xp, yp) : ( rotate_x( xp, yp, xc, yc, θ ), rotate_y( xp, yp, xc, yc, θ ) )
    distance( x0::AbstractFloat, y0::AbstractFloat, x1::AbstractFloat, y1::AbstractFloat ) = √( ( x1 - x0 )^2 + ( y1 - y0 )^2 )
    distance( p0::Tuple{T, T}, p1::Tuple{T, T} ) where {T <: AbstractFloat} = √( ( p1[1] - p0[1] )^2 + ( p1[2] - p0[2] )^2 )

    function toCoords(dem, r, c)
        return Tuple{Float64, Float64}(ga.coords(dem, [r, c]))
    end

    function toIndexes(dem, x, y)
        return Tuple{Int64, Int64}(ga.indices(dem, [x, y]))
    end

    # VERSIONE CON "GeoArrays"
    function viewshed( dtm::GeoArray{Float32}, x0::Float64, y0::Float64, h0::Float32 )
        # Source cell
        r0, c0 = toIndexes(dtm, x0, y0)
        # Final point of the right horizontal profile 
        rl, cl = size(dtm)
        # Total height of the source accounting for terrain
        th0 = dtm[r0, c0][1] + h0
        # Visibility matrix
        vizmat = falses(size(dtm)...)
        # The source is always visible
        vizmat[r0, c0] = true
        # Check the visibility along 360 lines radially expanding from the source at an angle of 1° between two adjacent.
        for α in 1:89, β in [0, 90, 180, 270]
            # Final point of the profile of the line with agle `α + β`
            rn, cn = rotate_point( rl, c0, r0, c0, α + β )
            Δr = rn - r0
            Δc = cn - c0
            # Values to add to row and column of a cell on the line to pass to another cell of the line
            r_inc, c_inc = (Δr, Δc) ./ max( abs(Δr), abs(Δc) )
            # Indexes of the first cell after the source cell
            r = r0 + r_inc
            c = c0 + c_inc
            r1 = round(Int64, r)
            c1 = round(Int64, c)
            # The first cell after the source is always visible
            vizmat[r1, c1] = true
            # Slope between the source and the first cell
            slope = ( dtm[r1, c1] - dtm[r0, c0] ) / distance( toCoords(dtm, r1, c1), (0.0, 0.0) )
            # Iterate over each cell of the profile on the line
            while r <= rl && r <= rn && c <= cl && c <= cl
                # Calculate the precise indexes of the cell
                rint, cint = round.(Int64, [r, c])
                # Skip cell already known to be visible
                !vizmat[rint, cint] && continue
                # Compute the slope of the new cell 
                new_slope = (dtm[rint, cint] - dtm[r0, c0]) / distance( toCoords(dtm, r1, c1), (0.0, 0.0) )
                # If the new slope is greater than the original one the cell is visible
                if new_slope >= slope
                    vizmat[rint, cint] = true
                    slope = new_slope
                end
                # If the cell is higher than the source al cell beyond the current ne will be hidden
                dtm[rint, cint][1] > th0 && break
                r += r_inc
                c += c_inc
            end
        end
        return vizmat
    end

    function run_light( dem, source::Tuple{Float64, Float64}, resolution::Integer, intensity::Real, source_height::Real=0.0, observer_height::Real=1.75, rarefraction::Real=0.14286;
                        output_path::AbstractString=".\\light_intensity.tiff"  )
        if intensity < 0
            throw(DomainError(intenisty, "`intenisty` must be greater than 0."))
        end

        x_source, y_source = source
        # The assumption is that the lightsource can only project in a downward emisphere (like a streetlight) and so all terrain that is higher than the source will
        #   block the light for the terrain beyond.
        vis_mat = viewshed( dem, x_source, y_source, source_height )
        rows, cols = size(vis_mat)
        noDataValue = -9999.f0
        data = fill(noDataValue, rows, cols)
        @inbounds for r in 1:rows, c in 1:cols
            if dem[r, c] != noDataValue
                if !vis_mat[r, c]
                    data[r, c] = 0.f0
                else
                    # I₂ = ( d₁^2 / d₂^2 ) * I₁
                    #  d₁ = 1.0 -> I₂ = I₁ * ( 1 / d₂^2 ) = I₁ / d₂^2
                    data[r, c] = intensity / distance( toCoords(dem, r, c), (x_source, y_source) )^2
                end
            end
        end
        return data
        #   Functions.writeRaster( data, agd.getdriver("GTiff"), agd.getgeotransform(dem), resolution, refsys, noDataValue, output_path, false )
    end

    dtmga = ga.read(split( @__DIR__ , "\\Porting\\")[1] * "\\Mappe\\DTM_32.tiff")
    #   dtm = ga.read(split( @__DIR__ , "\\Porting\\")[1] * "\\Mappe\\DTM_wgs84.tiff")
    src = (726454.9302346368, 5.025993899219433e6)
    #   src = (11.930065824163105,45.425861311724816)
    run_light( dtm, src, 25, 15.0, 1.0f0, output_path="C:\\Users\\DAVIDE-FAVARO\\Desktop\\test_light.tiff" )
=#


# CON "Rasters"



#   QUESTE FORMULE NON TENGONO CONTO CHE LA "Y" E' DECRESCENTE, MA FORSE NON CAMBIA PERCHE' L'ANGOLO E' LO STESSO MA SPOSTATO/INVERTITO ?
#       forse ho sistemato
rotate_x( xp::Float64, yp::Float64, xc::Float64, yc::Float64, θ::Int64 ) = xc + ( (xp - xc)cos(deg2rad(θ)) - (yp - yc)sin(deg2rad(θ)) )
rotate_y( xp::Float64, yp::Float64, xc::Float64, yc::Float64, θ::Int64 ) = yc - ( (xp - xc)sin(deg2rad(θ)) + (yp - yc)cos(deg2rad(θ)) )
rotate_point( xp::Float64, yp::Float64, xc::Float64, yc::Float64, θ::Int64 ) = θ == 0 ? (xp, yp) : ( rotate_x( xp, yp, xc, yc, θ ), rotate_y( xp, yp, xc, yc, θ ) )
distance( x0::Float64, y0::Float64, x1::Float64, y1::Float64 ) = √( ( x1 - x0 )^2 + ( y1 - y0 )^2 )
distance( p0::Tuple{T, T}, p1::Tuple{T, T} ) where {T <: Float64} = √( ( p1[1] - p0[1] )^2 + ( p1[2] - p0[2] )^2 )
Base.getindex( collection::Raster, x::Float64, y::Float64 ) = collection[X(Near(x)), Y(Near(y))][1]
Base.setindex!( collection::Raster, v, x::Float64, y::Float64, ) = collection[X(Near(x)), Y(Near(y))] .= v
condition( x::Float64, xn::Float64, y::Float64, yn::Float64, θ::Int64 ) = θ == 0 ? (x <= xn && y >= yn) :
                                                                              θ == 90 ? (x >= xn && y >= yn) :
                                                                                  θ == 180 ? (x >= xn && y <= yn) :
                                                                                      θ == 270 ? (x <= xn && y <= yn) : throw(DomainError(β, "`β` must either be 0, 90, 180 or 270."))

function viewshed( dtm::Raster{Float32, 3}, x0::Float64, y0::Float64, h0::Float32 )
    # Dimensions of the raster.
        # N.B. Y axis is in decreasing order.
    xmin = dtm.dims[1][1]
    xmax = dtm.dims[1][end]
    ymax = dtm.dims[2][1]
    ymin = dtm.dims[2][end]
    # Displacement from a cell to the next alogn the two axis.
    xstep = dtm.dims[1][2] - xmin
    ystep = dtm.dims[2][2] - ymin
    # Total height of the source accounting for terrain.
    th0 = dtm[x0, y0] + h0
    # Visibility matrix.
    vizmat = Raster(falses(size(dtm)...), dtm.dims)
    # The source is always visible.
    vizmat[x0, y0] = true
    # Check the visibility along 360 lines radially expanding from the source at an angle of 1° between two adjacent ones.
        # Iterating over the four quadrants and then on the angles in each of them makes easier to check the coordinates.
    for α in [0, 90, 180, 270], β in 0:89
        # Final point of the profile of the line with agle `α + β`.
        xn, yn = rotate_point(xmax, y0, x0, y0, α + β)
        # Diplacement between the center and the final point along the two axis.
        Δx = xn - x0
        Δy = yn - y0
        # Values to add to row and column of a cell on the line to pass to another cell of the line.
        x_inc, y_inc = ( (Δx, Δy) ./ max(abs(Δx), abs(Δy)) ) .* (xstep, ystep)
        # Coordinates of the first cell after the source cell.
        x = x0 + x_inc
        y = y0 + y_inc
        # The first cell after the source is always visible.
        vizmat[x, y] = true
        # Slope between the source and the first cell.
        slope = (dtm[x, y] - dtm[x0, y0]) / distance((x, y), (0.0, 0.0))
        # Go to the following cell in the "profile".
        x += x_inc
        y += y_inc
        # Iterate over each cell of the profile on the line.
        while xmin <= x <= xmax && ymin <= y <= ymax && condition(x, xn, y, yn, α)
            # Skip cells already known to be visible.
            vizmat[x, y] && continue
            # Compute the slope of the new cell.
            new_slope = (dtm[x, y] - dtm[x0, y0]) / distance((x, y), (0.0, 0.0))
            # If the new slope is greater than the original one the cell is visible.
            if new_slope >= slope
                vizmat[x, y] = true
                slope = new_slope
            end
            # If the cell is higher than the source, all cells beyond it will be hidden.
            dtm[x, y] > th0 && break
            # Go to the following cell.
            x += x_inc
            y += y_inc
        end
    end
    return vizmat
end

function run_light( dem::Raster{Float32, 3}, source::Tuple{Float64, Float64}, resolution::Integer, intensity::Real, source_height::Real=0.0, observer_height::Real=1.75, rarefraction::Real=0.14286;
                    output_path::AbstractString=".\\light_intensity.tiff"  )
    if intensity < 0
        throw(DomainError(intenisty, "`intenisty` must be greater than 0."))
    end

    x_source, y_source = source
    # The assumption is that the lightsource can only project in a downward emisphere (like a streetlight) and so all terrain that is higher than the source will
    #   block the light for the terrain beyond.
    vis_mat = viewshed(dem, x_source, y_source, source_height)
    rows, cols = size(vis_mat)
    noDataValue = -9999.f0
    data = fill(noDataValue, rows, cols)
    @inbounds for r in 1:rows, c in 1:cols
        if dem[r, c] != noDataValue
            if !vis_mat[r, c]
                data[r, c] = 0.0f0
            else
                # I₂ = ( d₁^2 / d₂^2 ) * I₁
                #  d₁ = 1.0 -> I₂ = I₁ * ( 1 / d₂^2 ) = I₁ / d₂^2
                data[r, c] = intensity / distance( toCoords(dem, r, c), (x_source, y_source) )^2
            end
        end
    end
    return data
    #   Functions.writeRaster( data, agd.getdriver("GTiff"), agd.getgeotransform(dem), resolution, refsys, noDataValue, output_path, false )
end





function run_light( dem::Raster{Float32, 3}, intensity::Float32, x0::Float64, y0::Float64, h0::Float32 )
    # Dimensions of the raster.
        # N.B. Y axis is in decreasing order.
    xmin = dem.dims[1][1]
    xmax = dem.dims[1][end]
    ymax = dem.dims[2][1]
    ymin = dem.dims[2][end]
    # Displacement from a cell to the next alogn the two axis.
    xstep = dem.dims[1][2] - xmin
    ystep = dem.dims[2][2] - ymin
    # Total height of the source accounting for terrain.
    th0 = dem[x0, y0] + h0
    # Light intensity matrix.
    data = Raster(fill(dem.missingval, size(dem)...), dem.dims)
    # The intensity of the light at the source is given
    data[x0, y0] = intensity
    # Check the visibility along 360 lines radially expanding from the source at an angle of 1° between two adjacent ones.
        # Iterating over the four quadrants and then on the angles in each of them makes easier to check the coordinates.
    for α in [0, 90, 180, 270], β in 0:89
        # Final point of the profile of the line with agle `α + β`.
        xn, yn = rotate_point(xmax, y0, x0, y0, α + β)
        # Diplacement between the center and the final point along the two axis.
        Δx = xn - x0
        Δy = yn - y0
        # Values to add to row and column of a cell on the line to pass to another cell of the line.
        x_inc, y_inc = ( (Δx, Δy) ./ max(abs(Δx), abs(Δy)) ) .* (xstep, ystep)
        # Coordinates of the first cell after the source cell.
        x = x0 + x_inc
        y = y0 + y_inc
        # Distance between source and first cell of the profile
        dist = distance((x, y), (0.0, 0.0))
        # Slope between the source and the first cell.
        slope = (dem[x, y] - dem[x0, y0]) / dist
        # Value of light intensity at the cell of coordinates (`x`, `y`).
            # I₂ = I₁ * ( d₁^2 / d₂^2 )  ->  d₁ = 1.0  ->  I₂ = I₁ * ( 1 / d₂^2 ) = I₁ / d₂^2.
        data[x, y] = intensity / dist^2
        # Go to the following cell in the "profile".
        x += x_inc
        y += y_inc
        # Iterate over each cell of the profile on the line.
        while xmin <= x <= xmax && ymin <= y <= ymax && condition(x, xn, y, yn, α)
            # Skip cells already known to be visible.
            ( ( data[x, y] != 0.f0 && data[x, y] != dem.missingval ) || dem[x, y] == dem.missingval ) && continue
            # Distance of the cell from the source.
            dist = distance((x, y), (0.0, 0.0))
            # Compute the slope of the new cell.
            new_slope = (dem[x, y] - dem[x0, y0]) / dist
            # If the new slope is greater than the original one the cell is visible.
            if new_slope >= slope
                data[x, y] = intensity / dist^2
                slope = new_slope
            else
                data[x, y] = 0.0f0
            end
            # If the cell is higher than the source, all cells beyond it will be hidden.
            dem[x, y] > th0 && break
            # Go to the following cell.
            x += x_inc
            y += y_inc
        end
    end
    return data
end














using Plots
using Rasters
dtmr = Raster(split( @__DIR__ , "\\Porting\\")[1] * "\\Mappe\\DTM_32.tiff")
src = (726454.9302346368, 5.025993899219433e6)

@code_warntype viewshed(dtmr, src[1], src[2], 1.0f0)
@time vizmat = viewshed(dtmr, src[1], src[2], 1.0f0)

@code_warntype run_light(dtmr, 15.0f0, src[1], src[2], 0.0f0)
@time intmat = run_light(dtmr, 15.0f0, src[1], src[2], 0.0f0)

@code_warntype run_light(dtmr, src, 25, 15.0, 1.0f0, output_path="C:\\Users\\DAVIDE-FAVARO\\Desktop\\test_light.tiff")
intmat = run_light(dtmr, src, 25, 15.0, 1.0f0, output_path="C:\\Users\\DAVIDE-FAVARO\\Desktop\\test_light.tiff")











# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



function run_light( dem, source, resolution::Integer, intenisty::Real, source_height::Real=0.0, observer_height::Real=1.75, rarefraction::Real=0.14286,
                    output_path::AbstractString=".\\light_intensity.tiff"  )

 """ CONTROLLO SUL RASTER dem
    if not self.dem.isValid():
        QMessageBox.warning(self,"Warning", "The dem file is not valid" )
        return
 """
    if any( i -> i < 0, intensity )
        throw(DomainError(intenisty, "`intenisty` must be greater than 0."))
    end

    src_layer = agd.getlayer(source, 0)
    src_geoms = agd.getgeom.(collect(src_layer))

    if any( geom -> agd.geomdim(geom) != 0, src_geoms )
        throw(DomainError(source, "`source` must be a point."))
    end

    refsys = agd.getspatialref(src_layer)

    if agd.importWKT(agd.getproj(dem)) != refsys
        throw(DomainError("The reference systems are not uniform. Aborting analysis."))
    end

    nfeature = 0
    features = agd.getgeom.( agd.getlayer(source, 0) )
    for feature in features
        x_source = agd.getx(feature, 0)
        y_source = agd.gety(feature, 0)


     # LA/E INTENSITA' LA/E INSERIAMO COME PARAMETRO (SINGOLO VALORE/ARRAY)
        #   idxlevel = self.source.fields().indexFromName('level')
        #   intensity=feature.attributes()[idxlevel]

        nfeature += 1

        #   start_time = time.time()

        # The assumption is that the lightsource can only project in a downward emisphere (like a streetlight) and so all terrain that is higher than the source will
        #   block the light for the terrain beyond.
        profiles = Viewshed.generate_profiles(dem, x_source, y_source, source_height)
        viewshed = Viewshed.viewshed(profiles)





     """ SERVE LA VIEWSHED COME RASTER CREDO """
        viewshed=QgsProject.instance().mapLayersByName(nameviewshed)[0]
     """                                     """

     # QUESTO O QUALCHE ALTRO MODO PER OTTENERE I LIMITI DEL RASTER VIEWSHED
        # Define the area of the `dtm` that coincides with the area of the `viewshed`
        minX, maxY, maxX, maxY = getSidesDistances()
        row_begin, col_begin = toIndexes(dtm, minX, minY)
        row_end, col_end = toIndexes(dtm, maxX, maxY)



        valNoData = -9999
        gtiff_driver = agd.getdriver("GTiff")
        target_ds = agd.create( output_path, gtiff_driver, row_end-row_begin, col_end-col_begin, 1, agd.GDAL.GDT_Float32 )
        agd.setgeotransform!(target_ds, [ minX, resolution, 0.0, maxY, 0.0, -resolution ])
        agd.setproj!( target_ds, refsys )
     """ NON SO QUALE SIA IL COMANDO PER SETTARE I METADATI CON `ArchGDAL`
        target_ds.SetMetadata(
            Dict(
                "credits" => "Envifate - Francesco Geri, Oscar Cainelli, Paolo Zatelli, Gianluca Salogni, Marco Ciolli - DICAM Università degli Studi di Trento - Regione Veneto",
                "modulo" => "Analisi inquinamento luminoso",
                "descrizione" => "Simulazione di inquinamento luminoso da sorgente puntuale singola o multipla",
                "srs" => refsys,
                "data" => today()
            )
        )
     """
        band1 = agd.getband(target_ds, 1)
        agd.setnodatavalue!( band1, Float64(valNoData) )
        agd.fillraster!(band1, valNoData)
        band = agd.read(band1) 

        
        # Compute the light intensity in said area of the dtm
        for row in row_begin:row_end
            for col in col_begin:col_end
                if viewshed[row, col] == 1
                    x, y = toCoords(dtm, row, col)
                    dist = √( (y - y_source)^2 + (x - x_source)^2 )
                    new_intensity = intensity / dist

                    res = new_intensity > 0 ? new_intensity : 0
                    if nfeature == 1
                        band[row, col] = res
                    else
                        band[row, col] += res
                    end
                else
                    band[row, col] = 0
                end
            end
        end
    end


    # MANCA SCRIVERE IL RASTER




end




function run_light( dem, source, resolution::Integer, srource_height::Real=0.0, observer_height::Real=1.75, rarefraction::Real=0.14286, processing_memory::Integer=500, output_path::AbstractString=".\\light_intensity.tiff"  )

    if agd.geomdim(source) != 0
        throw(DomainError(source, "`source` must be a point"))
    end

    refsys = agd.getspatialref(source)

    if agd.importWKT(agd.getproj(dem)) != refsys
        throw(DomainError("The reference systems are not uniform. Aborting analysis."))
    end


    nfeature = 0
    features = agd.getgeom.( agd.getfeature(source) )
    for feature in features
        x_source = agd.getx(geom, 0)
        y_source = agd.gety(geom, 0)


        idxlevel = self.source.fields().indexFromName('level')
        intensity=feature.attributes()[idxlevel]

        nfeature += 1


        
        


        grass_area=str(x_min)+','+str(x_max)+','+str(y_min)+','+str(y_max)+' ['+str(self.areastudio.crs().authid())+']'
        grass_coord=str(x_source)+','+str(y_source)+' ['+str(self.source.crs().authid())+']'

        nameviewshed='viewshedanalysis'+str(nfeature)

        params = { '-b' : True, '-c' : False, '-e' : False, '-r' : False, 'GRASS_RASTER_FORMAT_META' : '', 'GRASS_RASTER_FORMAT_OPT' : '',
                   'GRASS_REGION_CELLSIZE_PARAMETER' : 0, 'GRASS_REGION_PARAMETER' :grass_area, 'coordinates' : grass_coord,
                   'input' : self.dem.dataProvider().dataSourceUri(), 'max_distance' : -1, 'memory' : self.memory, 'observer_elevation' : self.hsource, 'output' : nameviewshed,
                   'refraction_coeff' : self.rarefraction, 'target_elevation' : self.htarget }
        viewshed_proc = processing.run('grass7:r.viewshed', params)




        #aggiungo per controllo la viewshed alla toc
        iface.addRasterLayer(viewshed_proc['output'])

        viewshed=QgsProject.instance().mapLayersByName(nameviewshed)[0]


        for row in 1:rows
            for col in 1:cols
                x = col*pixel_size+x_min+(pixel_size/2)
                y = row*pixel_size+y_min+(pixel_size/2)

                punto_controllo = QgsPointXY(x,y)

                for pol in polygons
                    poly = pol.geometry()
                    if poly.contains(punto_controllo)

                        cfr_viewshed=viewshed.dataProvider().identify(QgsPointXY(x, y),QgsRaster.IdentifyFormatValue)

                        if cfr_viewshed.results()[1] == 1
                            deltax=x-x_source
                            deltay=y-y_source
                            dist=math.sqrt(math.pow(deltay,2)+math.pow(deltax,2))

                            #new_intensity=(1/math.pow(dist,2))*intensity
                            new_intensity=(1/dist)*intensity

                            if nfeature == 1
                                if new_intensity > 0
                                    outData[row,col]=new_intensity
                                else
                                    outData[row,col]=0
                                end
                            else
                                if new_intensity > 0
                                    outData[row,col]=outData[row,col]+new_intensity
                                else
                                    outData[row,col]=outData[row,col]
                                end
                            end
                        else
                            outData[row,col]=0
                        end
                    end
                end
            end
        end
    end



    valNoData = -9999

    gtiff_driver = agd.getdriver("GTiff")
    target_ds = agd.create( output_path, gtiff_driver, round(Int64, x_res), round(Int64, y_res), 1, agd.GDAL.GDT_Float32 )
    agd.setgeotransform!(target_ds, [ x_min, resolution, 0.0, y_max, 0.0, -resolution ])
    agd.setproj!( target_ds, refsys )
 """ NON SO QUALE SIA IL COMANDO PER SETTARE I METADATI CON `ArchGDAL`
    target_ds.SetMetadata(
        Dict(
            "credits" => "Envifate - Francesco Geri, Oscar Cainelli, Paolo Zatelli, Gianluca Salogni, Marco Ciolli - DICAM Università degli Studi di Trento - Regione Veneto",
            "modulo" => "Analisi inquinamento luminoso",
            "descrizione" => "Simulazione di inquinamento luminoso da sorgente puntuale singola o multipla",
            "srs" => refsys,
            "data" => today()
        )
    )
 """
    band1 = agd.getband(target_ds, 1)
    agd.setnodatavalue!( band1, Float64(valNoData) )
    agd.fillraster!(band1, valNoData)
    xsize = agd.width(band1)
    ysize = agd.height(band1)
    band = agd.read(band1)
    # outData = deepcopy(band)
    outData = band
    rows=ysize-1
    cols=xsize-1


        self.label_status.setText("Preparazione output")
        self.label_status.setStyleSheet('color : #e8b445;font-weight:bold')

        outData_raster=outData[::-1]
        band.WriteArray(outData_raster)


    band= None
    target_ds = None



    base_raster_name=os.path.basename(self.path_output)
    raster_name=os.path.splitext(base_raster_name)[0]
    self.iface.addRasterLayer(self.path_output, raster_name)


    layer=None
    for lyr in list(QgsProject.instance().mapLayers().values()):
        if lyr.name() == raster_name
            layer = lyr
        end
    end


    functions.applystyle(layer,'gr',0.5)


    tempoanalisi=time.time() - start_time
    tempostimato=time.strftime("%H:%M:%S", time.gmtime(tempoanalisi))
    messaggio="---------------------------------\n"
    messaggio+="Fine modellazione\n"
    messaggio+="\nTempo di analisi: "+tempostimato+"\n"
    messaggio+="---------------------------------\n\n"
    self.console.appendPlainText(messaggio)

    self.label_status.setText("In attesa di dati")
    self.label_status.setStyleSheet('color : green; font-weight:bold')
end