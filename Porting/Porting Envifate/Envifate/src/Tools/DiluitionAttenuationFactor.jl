module DiluitionAttenuationfactor
"""
"""



using Dates


import ArchGDAL as agd


include("..\\Library\\Functions.jl")



mutable struct DAF
    secondary_source_concentration::Float64
    x::Float64
    y::Float64
    α_x::Float64
    α_y::Float64
    decay_coeff::Float64
    darcy_velocity::Float64
    kd
    soil_density::Float64
    tera_e
    orthogonal_extension::Float64
    time

    acquifer_flow_direction::Float64
    algorithm::Symbol
    option::Symbol
  
    R
    DAF
    DAF_tot
  
    DAF(secondary_source_concentration,x,y,α_x,α_y,decay_coeff,darcy_velocity,kd,soil_density,tera_e,orthogonal_extension,time,acquifer_flow_direction,algorithm,option) = new(secondary_source_concentration,x,y,α_x,α_y,decay_coeff,darcy_velocity,kd,soil_density,tera_e,orthogonal_extension,time,acquifer_flow_direction,algorithm,option)
end

mutable struct Leach
    h
    tera_w
    tera_a
    kd
    effective_infiltration::Float64
    soil_density::Float64
    source_thickness::Float64
    aquifer_depth::Float64
    darcy_velocity::Float64
    mixed_zone_depth::Float64
    orthogonal_width::Float64
  
    kw
    koc
    ldf
    sam
    leaching_factor
  
    Leach(h,tera_w,tera_a,kd,effective_infiltration,soil_density,source_thickness,aquifer_depth,darcy_velocity,mixed_zone_depth,orthogonal_width) = new(h,tera_w,tera_a,kd,effective_infiltration,soil_density,source_thickness,aquifer_depth,darcy_velocity,mixed_zone_depth,orthogonal_width)
end



# ================================================ DAF functions ====================================================================================

#   LE FUNZIONI DI `DAF` USANO LA FUNZIONE erf DI PYTHON IN JULIA TALE FUNZIONE SI TROVA NEL PACCHETTO `SpecialFunctions.jl`

function calc_R!(c::DAF)
    c.R = 1 + ( c.kd * ( c.ro_s / c.tera_e ) )
    return c.R
end  


function calc_DAF_ispra!(c::DAF)
    ######################## modello di domenico ###########################
    # vedere appendice C pagina 2 del documento Criteri metodologici per l'applicazione dell'analisi assoluta di rischio ai siti contaminati
    # la formula originale prevede la produttoria delle 3 componenti x,y,z moltiplicata per 1/4
    # eliminando la terza componente dell'asse z è necessario moltplicare per 1/2 (quindi 0.5)
    # per verifica vedere Domenico P.A. e Schwartz F.W. (1998), Physical and Chemical Hydrogeology, John Wiley and Sons, New York.
    # da pagina 642 a pag 644
    if c.α_x == 0
      c.α_x = 0.1c.x
    end
    if c.α_y == 0
      c.α_y = c.α_x / 3
    end
  
    R = 1 + ( c.kd * ( c.ro_s / c.tera_e ) )
    daf1 = 0.50ℯ^( ( c.x / 2c.α_x ) * ( 1 - √( 1 + ( ( 4c.decay_coeff * c.α_x * R ) / c.v_e ) ) ) )
    #daf1 = exp( ( c.x / ( 2c.α_x ) ) )
    #daf2 = erf( c.s_w / ( 4√( c.α_y * c.x ) ) )
    daf21 = erf( ( c.y + 0.5c.s_w ) / ( 2√( c.α_y * c.x ) ) )
    daf22 = erf( ( c.y - 0.5c.s_w ) / ( 2√( c.α_y * c.x ) ) )
    #daf_prova = erf( ( c.y + 0.5c.s_w ) / ( 2√( c.α_y * c.x ) ) )
    daf3 = daf21 - daf22
    DAF_tot = daf1 * daf3
  
    return DAF_tot
end


function calc_DAF_ispra2!(c::DAF)
    if c.α_x == 0
      c.α_x = 0.1c.x
    end
    if c.α_y == 0
      c.α_y = c.α_x / 3
    end
    
    #daf1 = ( c.x / 2c.α_x ) * ( 1 - √( 1 + ( ( 4c.decay_coeff * c.α_x * c.R ) / c.v_e ) ) )
    daf1 = exp( c.x / ( 2c.α_x ) * 0 )
    #daf1e = exp(daf1)
    daf2 = erf( c.s_w / ( 4√( c.α_y * c.x ) ) )
    c.DAF = daf1 * daf2
  
    return c.DAF
end


function calc_DAF!(c::DAF)
    if c.α_x == 0
      c.α_x = 0.1( c.x / 100 )
    end
    if c.α_y == 0
      c.α_y = c.α_x / 3
    end
  
    dx = c.α_x * c.v_e
    dy = c.α_y * c.v_e
    daf_a = c.secondary_source_concentration / ( 4c.tera_e * π * c.T * √(dx * dy) )
    daf_b = ℯ^( -( ( (( c.x - (c.v_e * c.T) )^2) / (4dx * c.T) ) + ((c.y^2) / ( 4dy *c.T )) ) )   
    c.DAF = daf_a * daf_b
  
    return c.DAF
end


function calc_DAF_uni!(c::DAF)
    if c.α_x == 0
      c.α_x = 0.1c.x
    end
  
    dx = c.α_x * c.v_e
    daf_a = c.secondary_source_concentration / ( 2c.tera_e * √( 4dx * π * c.T ) ) 
    daf_b = ℯ^( -(( ((c.x - (c.v_e * c.T))^2) / ( 4dx * c.T ) )) )
    c.DAF = daf_a * daf_b
  
    return c.DAF
end


function calc_DAF_c!(c::DAF)
    #continuous
    if c.α_x == 0
      c.α_x =  0.1c.x
    end
    if c.α_y == 0
      c.α_y = c.α_x / 3
    end
  
    dx = c.α_x * c.v_e
    dy = c.α_x * c.v_e
    r = √( (c.x^2) + ( (c.y^2) * ( dx / dy) ) )
    daf_a = c.secondary_source_concentration / ( 4c.tera_e * √(π) * √( c.v_e * r ) * √(dy) )
    daf_b = ℯ^( ( ( c.x - r ) * c.v_e ) / ( 2 * dx ) ) 
    c.DAF = daf_a * daf_b
  
    return c.DAF
end

function calcDAF!(d::DAF)
    if d.x <= 0
        return 0.0
    else
        concentration = 0.0
        if d.algorithm == :fickian
            if d.option == :pulse
                concentration = calc_DAF!(d)
            else
                concentration = calc_DAF_c!(d)
            end
        else
            concentration = d.secondary_source_concentration * calc_DAF_ispra!(d)
        end
        return concentration
    end
end



# =============================================== Leaching functions ================================================================================

function calc_kw!( l::Leach )
    l.koc = l.kd
    l.kd = 0.01l.koc
    l.kw = l.soil_density / ( l.tera_w + ( l.kd * l.soil_density ) + ( l.h * l.tera_a ) )
    return l.kw
end


function calc_ldf!( l::Leach )
    darcy = l.darcy_velocity * 100.0 * 86400.0 * 365.0
    l.ldf = 1 + ( darcy * ( l.mixed_zone_depth / ( l.effective_infiltration * l.W ) ) )
    return l.ldf
end


function calc_sam!( l::Leach )
    l.sam = l.dz/l.lf
    return l.sam    
end
  

function calc_LF!( l::Leach )
    l.leaching_factor = ( l.kw * l.sam ) / l.ldf
    return l.leaching_factor
end



"""
    compute_result!( dtm::AbstractArray, r0::Integer, c0::Integer, ri::Integer, ci::Integer, daf::DAF )

Given the raster `dtm` and the indexes (`r0`, `c0`) of the source, modify the postion values of object `daf` and return the concentration at indexes (`ri`, `ci`)
"""
function compute_result!( dtm::AbstractArray, r0::Integer, c0::Integer, ri::Integer, ci::Integer, daf::DAF )
    daf.x, daf.y = Functions.compute_position(dtm, r0, c0, ri, ci, daf.acquifer_flow_direction)
    return calcDAF!(daf)
end



"""
    function leach( source, contaminants, concentrations, aquifer_depth, acquifer_flow_direction, mean_rainfall, texture, resolution::Integer, time::Integer=1,
                    orthogonal_extension::Real=10000.0, soil_density::Real=1.70, source_thickness::Real=1.0, darcy_velocity::Real=0.000025, mixed_zone_depth::Real=1.0,
                    decay_coeff::Real=0.0, algorithm::Symbol=:fickian, option::Symbol=:continuous, output_path::AbstractString=".\\output_model_daf.tiff" )

Run the simulation of leaching and dispersion of contaminants in an aquifier, returning a map of the possible worst case spreading of the contaminants

# Arguments
- `source`: source point of the contaminants.
- `contaminants`: type of substance.
- `concentrations`: concentration of the contaminants at the source.
- `aquifer_depth`: depth of the aquifier in meters.
- `acquifer_flow_direction`: angle of direction of the flow in degrees.
- `mean_rainfall`: average rainfall volume.
- `texture`: type of terrain at the source.
- `resolution::Integer`: dimension of a cell for the analysis.
- `time::Integer=1`: starting time.
- `orthogonal_extension::Real=10000.0`: X
- `soil_density::Real=1.70`: density of the terrain.
- `source_thickness::Real=1.0`: thickness of the terrain layer at the source.
- `darcy_velocity::Real=0.000025`: X
- `mixed_zone_depth::Real=1.0`: X
- `decay_coeff::Real=0.0`: X
- `algorithm::Symbol=:fickian`: type of algorithm to be used.
- `option::Symbol=:continuous`: second option to define the kind o algorithm to use.
- `output_path::AbstractString=".\\output_model_daf.tiff": output file path. 
"""
function leach( source, contaminants, concentrations, aquifer_depth, acquifer_flow_direction, mean_rainfall, texture, resolution::Integer, time::Integer=1,
                orthogonal_extension::Real=10000.0, soil_density::Real=1.70, source_thickness::Real=1.0, darcy_velocity::Real=0.000025, mixed_zone_depth::Real=1.0,
                decay_coeff::Real=0.0, algorithm::Symbol=:fickian, option::Symbol=:continuous, output_path::AbstractString=".\\output_model_daf.tiff" )

    if algorithm ∉ [:fickian, :domenico]
        throw(DomainError(algorithm, "`algorithm` must either be `:fickian` or `:domenico`"))
    end

    if option ∉ [:pulse, :continuous]
        throw(DomainError(option, "`option` must either be `:continuous` or `:pulse`"))
    end

    geom = agd.getgeom(collect(agd.getlayer(source, 0))[1])

    if agd.geomdim(geom) != 0
        throw(DomainError(source, "`source` must be a point"))
    end

    refsys = agd.getproj(dem)

    if agd.importWKT(refsys) != agd.getspatialref(geom)
        throw(DomainError("The reference systems are not uniform. Aborting analysis."))
    end

    effective_infiltration *= (mean_rainfall / 10.2)^2

    x_source = agd.getx(geom, 0)
    y_source = agd.gety(geom, 0)
    r_source, c_source = toIndexes(dtm, x_source, y_source)

    tera_a, tera_w, effective_infiltration, tera_e, grain = Functions.texture_extract( texture, ["tot_por", "c_water_avg", "effective_infiltration", "por_eff", "grain"], ".\\..\\library\\" )
    if all(isempty.([tera_a, tera_w, effective_infiltration, tera_e, grain]))
        throw(DomainError("Analysis error, check input parameters"))
    end

    gtiff_driver = agd.getdriver("GTiff")

    points = []
    values = []
    for i in 1:length(contaminants)
        path = output_path * "output_model_$(contaminants[i]).tiff"

        push!( points, [(r_source, c_source)] )
        push!( values, [concentrations[i]] )

        h, kd = Functions.substance_extract( contaminants[i], ["c_henry", "koc_kd"], ".\\..\\library\\" )
        if isempty(h) && isempty(kd)
            throw(DomainError("Analysis error, check input parameters"))
        end

        #                                       ief,                    ro,           dz,               lf,            ve,             dgw               sw
        element = Leach( h, tera_w, tera_a, kd, effective_infiltration, soil_density, source_thickness, aquifer_depth, darcy_velocity, mixed_zone_depth, orthogonal_extension ) 
        calc_kw!(element)
        calc_ldf!(element)
        calc_sam!(element)
        leaching_factor = calc_LF!(element)
        secondary_source_concentration = concentrations[i] * leaching_factor
        daf = DAF( secondary_source_concentration, x_source, y_source, 0, 0, decay_coeff, darcy_velocity, kd, soil_density, tera_e, orthogonal_extension, time, acquifer_flow_direction, algorithm, option )
        # Fill points with the indexes of each point that will compose the result raster and values with the concentrations in the respective points 
        expand!( points[i], values[i], dtm, r_source, c_source, daf )

        maxR = maximum( point -> point[1], points )
        minR = minimum( point -> point[1], points )
        maxC = maximum( point -> point[2], points )
        minC = minimum( point -> point[2], points )

 #= SENZA FUNZIONI
        rows = maxR - minR
        cols = maxC - minC
        minX, maxY = toCoords(dtm, minX, maxY)

        gtiff_driver = agd.getdriver("GTiff")
        target_ds = agd.create( path, gtiff_driver, rows, cols, 1, agd.GDAL.GDT_Float32 )
     # NON SONO CERTO CHE IL GEOTRASFORM VADA BENE
        agd.setgeotransform!( target_ds, [ minX, resolution, 0.0, maxY, 0.0, -resolution ] )
        agd.setproj!( target_ds, refsys )
        valNoData = -9999.0
        band1 = agd.getband( target_ds, 1 )
        agd.setnodatavalue!(band1, valNoData)
        agd.fillraster!(band1, valNoData)
        band = agd.read(band1)

        for (point, value) in zip(points[i], values[i])
            r, c = point - (minR, minC)
            band[r, c] = value
        end
    end
 =#

        geotransform = agd.getgeotransform(dem)
        geotransform[[1, 4]] .+= (minR - 1, maxC - 1) .* geotransform[[2, 6]]

        #   data = [ isnothing( findfirst(p -> p == (r, c), points) ) ? noData : values[findfirst(p -> p == (r, c), points)] for r in minR:maxR, c in minC:maxC ]
        noDataValue = -9999.f0
        data = fill(noDataValue, maxR-minR, maxC-minC)
        for r in minR:maxR, c in minC:maxC
            match = findfirst(p -> p == (r, c), points)
            if !isnothing(match)
                data[r-minR+1, c-minC+1] = values[match]
            end
        end
        Functions.writeRaster(data, gtiff_driver, geotransform, resolution, refsys, noDataValue, path, false)
    end


 """" AGGIUNTA DI UN LAYER AL RASTER FINALE
    band = nothing
    target_ds = nothing

    base_raster_name=os.path.basename(ef.path_output+str(sostanza))
    raster_name=os.path.splitext(base_raster_name)[0]
    ef.iface.addRasterLayer(ef.path_output, raster_name)

    contatore_sostanza=0
    ef.list_result=[]
 """


end



end # module





















#------------------------------------------------ TESTING------------------------------------------------------------------------------





import ArchGDAL as agd
using Plots
# Leggi il vettoriale
sat_file = split( @__DIR__, "\\Tirocinio\\")[1] * "\\Tirocinio\\Mappe\\sat\\sette_sorelle.shp"
sat = agd.read(sat_file)
sat_layer = agd.getlayer(sat, 0)
plot(sat_layer)
agd.imread(sat)















import ArchGDAL as agd
using Plots



# Leggi il vettoriale
cmn_file = split( @__DIR__, "\\Tirocinio\\")[1] * "\\Tirocinio\\Mappe\\c0104011_Comuni\\c0104011_Comuni.shp"
cmn = agd.read(cmn_file)

# Ottenere il layer; un layer raccoglie le features
layer = agd.getlayer(cmn, 0)

# Ottenere tutti i layers del vettoriale
layers_num = agd.nlayer(cmn)
layers = [ agd.getlayer(cmn, i) for i in 0:layers_num-1 ]

# Per fare indexing sulle features bisogna fare collect sul layer
features = collect(layer)
# Altrimenti si può semplicemente iterare sul layer
#   for feature in layer
#       println(feature)
#   end

# Ritorna il campo geometry della feature
geometry = agd.getgeom(features[1])

# Plottare le geometrie del vettoriale
 # Per i dati del vettoriale immagino che si debba selezionare il campo di
  # interesse e plottare quello
plot(geometry)
for feature in features[2:end]
    plot!( agd.getgeom(feature))
end
current()

# Ritorna il primo campo della feature
field0 = agd.getfield( features[1], 0 )

# Per trovare l'indice dato il nome del campo
fiel_index = agd.findfieldindex.( Ref(features[1]), [:id_stazion, :comune] )

# Come sopra ma con simboli (nomi dei campi invece che loro indice)
field0 = agd.getfield( features[1], :CODISTAT )

# Testare se il layer di un vettore ha un certo tipo di geometria
agd.getgeomtype(layer) == agd.wkbPolygon


agd.getspatialref(layer)


agd.create( "D:\\Z_Tirocinio_Dati\\source.shp", driver=agd.getdriver("ESRI Shapefile") ) do ds
    agd.createlayer( geom=GDAL.wkbPoint, spatialref=agd.importEPSG(4326) ) do layer
        agd.createfeature(layer) do feature
            agd.setgeom!( feature, agd.createpoint(lon, lat) )
        end
        agd.copy(layer, dataset=ds)
    end
end














import ArchGDAL as agd

dtm_file = split( @__DIR__ , "\\Porting\\")[1] * "\\Mappe\\DTM_wgs84.tiff"
dtm = agd.read(dtm_file)


# Accedere al dataset e farci indexing direttamente
dtm = agd.readraster(dtm_file)
dtm[4000, 6000]

# Per modificare è necessario leggere la band che si vuole modificare e poi leggerla di nuovo per ottenere una matrice modificabile
 # getband lo apre in read only
band1 = agd.getband(dtm, 1)
band = agd.read(band1)

# Per ottenere lo spatial reference system del dtm
agd.importWKT(agd.getproj(dtm))







import ArchGDAL as agd
import GDAL as gdl
dtm_file = split( @__DIR__ , "\\Porting\\")[1] * "\\Mappe\\DTM_32.tiff"
dtm = agd.read(dtm_file)


split( agd.gdalinfo( dtm ), "\n" , keepempty=false )
str = pointer("Center")
opt = Base.unsafe_convert(Ptr{Cstring}, str)
option = gdl.gdalinfooptionsnew(opt, C_NULL)
gdl.gdalinfo(dtm, option)

a = String["Center"]
str = pointer(a)








import ArchGDAL as agd
import GeoArrays as ga



vntst_file = "C:\\Users\\Lenovo\\Documents\\GitHub\\Tirocinio\\Mappe\\Stazioni_Veneto\\Stazioni_Veneto.shp"
vntst = agd.read(vntst_file)
layer = agd.getlayer(vntst, 0)
features = collect(layer)










# Raster
dtm_file = split( @__DIR__ , "\\Porting\\")[1] * "\\Mappe\\DTM_32.tiff"
dtm = agd.readraster(dtm_file)
gdtm = ga.read(dtm_file)

x, y = ga.coords(gdtm, [4000, 6000])
x1, y1 = toCoords(dtm, 4000, 6000)




dtm_band = agd.getband(dtm, 1)
agd.toPROJ4(agd.importWKT(agd.getproj(dtm)))

agd.getgeotransform(dtm)

band = agd.getband(dtm, 1)
band1 = agd.read(band)


# Vector


cmn_file = split( @__DIR__ , "\\Porting\\")[1] * "\\Mappe\\c0104011_Comuni\\c0104011_Comuni.shp"
cmn = agd.read(cmn_file)

cmn_layer = agd.getlayer(cmn, 0)

agd.envelope(cmn_layer)








# crs del dtm
dtm_file = split( @__DIR__ , "\\Porting\\")[1] * "\\Mappe\\DTM_32.tiff"
dtm = agd.read(dtm_file)
crs_dtm = agd.importWKT(agd.getproj(dtm))


band = agd.getband(dtm, 1)
agd.getspatialref(band)


# crs delle stazioni
sts_file = split( @__DIR__ , "\\Porting\\")[1] * "\\Mappe\\Stazioni_Veneto\\Stazioni_Veneto.shp"
sts = agd.read(sts_file)
sts_layers = collect(agd.getlayer(sts, 0))
crs_src = agd.getspatialref(agd.getlayer(sts, 0))

# crs della sorgente
point = agd.getgeom(sts_layers[1])
agd.geomdim(point)
crs_point = agd.getspatialref(point)

















import ArchGDAL as agd
using Rasters
using Shapefile

dtm_file = split( @__DIR__ , "\\Porting\\")[1] * "\\Mappe\\DTM_wgs84.tiff"
#   con_file = "C:\\Users\\DAVIDE-FAVARO\\Desktop\\Connectivity Data\\connectivity.tiff"
sat_file = split( @__DIR__ , "\\Porting\\")[1] * "\\Mappe\\sat\\sette_sorelle.shp"

dtmr = Raster(dtm_file)
dtma = agd.read(dtm_file)

conr = Raster(con_file)
cona = agd.read(con_file)

sats = Shapefile.Table(sat_file)
sata = agd.read(sat_file)














#---------------------------------------------------------------------------------------------------------------------------------------