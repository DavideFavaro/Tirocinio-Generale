module DiluitionAttenuationfactor

# -*- coding: utf-8 -*-
#=
/***************************************************************************
 Envifate
                                 A QGIS plugin
 Envifate: Open source tool for environmental risk analysis
                              -------------------
        begin                : 2016-07-15
        git sha              : $Format:%H$
        copyright            : (C) 2016 by Francesco Geri
        email                : francescogeri@tim.com
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
 # AGGIUNTE ALL'ORIGINALE PER SEMPLIFICARE I CALCOLI
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

""" LE FUNZIONI DI `DAF` USANO LA FUNZIONE erf DI PYTHON IN JULIA TALE FUNZIONE SI TROVA NEL PACCHETTO `SpecialFunctions.jl`"""

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



#  FUNZIONI CHE PROBABILMENTE POSSONO ESSERE SPOSTATE IN UN FILE A PARTE (IDEALMENTE Functions)
"""
Recursively compute the concentration of each point and add the value and its indexes to positions
"""
function expand!( positions::AbstractVector, results::AbstractVector, dtm, indx_x::Integer, indx_y::Integer, daf::DAF )
    if (indx_x, indx_y) in positions
        xs = [ indx_x+1, indx_x, indx_x-1, indx_x ]
        ys = [ indx_y, indx_y+1, indx_y, indx_y-1 ]
        expand!.( Ref(positions), Ref(results), Ref(dtm), xs, ys, daf )
        return nothing
    else
        Δx, Δy = toCoords(dtm, positions[1][1], positions[1][2]) - toCoords(dtm, indx_x, indx_y)
        dir = deg2rad(daf.acquifer_flow_direction)
        cosdir = cos(dir)
        sindir = sin(dir)
        x = Δx * cosdir - Δy * sindir
        y = Δy * sindir + Δy * cosdir
        daf.x = x
        daf.y = y
        concentration = calcDAF!(daf)
        if round(concentration, digits=5) > 0
            push!( positions, (ind_x, ind_y) )
            push!( results, concentration )
            xs = [ indx_x+1, indx_x, indx_x-1, indx_x ]
            ys = [ indx_y, indx_y+1, indx_y, indx_y-1 ]
            expand!.( Ref(positions), Ref(results), Ref(dtm), xs, ys, daf )
        end
        return nothing
    end
end



function leach( source, contaminants, concentrations, aquifer_depth, acquifer_flow_direction, mean_rainfall, texture, resolution::Integer, time::Integer=1,
                orthogonal_extension::Real=10000.0, soil_density::Real=1.70, source_thickness::Real=1.0, darcy_velocity::Real=0.000025, mixed_zone_depth::Real=1.0,
                decay_coeff::Real=0.0, algorithm::Symbol=:fickian, option::Symbol=:continuous, output_path::AbstractString="" )

    if algorithm ∉ [:fickian, :domenico]
        throw(DomainError(algorithm, "`algorithm` must either be `:fickian` or `:domenico`"))
    end

    if option ∉ [:pulse, :continuous]
        throw(DomainError(option, "`option` must either be `:continuous` or `:pulse`"))
    end

    if agd.geomdim(source) != 0
        throw(DomainError(source, "`source` must be a point"))
    end

    refsys = agd.getspatialref(source)
    effective_infiltration *= (mean_rainfall / 10.2)^2

    feature = collect(agd.getfeature(source))
    geom = agd.getgeom(feature[1])
    x_source = agd.getx(geom, 0)
    y_source = agd.gety(geom, 0)
    r_source, c_source = toIndexes(dtm, x_source, y_source)

    tera_a, tera_w, effective_infiltration, tera_e, grain = Functions.texture_extract( texture, ["tot_por", "c_water_avg", "effective_infiltration", "por_eff", "grain"], ".\\..\\library\\" )
    if all(isempty.([tera_a, tera_w, effective_infiltration, tera_e, grain]))
        throw(DomainError("Analysis error, check input parameters"))
    end

    points = []
    values = []
    for i in 1:length(contaminants)
        path = output_path * "output_model_$(contaminants[i]).tiff"

        push!( points, [ (r_source, c_source) ] )
        push!( values, [ concentrations[i] ] )

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

        rows = maxR - minR
        cols = maxC - minC
        minX, maxY = toCoords(dtm, minX, maxY)

        gtiff_driver = agd.getdriver("GTiff")
        target_ds = agd.create( path, gtiff_driver, rows, cols, 1, agd.GDAL.GDT_Float32 )
     # NON SONO CERTO CHE IL GEOTRASFORM VADA BENE
        agd.setgeotransform!( target_ds, [ minX, resolution, 0.0, maxY, 0.0, -resolution ] )
        agd.setproj!( target_ds, refsys )
     """ NON SO QUALE SIA IL COMANDO PER SETTARE I METADATI CON `ArchGDAL`
        target_ds.SetMetadata(
            Dict(
                "credits" => "Envifate - Francesco Geri, Oscar Cainelli, Paolo Zatelli, Gianluca Salogni, Marco Ciolli - DICAM Università degli Studi di Trento - Regione Veneto",
                "modulo" => "Dispersione in falda",
                "descrizione" => "Simulazione di dispersione inquinante in falda",
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

        for (point, value) in zip(points[i], values[i])
            r, c = point - (minR, minC)
            band[r, c] = value
        end
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
cmn_file = "C:\\Users\\Lenovo\\Documents\\GitHub\\Tirocinio\\Mappe\\c0104011_Comuni\\c0104011_Comuni.shp"
cmn = agd.read(cmn_file)

# Ottenere il layer; un layer raccoglie le features
layer = agd.getlayer(cmn, 0)

# Ottenere tutti i layers del vettoriale
layers_num = agd.nlayer(cmn)
layers = [ agd.getlayer(cmn, i) for i in 0:layers_num-1 ]

# Per fare indexind sulle features bisogna fare collect sul layer
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




import ArchGDAL as agd

dtm_file = split( @__DIR__ , "\\Porting\\")[1] * "\\Mappe\\DTM_32.tiff"

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




#---------------------------------------------------------------------------------------------------------------------------------------