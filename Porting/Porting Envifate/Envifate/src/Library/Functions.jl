module Functions

export substance_extract, texture_extract, air_extract, cn_extract, cn_list_extract,
       array2raster!, writeRaster!, applystyle

import ArchGDAL as agdal
import DBInterface as dbi
import SQLite as sql


function substance_extract( substance_id, fields, dbloc = "" )
    #estrazione valori sostanze
    db = sql.DB(dbloc*"substance.db")
    sql_fields = join( fields, "," )
    query_substance = sql.Stmt( db, "SELECT ? FROM substance WHERE id = ?AAA" )
    sql.bind!( query_substance, [ sql_fields, substance_id ] )
    results = dbi.execute( query_substance )
    res_fields = [ x for x in results ]
    return res_fields
end


function texture_extract( texture_name, fields, dbloc = "" )
    #estrazione valori sostanze
    db = sql.DB(dbloc*"substance.db")
    sql_fields = join( fields, "," )
    query_texture = sql.Stmt( db, "SELECT ? FROM texture WHERE nome LIKE ?" )
    sql.bind!( query_texture, [ sql_fields, texture_name ] )
    results = dbi.execute( query_texture ) 
    res_fields = [ x  for x in results ]
    return res_fields
end

function air_extract( stability_class, outdoor, dbloc::AbstractString=*( @__DIR__, "\\") )
    db = sql.DB(dbloc*"substance.db")
    query_texture = sql.Stmt( db, "SELECT sigmay1, sigmay2, sigmayexp, sigmaz1, sigmaz2, sigmazexp FROM air_stability WHERE class LIKE ?NNN AND outdoor LIKE ?NNN" )
    sql.bind!( query_texture, [ stability_class, outdoor ] )
    results = dbi.execute( query_texture )
    res_fields = [ x for x in results ]
    return res_fields
end


function cn_extract( cnl, soil, dbloc::AbstractString=*( @__DIR__, "\\") )
    db = sql.DB(dbloc*"substance.db")
    classecn = "cn_"*String(cnl)
    query_cn = sql.Stmt( db, "SELECT ? FROM cn WHERE id = ?AAA" )
    sql.bind!( query_cn, [ classecn, soil ] )
    results = dbi.execute(query_cn)
    res_fields = [ x for x in results ]
	return res_fields
end


function cn_list_extract( dbloc::AbstractString=*( @__DIR__, "\\") )
	db = sql.DB(dbloc*"substance.db")
    query_cn = sql.Stmt( db, "SELECT * FROM cn" )
    results = dbi.execute(query_cn)
#    listaclc = Dict()
#    for row in results
#        lista_soil = [ x for x in row ]
#        listaclc[ row[5] ] = lista_soil
#    end
    listaclc = Dict( row[5] => [ x for x  in row ] for row in results )
    return listaclc
end


function array2raster!( newRasterfn, xmin, ymin, pixelWidth, pixelHeight, xsize, ysize, array )
	# vedi https://pcjericks.github.io/py-gdalogr-cookbook/raster_layers.html#create-raster-from-array
    # cols = array.shape[1]
    # rows = array.shape[0]
    cols = xsize
    rows = ysize
    originX = xmin
    originY = ymin

    driver = agdal.getdriver("GTiff")
    outRaster = agdal.create( newRasterfn, driver, cols, rows, 1, ArchGDAL.GDT_Byte )
    agdal.setgeotransform!( outRaster, [ originX, pixelWidth, 0, originY, 0, pixelHeight ] )
    outband = agdal.getband(1)
    agdal.write!( outband, array )
    #outRasterSRS = agdal.importEPSG(4326)
    #agdal.setproj!( outRaster,  agdal.toWKT(outRasterSRS) )
"""
    outband.FlushCache()
"""
end


function writeRaster!( newRasterfn, xmin, ymin, pixelWidth, pixelHeight, xsize, ysize, array )
    reversed_arr = reverse(array) # reverse array so the tif looks like the array
    array2raster!( newRasterfn, xmin, ymin, pixelWidth, pixelHeight, xsize, ysize, reversed_arr ) # convert array to raster
end



function applystyle( layer, colore, opacity ) end

end # module