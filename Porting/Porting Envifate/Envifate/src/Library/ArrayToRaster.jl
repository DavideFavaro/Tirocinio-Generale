module ArrayToRaster

export array2raster!, writeRaster! #functions

#= imports
    import gdal, ogr, os, osr
    import numpy as np
=#

import ArchGDAL as agd


function array2raster!( newRasterfn, rasterOrigin, pixelWidth, pixelHeight, array )
    cols = size(array)[2]
    rows = size(array)[1]
    originX = rasterOrigin[1]
    originY = rasterOrigin[2]

    driver = agd.getdriver("GTiff")
    outRaster = agd.create( newRasterfn, driver, cols, rows, 1, ArchGDAL.GDT_Byte )
    agd.setgeotransform!( outRaster, [originX, pixelWidth, 0, originY, 0, pixelHeight] )
    outband = agd.getband( outRaster, 1 )
    agd.write!( outband, array )
    #outRasterSRS = agd.importEPSG(4326)
    #agd.setproj!( outRaster, agd.toWKT(outRasterSRS) ) 
"""
    outband.FlushCache()
"""
end


function writeRaster!( newRasterfn, rasterOrigin, pixelWidth, pixelHeight, array )
    reversed_arr = reverse(array) # reverse array so the tif looks like the array
    array2raster!( newRasterfn, rasterOrigin, pixelWidth, pixelHeight, reversed_arr ) # convert array to raster
end



if abspath(PROGRAM_FILE) == @__FILE__
    rasterOrigin = ( -123.25745, 45.43013 )
    pixelWidth = 10
    pixelHeight = 10
    newRasterfn = "test.tif"
    array = [ 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5;
              5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5;
              5 0 0 0 0 5 0 0 0 0 5 0 0 0 5 0 5 5 5;
              5 0 5 5 5 5 5 0 5 0 5 0 5 0 5 0 5 5 5;
              5 0 5 0 0 5 5 0 5 0 5 0 0 0 5 0 5 5 5;
              5 0 5 5 0 5 5 0 5 0 5 0 5 0 5 0 5 5 5;
              5 0 0 0 0 5 0 0 0 0 5 0 5 0 5 0 0 0 5;
              5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5;
              5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5;
              5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 ]

    writeRaster!( newRasterfn, rasterOrigin, pixelWidth, pixelHeight, array )
end

end # module