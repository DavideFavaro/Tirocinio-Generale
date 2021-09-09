import gdal, ogr, os, osr
import numpy as np


def array2raster(newRasterfn,rasterOrigin,pixelWidth,pixelHeight,array):

    cols = array.shape[1]
    rows = array.shape[0]
    originX = rasterOrigin[0]
    originY = rasterOrigin[1]

    driver = gdal.GetDriverByName('GTiff')
    outRaster = driver.Create(newRasterfn, cols, rows, 1, gdal.GDT_Byte)
    outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
    outband = outRaster.GetRasterBand(1)
    outband.WriteArray(array)
    #outRasterSRS = osr.SpatialReference()
    #outRasterSRS.ImportFromEPSG(4326)
    #outRaster.SetProjection(outRasterSRS.ExportToWkt())
    outband.FlushCache()


def writeraster(newRasterfn,rasterOrigin,pixelWidth,pixelHeight,array):
    reversed_arr = array[::-1] # reverse array so the tif looks like the array
    array2raster(newRasterfn,rasterOrigin,pixelWidth,pixelHeight,reversed_arr) # convert array to raster


if __name__ == "__main__":
    rasterOrigin = (-123.25745,45.43013)
    pixelWidth = 10
    pixelHeight = 10
    newRasterfn = 'test.tif'
    array = np.array([[ 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
                      [ 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
                      [ 5, 0, 0, 0, 0, 5, 0, 0, 0, 0, 5, 0, 0, 0, 5, 0, 5, 5, 5],
                      [ 5, 0, 5, 5, 5, 5, 5, 0, 5, 0, 5, 0, 5, 0, 5, 0, 5, 5, 5],
                      [ 5, 0, 5, 0, 0, 5, 5, 0, 5, 0, 5, 0, 0, 0, 5, 0, 5, 5, 5],
                      [ 5, 0, 5, 5, 0, 5, 5, 0, 5, 0, 5, 0, 5, 0, 5, 0, 5, 5, 5],
                      [ 5, 0, 0, 0, 0, 5, 0, 0, 0, 0, 5, 0, 5, 0, 5, 0, 0, 0, 5],
                      [ 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
                      [ 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
                      [ 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5]])


    writeraster(newRasterfn,rasterOrigin,pixelWidth,pixelHeight,array)