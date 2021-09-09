#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import print_function
from builtins import str
from qgis.PyQt import QtCore, QtGui
from PyQt5.QtGui import QColor
from qgis.core import *
import numpy as np
try:
  import sqlite3
except:
  # fix_print_with_import
  print("librerie per la connessione al database sqlite non trovate")

import os

from osgeo import gdal,ogr

def substance_extract(id_s,fields,dbloc=None):
	#estrazione valori sostanze
	if dbloc==None:
		dbloc=''
	conn = sqlite3.connect(dbloc+"substance.db")
	cursor=conn.cursor()

	sql_fields=""
	i=0
	for x in fields:
	  i+=1
	  sql_fields+=x
	  if i<len(fields):
	    sql_fields+=","
	query_substance="select "+sql_fields+" from substance where id="+str(id_s)

	cursor.execute(query_substance)

	sql_fetch=cursor.fetchall()

	res_fields=[]

	for row in sql_fetch:
	  for x in row:
	    res_fields.append(x)

	conn.close()

	return res_fields




def texture_extract(texture,fields,dbloc=None):
	#estrazione valori sostanze
	if dbloc==None:
		dbloc=''
	conn = sqlite3.connect(dbloc+"substance.db")
	cursor=conn.cursor()
	sql_fields=""
	i=0
	for x in fields:
	  i+=1
	  sql_fields+=x
	  if i<len(fields):
	    sql_fields+=","

	query_texture="select "+sql_fields+" from texture where nome like lower('"+str(texture)+"')"

	cursor.execute(query_texture)

	sql_fetch=cursor.fetchall()

	res_fields=[]

	for row in sql_fetch:
	  for x in row:
	    res_fields.append(x)

	conn.close()

	return res_fields



def air_extract(c_stability,outdoor,dbloc=None):
	if dbloc==None:
		dbloc=os.path.dirname(os.path.abspath(__file__))+"/"
	#print os.path.dirname(os.path.abspath(__file__))
	conn = sqlite3.connect(dbloc+"substance.db")
	cursor=conn.cursor()
	query_texture="select sigmay1,sigmay2,sigmayexp,sigmaz1,sigmaz2,sigmazexp from air_stability where class like lower('"+str(c_stability)+"') and outdoor like lower('"+str(outdoor)+"')"

	cursor.execute(query_texture)

	sql_fetch=cursor.fetchall()

	res_fields=[]

	for row in sql_fetch:
	  for x in row:
	    res_fields.append(x)

	conn.close()

	return res_fields



def cn_extract(cnl,soil,dbloc=None):
	if dbloc==None:
		dbloc=os.path.dirname(os.path.abspath(__file__))+"/"
	#print os.path.dirname(os.path.abspath(__file__))
	conn = sqlite3.connect(dbloc+"substance.db")
	cursor=conn.cursor()
	
	classecn="cn_"+str(cnl)
	query_cn="select "+classecn+" from cn where id = "+str(soil)
	#import pdb; pdb.set_trace()
	cursor.execute(query_cn)

	sql_fetch=cursor.fetchall()

	res_fields=[]

	for row in sql_fetch:
	  for x in row:
	    res_fields.append(x)

	conn.close()

	return res_fields


def cn_list_extract(dbloc=None):
	if dbloc==None:
		dbloc=os.path.dirname(os.path.abspath(__file__))+"/"
	#print os.path.dirname(os.path.abspath(__file__))
	conn = sqlite3.connect(dbloc+"substance.db")
	cursor=conn.cursor()

	listaclc={}

	query_cn="select * from cn"
	#import pdb; pdb.set_trace()
	cursor.execute(query_cn)

	sql_fetch=cursor.fetchall()


	for row in sql_fetch:

		lista_soil=[]
		for x in row:
			lista_soil.append(x)
		
		listaclc[row[5]]=lista_soil
	conn.close()

	return listaclc


def array2raster(newRasterfn,xmin,ymin,pixelWidth,pixelHeight,xsize,ysize,array):
	# vedi https://pcjericks.github.io/py-gdalogr-cookbook/raster_layers.html#create-raster-from-array
    # cols = array.shape[1]
    # rows = array.shape[0]
    cols = xsize
    rows = ysize
    originX = xmin
    originY = ymin

    driver = gdal.GetDriverByName('GTiff')
    outRaster = driver.Create(newRasterfn, cols, rows, 1, gdal.GDT_Byte)
    outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
    outband = outRaster.GetRasterBand(1)
    outband.WriteArray(array)
    #outRasterSRS = osr.SpatialReference()
    #outRasterSRS.ImportFromEPSG(4326)
    #outRaster.SetProjection(outRasterSRS.ExportToWkt())
    outband.FlushCache()


def writeraster(newRasterfn,xmin,ymin,pixelWidth,pixelHeight,xsize,ysize,array):
    reversed_arr = array[::-1] # reverse array so the tif looks like the array
    array2raster(newRasterfn,xmin,ymin,pixelWidth,pixelHeight,xsize,ysize,reversed_arr) # convert array to raster





def applystyle(layer,colore,opacity):
	provider = layer.dataProvider()
	ext = layer.extent()
	stats = provider.bandStatistics(1,QgsRasterBandStats.All,ext,0)
	bandmin=stats.minimumValue
	bandmax=stats.maximumValue
	interval_lst=np.linspace(bandmin,bandmax,10)

	lst_color_gr=['#3be607','#92db2d','#e9cf52','#feb751','#fe9b43','#fb7b35','#f55629','#eb3420','#d41a23','#bd0026']
	lst_color_viridis=['#440154','#472878','#3e4a89','#31688e','#25838e','#1e9e89','#35b779','#6cce59','#b5de2c','#fde725']
	lst_color_magma=['#000004','#180f3e','#440f76','#721f81','#9e2f7f','#cd3f71','#f1605d','#fd9567','#fec98d','#fcfdbf']

	dict_color={'gr':lst_color_gr,'viridis':lst_color_viridis,'magma':lst_color_magma}

	fcn = QgsColorRampShader()
	fcn.setColorRampType(QgsColorRampShader.Interpolated)
	lst=[]
	for x in range(0,10):
		lst.append(QgsColorRampShader.ColorRampItem(interval_lst[x], QColor(dict_color[colore][x])))
	fcn.setColorRampItemList(lst)
	shader = QgsRasterShader()
	shader.setRasterShaderFunction(fcn)
	renderer = QgsSingleBandPseudoColorRenderer(layer.dataProvider(), 1,shader)
	layer.setRenderer(renderer)

	layer.renderer().setOpacity(opacity)


