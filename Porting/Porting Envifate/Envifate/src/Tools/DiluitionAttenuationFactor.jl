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

#= imports
    from __future__ import print_function
    from builtins import str
    from builtins import range
    from qgis.PyQt import QtCore, QtGui

    from PyQt5.QtCore import QSettings, QTranslator, QCoreApplication, Qt, QObject, pyqtSignal, pyqtRemoveInputHook
    from PyQt5.QtGui import QIcon, QStandardItem,QStandardItemModel
    from PyQt5.QtWidgets import QAction, QDialog, QFormLayout, QMenu, QComboBox, QTableWidgetItem, QHBoxLayout, QLineEdit, QPushButton, QWidget, QSpinBox, QTableWidget, QTableWidgetItem, QMessageBox, QFileDialog,QListView

    import sys

    # Initialize Qt resources from file resources.py
    #import resources


    # Import the code for the dialog

    #from open_risk_dialog import OpenRiskDialog
    import os.path
    try:
    import sqlite3
    except:
    # fix_print_with_import
    print("librerie per la connessione al database sqlite non trovate")

    import qgis
    from qgis.core import *
    from qgis.gui import *
    from qgis.utils import iface


    import numpy as np
    import math
    import struct
    import time, datetime
    import platform


    from osgeo import gdal,ogr,osr

    from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar

    import matplotlib.pyplot as plt


    import pdb


    from envifate_dialog import EnviDialog

    from configuration_dialog import ConfigurationDialog

    sys.path.append( os.path.dirname(__file__)+"/../library" )

    import functions, leaching, daf, do_setting
=#

import ArchGDAL as agd
using Dates

include("../Library/Daf.jl")
include("../Library/Leaching.jl")


Base.:-( x::Tuple{Number, Number}, y::Tuple{Number, Number} ) = ( x[1] - y[1], x[2] - y[2] )
Base.:-( x::Vector{T}, y::Tuple{T, T} ) where {T <: Number} = length(x) == length(y) ? [ e1 - e2 for (e1, e2) in zip(x, y) ] : throw(ArgumentError("`x` and `y` must have the same size"))
Base.:-( x::Tuple{T, T}, y::Vector{T} ) where {T <: Number} = length(x) == length(y) ? Tuple( e1 - e2 for (e1, e2) in zip(x, y) ) : throw(ArgumentError("`x` and `y` must have the same size"))
Base.:+( x::Tuple{Number, Number}, y::Tuple{Number, Number} ) = ( x[1] + y[1], x[2] + y[2] )
Base.:*( x::Tuple{Number, Number}, y::Number ) = ( x[1] * y, x[2] * y )
Base.:*( x::Number, y::Tuple{Number, Number} ) = y * x
Base.:*( x::Tuple{Number, Number}, y::Tuple{Number, Number} ) = ( x[1] * y[1], y[1] * y[2] )
Base.:/( x::Tuple{Number, Number}, y::Number ) = ( x[1] / y, x[2] / y )
Base.:/( x::Number, y::Tuple{Number, Number} ) = y / x
Base.:/( x::Tuple{Number, Number}, y::Tuple{Number, Number} ) = ( x[1] / y[1], x[2] / y[2] )
Base.:^( x::Tuple{Number, Number}, y::Number ) = ( x[1]^y, x[2]^y )

srid = [ 3003, 3004, 32632, 32633, 3857, 4326 ]

classiwind = Dict(
    "N" => 0,
    "NE" => 45,
    "E" => 90,
    "SE" => 135,
    "S" => 180,
    "SW" => 225,
    "W" => 270,
    "NW" => 315
)

classiwind2 = Dict(
    0 => 180,
    45 => 225,
    90 => 270,
    135 => 315,
    180 => 0,
    225 => 45,
    270 => 90,
    315 => 135
)

mutable struct DAF
    c0
    x
    y
    α_x::Float64
    α_y::Float64
    λ1
    darcy::Float64
    kd
    soild::Float64
    tera_e
    s_w
    T
  
    R
    DAF
    DAF_tot
  
    ClassDaf(c0,x,y,α_x,α_y,λ1,darcy,kd,soild,tera_e,s_w,T) = new(c0,x,y,α_x,α_y,λ1,darcy,kd,soild,tera_e,s_w,T)
end

mutable struct Leach
    h
    tera_w
    tera_a
    kd
    ief::Float64
    soild::Float64
    sourcethick::Float64
    aquifer_depth::Float64
    darcy::Float64
    Δgw::Float64
    W::Float64
  
    kw
    koc
    ldf
    sam
    LF
  
    Leach(h,tera_w,tera_a,kd,ief,ro_s,dz,lf,v_e,Δgw,W) = new(h,tera_w,tera_a,kd,ief,ro_s,dz,lf,v_e,Δgw,W)
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
    daf1 = 0.50ℯ^( ( c.x / 2c.α_x ) * ( 1 - √( 1 + ( ( 4c.λ1 * c.α_x * R ) / c.v_e ) ) ) )
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
    
    #daf1 = ( c.x / 2c.α_x ) * ( 1 - √( 1 + ( ( 4c.λ1 * c.α_x * c.R ) / c.v_e ) ) )
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
    daf_a = c.c0 / ( 4c.tera_e * π * c.T * √(dx * dy) )
    daf_b = ℯ^( -( ( (( c.x - (c.v_e * c.T) )^2) / (4dx * c.T) ) + ((c.y^2) / ( 4dy *c.T )) ) )   
    c.DAF = daf_a * daf_b
  
    return c.DAF
end
  
function calc_DAF_uni!(c::DAF)
    if c.α_x == 0
      c.α_x = 0.1c.x
    end
  
    dx = c.α_x * c.v_e
    daf_a = c.c0 / ( 2c.tera_e * √( 4dx * π * c.T ) ) 
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
    daf_a = c.c0 / ( 4c.tera_e * √(π) * √( c.v_e * r ) * √(dy) )
    daf_b = ℯ^( ( ( c.x - r ) * c.v_e ) / ( 2 * dx ) ) 
    c.DAF = daf_a * daf_b
  
    return c.DAF
end


# =============================================== Leaching functions ================================================================================

function calc_kw!( l::Leach )
    l.koc = l.kd
    l.kd = 0.01l.koc
    l.kw = l.soild / ( l.tera_w + ( l.kd * l.soild ) + ( l.h * l.tera_a ) )
    return l.kw
end

function calc_ldf!( l::Leach )
    v_darcy = l.darcy * 100.0 * 86400.0 * 365.0
    l.ldf = 1 + ( v_darcy * ( l.Δgw / ( l.ief * l.W ) ) )
    return l.ldf
end

function calc_sam!( l::Leach )
    l.sam = l.dz/l.lf
    return l.sam    
end
  
function calc_LF!( l::Leach )
    l.LF = ( l.kw * l.sam ) / l.ldf
    return l.LF
end


#------------------------------------------------ TESTING------------------------------------------------------------------------------

import ArchGDAL as agd

dtm_file = split( @__DIR__ , "\\Porting\\")[1] * "\\Mappe\\DTM_32.tiff"
dtm = agd.read(dtm_file)

agd.getgeotransform(dtm)

band = agd.getband(dtm, 1)
band1 = agd.read(band)

#---------------------------------------------------------------------------------------------------------------------------------------


function leach( source, area, contaminants, conc, aquifer_depth, dirwind, dirfalda, p,  texture, resolution::Int64, time::Int64=1, sw::Real=10000.0, soild::Real=1.70,
                sourcethick::Int64=1, darcy::Real=0.000025, Δgw::Int64=1, λ1::Real=0.0, algorithm::Symbol=:fickian, option::Symbol=:continuous, output_path::AbstractString="" )
    if algorithm ∉ [:fickian, :domenico]
        throw(DomainError(algorithm, "`algorithm` must either be `:fickian` or `:domenico`"))
    end
    if option ∉ [:pulse, :continuous]
        throw(DomainError(option, "`option` must either be `:continuous` or `:pulse`"))
    end

    x_w = azimut = dirfalda - 90

    # attenzione in realtà la massa è in grammi non in kg: correggere anche sopra alla linea 356
    #ef.C=ef.C_kg*1000


 """
    # wkbType: 1:point, 6:multipolygon, 2: Linestring

    if ef.source.wkbType()!=1:
        QMessageBox.warning(ef,"Warning", "La sorgente deve avere geometria puntuale" )
        return

    ef.areastudio=ef.listalayers[ef.text_area]


    if ef.areastudio.wkbType()!=6:
        QMessageBox.warning(ef,"Warning", u"L'area di studio deve avere geometria poligonale" )
        return
 """

    if agd.getspatialref(area) != agd.getspatialref(source)
        throw(DomainError("Warning", "Errore: i sistemi di riferimento non sono uniformi. Impossibile continuare con l'analisi." ))
    end
    refsys = agd.getspatialref(source)

    for (substance, concentration) in zip(contaminants, conc)
        if output_path == ""
            output_path = ".\\output_model_$sostanza.tiff"
        end

        res_fields = Functions.substance_extract( substance, ["c_henry", "koc_kd"], *( @__DIR__, "/../library/" ) )
        h, kd = res_fields
     # CHE `texture` STA USANDO? QUELLO PASSATO COME PARAMETRO? SE SI IL CALCOLO SI POTREBBE PORTARE FUORI DAL `for`
        res_texture = Functions.texture_extract( texture, ["tot_por", "c_water_avg", "ief", "por_eff", "grain"], *( @__DIR__, "/../library/" ) )
        tera_a, tera_w, ief, tera_e, grain = res_texture
     # DA DOVE SALTA FUORI `pioggia`
        ief *= (pioggia/10.2)^2

     """ PRINT DI COSE
        messaggio="Inizio elaborazione Acquifero Envifate\n"
        messaggio+="---------------------------\n\n"
        messaggio+="FILE DI INPUT:\n"
        messaggio+="Vettoriale sorgente: "+str(ef.text_vector)+"\n"
        messaggio+="Vettoriale confine: "+str(ef.text_area)+"\n"

        messaggio+="VARIABILI:\n"
        messaggio+="Sostanza inquinante: "+str(sostanza)+"\n"
        messaggio+=u"Concentrazione inquinante: "+str(ef.listaconc[contatore_sostanza])+"\n"
        messaggio+="Tessitura del suolo: "+str(ef.text_texture)+"\n\n"
        messaggio+=u"Densità del suolo: "+str(ef.ro)+"\n"
        messaggio+="Spessore sorgente: "+str(ef.dz)+"\n"
        messaggio+=u"Velocità Darcy: "+str(ef.ve)+"\n"
        messaggio+=u"Profondità mixed zone: "+str(ef.dgw)+"\n"
        messaggio+=u"Estensione orizzontale: "+str(ef.sw)+"\n"
        messaggio+="Coefficiente di decadimento: "+str(ef.lambda1)+"\n"
        messaggio+=u"Profondità acquifero: "+str(ef.lf)+"\n"

        messaggio+=u"Precipitazione media annua: "+str(ef.pioggia)+"\n"
        messaggio+=u"Direzione corrente falda: "+str(ef.dirfalda)+"\n"
        messaggio+="Risoluzione: "+str(ef.res)+"\n\n"
        messaggio+="Nome mappa prodotta: "+str(os.path.basename(ef.path_output+str(sostanza)))+"\n\n"
        messaggio+="ALGORITMI UTILIZZATI: \n"
        messaggio+='Leaching: Calcolo del fattore di Lisciviazione (Agenzia per la Protezione dell’Ambiente. "Criteri metodologici per l\'applicazione dell\'analisi assoluta di rischio ai siti contaminati." (2008).\n\n'
        messaggio+='DAF: '+str(ef.algoritmo)+' (Domenico P.A. e Schwartz F.W. (1998), Physical and Chemical Hydrogeology, John Wiley and Sons, New York)\n\n'
        messaggio+="---------------------------\n\n"
     """

     """ CHECK SUI VALORI DI `res_fields` E `res_texture`
        #controllo se i valori sono presenti
        if '' in res_fields or '' in res_texture:
            messaggio="Analisi non effettuata - modificare i parametri"
            ef.console.appendPlainText(messaggio)
            return
     """
        #                                            ro,    dz,          lf,            ve,    dgw
        element = Leach( h, tera_w, tera_a, kd, ief, soild, sourcethick, aquifer_depth, darcy, Δgw, sw ) 
        kw = calc_kw!(element)
        ldf = calc_ldf!(element)
        sam = calc_sam!(element)
        LF = calc_LF!(element)

     """ PRINT DI COSE
        risultato="RISULTATI INTERMEDI:\n\n"
        risultato+="kw = %f \n" % ef.kw

        risultato+="ldf = %f \n" % ef.ldf

        #risultato+="sam = %f \n" % ef.sam
     """
        c0 = concentration * LF
     """ PRINT DI COSE
        risultato+="fattore di lisciviazione : %f \n" % ef.LF
        risultato+="concentrazione sorgente secondaria : %f mg/l \n" % ef.c0
     """

     # NON SO COME FARE LA STESSA OPERAZIONE SUVETTORIALI CON ArchGDAL
        path_layer = area.dataProvider().dataSourceUri()
        path = path_layer.split("|")
        source_ds = ogr.Open(path[0])
        area_layer = source_ds.GetLayer()
        x_min, y_min, x_max, y_max = round.( Int64, area_layer.GetExtent() )

     # DOVE VIENE USATA QUESTA VARIABILE?
        mem_driver = agd.getdriver("MEM")
        valNoData = -9999

        # Create the destination data source
        x_res = (x_max - x_min) / resolution
        y_res = (y_max - y_min) / resolution

        gtiff_driver = agd.getdriver("GTiff")
        target_ds = agd.create( output_path, gtiff_driver, round(Int64, x_res), round(Int64, y_res), 1, agd.GDAL.GDT_Float32 )
        agd.setgeotransform!(target_ds, [ x_min, resolution, 0.0, y_max, 0.0, -resolution ])
     # DOVE VIENE USATA QUESTA VARIABILE?
        projectionfrom = agd.getproj(target_ds)
        agd.setproj!( target_ds, agd.importEPSG(agd.fromWKT(refsys)) )
     """ NON SO QUALE SIA IL COMANDO PER SETTARE I METADATI CON `ArchGDAL`
        target_ds.SetMetadata({'credits':'Envifate - Francesco Geri, Oscar Cainelli, Paolo Zatelli, Gianluca Salogni, Marco Ciolli - DICAM Università degli Studi di Trento - Regione Veneto',
                                'modulo':'Dispersione in falda',
                                'descrizione':'Simulazione di dispersione inquinante in falda',
                                'srs':ef.source.crs().authid(),
                                'data':datetime.datetime.now().strftime("%d-%m-%y")})
     """
        band1 = agd.getband(target_ds, 1)
        agd.setnodatavalue!( band1, Float64(valNoData) )
        band = agd.read(band)
        fill!(band, valNoData)
        xsize = agd.width(band)
        ysize = agd.height(band)
        # outData = deepcopy(band)
        outData = band
     """ OPERAZIONI SUL VETTORIALE `area`
        polygons = [feature for feature in ef.areastudio.getFeatures()]

        feature = next(ef.source.getFeatures())
        geom = feature.geometry().asPoint()
        x_source=geom[0]
        y_source=geom[1]
     """

        #   original_point = QgsPoint(x_source,y_source)
        original_point = ( x_source, y_source )
        list_result = []
        sw /= 1000
        
        # DA DOVE SI OTTIENE `max`
        for length in resolution:resolution:1000
            cosa = sin(deg2rad(azimut))
            cosb = cos(deg2rad(azimut))

         # DOVE VIENE USATA?
            #   end_point = QgsPoint(original_point.x()+(length*cosa), original_point.y()+(length*cosb))
            end_point = original_point + length * (cosa, cosb)

            calcolo_daf = DAF( c0, 100length, 0, 0, 0, λ1, darcy, kd, soild, tera_e, 1, time )
            cfinal_p = c0 * calc_DAF_ispra!(calcolo_daf)
            push!( list_result, cfinal_p )
        end

        for row in 1:ysize
            for col in 1:xsize
                x, y = (col, row) * resolution + (x_min, y_min) + (resolution/2)
                punto_controllo = ( x, y )

                """ -- QUESTA PARTE NON E' ANCORA SISTEMATA -- """
                for pol in polygons:
                    poly = pol.geometry()
                    if poly.contains(punto_controllo):
                """ ------------------------------------------ """
                        Δx = x - x_source
                        Δy = y - y_source
                        dist = √( Δy^2 + Δx^2 )
                        xvero = Δx * cos(deg2rad(azimut)) - Δy * sin(deg2rad(azimut))
                        yvero = Δx * sin(deg2rad(azimut)) + Δy * cos(deg2rad(azimut))

                        if xvero > 0
                            element_daf = DAF( c0, xvero, yvero, 0, 0, λ1, darcy, kd, soild, tera_e, sw, time )
                            if algorithm == :fickian
                                if opzione == :pulse
                                    cfinal = calc_DAF!(element_daf)
                                else opzione == :continuous
                                    cfinal = calc_DAF_c!(element_daf)
                                end
                            else
                                cfinal = c0 * calc_DAF_ispra!(element_daf)
                            end
                            outData[row, col] = cfinal
                        else
                            outData[row, col] = 0
                        end
                    end
                end
            end
            #   outData_raster = outData [::-1]
            #   band.WriteArray(outData_raster)
            #   astats=band.GetStatistics(0, 1)
            
            # band = deepcopy(outData)
        end

        band = nothing
        target_ds = nothing

        base_raster_name=os.path.basename(ef.path_output+str(sostanza))
        raster_name=os.path.splitext(base_raster_name)[0]
        ef.iface.addRasterLayer(ef.path_output, raster_name)

        contatore_sostanza=0
        ef.list_result=[]

     """
        tempostimato=time.strftime("%H:%M:%S", time.gmtime(tempoanalisi))
        messaggio="---------------------------------\n"
        messaggio+="Fine modellazione\n"
        messaggio+="\nTempo di analisi: "+tempostimato+"\n"
        messaggio+="---------------------------------\n\n"
        messaggio+="ANALISI STATSTICHE DI BASE\nvalore minimo: "+str(astats[0])+"\n"+"valore massimo: "+str(astats[1])+"\n"+"valore medio: "+str(astats[2])+"\n"+"deviazione standard: "+str(astats[3])
     """
    if ef.multiplesubstance_control==False:
        ax1f1 = ef.figure.add_subplot(111)
        ax1f1.plot(ef.list_result)

        ef.canvas_mat.draw()

        ef.label_status.setText("In attesa di dati")
        ef.label_status.setStyleSheet('color : green; font-weight:bold')

end # module