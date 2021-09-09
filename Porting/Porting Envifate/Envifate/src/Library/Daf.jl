module Daf

export ClassDaf, #struct
       calc_R!, calc_DAF_ispra!, calc_DAF_ispra2!, calc_DAF!, calc_DAF_uni!, calc_DAF_c!, #functions associated with ClassDaf
       fuori #function

#!/usr/bin/python
# -*- coding: utf-8 -*-

#importazione moduli
"""imports
  from __future__ import print_function
  from builtins import object
  import os, osr, sys, argparse,math

  try:
    import gdal, ogr
  except:
    # fix_print_with_import
    print("librerie gdal/ogr non trovare")

  try:
      import numpy
  except:
    # fix_print_with_import
    print("librerie numpy non trovare")


  import functions
"""

using ArchGDAL
using ArgParse

include("Functions.jl")

#fine importazione moduli


mutable struct ClassDaf
  """docstring for element"""
  c0
  x
  y
  α_x::Float64
  α_y::Float64
  λ1
  v_e::Float64
  kd
  ro_s
  tera_e
  s_w
  T

  R
  DAF
  DAF_tot

  ClassDaf(c0,x,y,α_x,α_y,λ1,v_e,kd,ro_s,tera_e,s_w,T) = new(c0,x,y,float(α_x),float(α_y),λ1,float(v_e),kd,ro_s,tera_e,s_w,T)
end


function calc_R!(c::ClassDaf)
  c.R = 1 + ( c.kd * ( c.ro_s / c.tera_e ) )
  return c.R
end  


function calc_DAF_ispra!(c::ClassDaf)

  ######################## modello di domenico ###########################
  # vedere appendice C pagina 2 del documento Criteri metodologici per l'applicazione dell'analisi assoluta di rischio ai siti contaminati
  # la formula originale prevede la produttoria delle 3 componenti x,y,z moltiplicata per 1/4
  # eliminando la terza componente dell'asse z è necessario moltplicare per 1/2 (quindi 0.5)
  # per verifica vedere Domenico P.A. e Schwartz F.W. (1998), Physical and Chemical Hydrogeology, John Wiley and Sons, New York.
  # da pagina 642 a pag 644
  
  if c.α_x == 0
    c.α_x = c.x * 0.1
  end
  if c.α_y == 0
    c.α_y = c.α_x / 3
  end

  R = 1 + ( c.kd * ( c.ro_s / c.tera_e ) )

  
  daf1 = 0.50 * exp( ( c.x / 2 * c.α_x ) * ( 1 - √( 1 + ( ( 4 * c.λ1 * c.α_x * R ) / c.v_e ) ) ) )
  #daf1 = exp( ( c.x / ( 2 * c.α_x ) ) )
  #daf2 = erf( c.s_w / ( 4 * √( c.α_y * c.x ) ) )
  daf21 = erf( (c.y + 0.5 * c.s_w ) / ( 2 * √( c.α_y * c.x ) ) )
  daf22 = erf( (c.y - 0.5 * c.s_w ) / ( 2 * √( c.α_y * c.x ) ) )
  #daf_prova = erf( ( c.y + 0.5 * c.s_w ) / ( 2 * √( c.α_y * c.x ) ) )
  daf3 = daf21 - daf22
  
  #import ipdb; ipdb.set_trace()

  DAF_tot = daf1 * daf3

  return DAF_tot
end


function calc_DAF_ispra2!(c::ClassDaf)

  #import ipdb; ipdb.set_trace()

  if c.α_x == 0
    c.α_x = c.x * 0.1
  end
  if c.α_y == 0
    c.α_y = c.α_x / 3
  end
  
  #daf1 = ( c.x / 2 * c.α_x ) * ( 1 - √( 1 + ( ( 4 * c.λ1 * c.α_x * c.R ) / c.v_e ) ) )
  daf1 = exp( c.x / ( 2 * c.α_x ) * 0 )
  #daf1e = exp(daf1)
  daf2 = erf( c.s_w / ( 4 * √( c.α_y * c.x ) ) )
  
  #import pdb; pdb.set_trace() 
  c.DAF = daf1 * daf2

  return c.DAF
end


function calc_DAF!(c::ClassDaf)
  #import ipdb; ipdb.set_trace()

  if c.α_x == 0
    c.α_x = ( c.x / 100 ) * 0.1
  end
  if c.α_y == 0
    c.α_y = c.α_x / 3
  end

  dx = c.α_x * c.v_e
  dy = c.α_y * c.v_e

  daf_a = c.c0 / ( 4 * c.tera_e * π * c.T * √(dx * dy) )

  daf_b = exp( -( ( (( c.x - (c.v_e * c.T) )^2) / (4 * dx * c.T) ) + ((c.y^2) / ( 4 * dy *c.T )) ) )   
    
  c.DAF = daf_a * daf_b

  return c.DAF
end


function calc_DAF_uni!(c::ClassDaf)

  if c.α_x == 0
    c.α_x = c.x * 0.1
  end

  dx = c.α_x * c.v_e

  daf_a = c.c0 / ( 2 * c.tera_e * √( 4 * dx * π * c.T ) ) 

  daf_b = exp( -(( ((c.x - (c.v_e * c.T))^2) / ( 4 * dx * c.T ) )) )
      
  c.DAF = daf_a * daf_b

  return c.DAF
end


function calc_DAF_c!(c::ClassDaf)
  #continuous

  if c.α_x == 0
    c.α_x = c.x * 0.1
  end
  if c.α_y == 0
    c.α_y = c.α_x / 3
  end

  dx = c.α_x * c.v_e
  dy = c.α_x * c.v_e

  r = √( (c.x^2) + ( (c.y^2) * ( dx / dy) ) )

  daf_a = c.c0 / ( 4 * c.tera_e * √(π) * √( c.v_e * r ) * √(dy) )

  daf_b = exp( ( ( c.x - r ) * c.v_e ) / ( 2 * dx ) ) 
      
    
  c.DAF = daf_a * daf_b

  return c.DAF
end

####### end class daf 


function fuori()
  # fix_print_with_import
  print("hello")
end


if abspath(PROGRAM_FILE) == @__FILE__

  #input variables

  aps = ArgParseSettings(description="Two-dimensional DAF (dilution attenuation factor) model.")
  @add_arg_table! aps begin
    "-s", "--substance"
      help = "Id of pollution element"
      required = true
    "-t", "--texture"
      help = "Soil texture"
      required = true
    "-ro", "--density"
      help = "Soil density"
      default = 1.7
      required = false
    "-v", "--darcy"
      help = "Darcy velocity"
      default = 0.000025
      required = false
    "-ax", "--α_x"
      help = "Longitudinal dispersion index"
      default = 0.00
      required = false
    "-ay", "--α_y"
      help = "Trasversal dispersion index"
      default = 0.00
      required = false
    "-x", "--x"
      help = "Coordinate x"
      arg_type = Int
      required = true
    "-y", "--y"
      help = "Coordinate y"
      arg_type = Int
      required = true
    "-c0", "--concentration_source"
      help = "Source concentration"
      arg_type = Float64
      required = true
    "-l", "--λ1"
      help = "First order decadiment index"
      arg_type = Float64
      default = 0.0
      required = false
    "-sw", "--sw"
      help = "Extension x"
      arg_type = Float64
      default = 1.0
      required = false
    "-T", "--Time"
      help = "Time"
      arg_type = Float64
      required = true
    "-opz", "--option"
      help = "0: Domenico-Schwartz; 1: Fickian"
      default = 0
      required = false
    "-lev", "--level"
      help = "0: input mass; 1: continuos; 2: one-dimensional"
      default = 0
      required = false
  end
  args = parse_args( aps )
   
  ## show values ##
  #println("Sostanza: $(args.substance)" )
  

  lst_fields = ["koc_kd"]
  res_fields = Functions.substance_extract( args.substance, lst_fields )

  kd=res_fields[0]

  lst_fields_t = ["por_eff","grain"]
  res_texture = Functions.texture_extract( args.texture, lst_fields_t )



  tera_e = res_texture[0]
  grain = res_texture[1]

  opz = args.option
  lev = args.level

  #controllo se i valori sono presenti
  if isempty(res_fields) || isempty(res_texture)
    # fix_print_with_import
    println("daf not computable")
    exit()
  end

  #def __init__(self, c0,x,y,α_x,α_y,lamba,v_e,kd,ro_s,tera_e,s_w):
  element = ClassDaf( args.concentration_source, args.x,args.y, args.α_x, args.α_y,
                      args.λ1, args.darcy, kd, args.density, tera_e, args.sw,args.Time )

  

  if opz == "1"  
    # fix_print_with_import
    println("Algoritmo: Fickian model")

  """ Extended version
    if lev == "0"
      daf = calc_DAF!(element)
    elseif lev == "1"
      daf = calc_DAF_c!(element)
    else
      daf = calc_DAF_uni!(element)
    end
  """
    daf = (lev == 0) ? calc_DAF!(element) : (lev == 1) ? calc_DAF_c!(element) : calc_DAF_uni!(element)

  else
    # fix_print_with_import
    println("Algoritmo: Domenico-Schwartz")
    daf = calc_DAF_ispra!(element)
  end
  
  #import ipdb; ipdb.set_trace()

  # fix_print_with_import
  println( "concentrazione punto conformità: $( daf * element.c0 ) ")


  #kw = calc_kw(element, args.density, tera_w, kd, h, tera_a )
end

end # module