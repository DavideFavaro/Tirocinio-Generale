module Leaching

export Leach, #struct
       calc_kw!, calc_ldf!, calc_LF!, calc_sam!, #functions associated with "Leach"
       fuori #function

#!/usr/bin/python
# -*- coding: utf-8 -*-

#importazione moduli
""" Imports
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

#using Pkg
#Pkg.add("ArgParse")

using ArchGDAL
using ArgParse

include("Functions.jl")

#fine importazione moduli

#class Leach
mutable struct Leach
  h
  tera_w
  tera_a
  kd
  ief
  ro_s
  dz
  lf
  v_e
  delta_gw
  W

  kw
  koc
  ldf
  sam
  LF

  Leach(h,tera_w,tera_a,kd,ief,ro_s,dz,lf,v_e,delta_gw,W) = new(h,tera_w,tera_a,kd,float(ief),ro_s,dz,lf,v_e,float(delta_gw),float(W))
end

#Class methods
#method calc_kw
function calc_kw!( l::Leach )
  foc = 0.01
  l.koc = l.kd
  l.kd = l.koc * foc
  l.kw = l.ro_s / ( l.tera_w + ( l.kd * l.ro_s ) + ( l.h * l.tera_a ) )
  return l.kw
end

#method clac_ldf 
function calc_ldf!( l::Leach )
  v_darcy = float(self.v_e) * 100 * 86400 * 365
  l.ldf = 1 + ( v_darcy * ( l.delta_gw / ( l.ief * l.W ) ) )
  #import ipdb; ipdb.set_trace()
  #import pdb; pdb.set_trace()
  return l.ldf
end

#method calc_sam
function calc_sam!( l::Leach )
  #import ipdb; ipdb.set_trace()
  l.sam = float(l.dz) / float(l.lf)
  return l.sam    
end

#method calc_LF
function calc_LF!( l::Leach )
    #l.LF = (l.kw) / l.ldf
    l.LF = ( l.kw * l.sam ) / l.ldf
    return l.LF
end
# end of class methods

#function fuori
function fuori()
  # fix_print_with_import
  println("hello");
end


if abspath(PROGRAM_FILE) == @__FILE__

  #check input variables
  # -s : id inquinante
  # -ro : densità suolo
  # -t : tessitura
  # -dz : spessore sorgente
  # -lf : profondità falda
  # -v : velocità Darcy
  # -dgw : spessore zona miscelazione falda
  # -W : estensione sorgente

  #parsing possibili valori da command line
  aps = ArgParseSettings( description="One-dimensional leaching script" )
  @add_arg_table aps begin
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
    "-dz", "--depth_source"
      help = "Source depth"
      default = 1
      required = false
    "-lf", "--depth"
      help = "Aquifer depth"
      required = true
    "-v", "--darcy"
      help = "Darcy velocity"
      default = 0.000025
      required = false
    "-dgw", "--depth_mixing"
      help = "Depth of mixing zone"
      default = 1
      required = false
    "-w", "--extent"
      help = "Source extent"
      default = 10000
      required = false
    "-p", "--rainfall"
      help = "Annual rainfall average (mm/d)"
      required = true
  end
  args = parse_args( aps )
  
  ####### importante########
  # inserire pioggia (attenzione unità di misura)


  ## show values ##
  #print ("Sostanza: %s" % args.substance )
  
  pioggia = float(args.rainfall) / 10

  # import ipdb; ipdb.set_trace()


  lst_fields = [ "c_henry","koc_kd" ]
  res_fields = Functions.substance_extract( args.substance, lst_fields )

  h = res_fields[0]
  kd = res_fields[1]



  lst_fields_t = [ "tot_por", "c_water_avg", "ief" ]
  res_texture = Functions.texture_extract(args.texture,lst_fields_t)
    
  #inserire pioggia Ief=ief_tabulato*pioggia(al quadrato)


  tera_a=res_texture[0]
  tera_w=res_texture[1]
  #import ipdb; ipdb.set_trace()
  ief= res_texture[2] * pioggia^2 

  

  #controllo se i valori sono presenti
  #if '' in res_fields or '' in res_texture:
  if isempty(res_fields) || isempty(res_texture)
    # fix_print_with_import
    println("leaching not computable")
    
    #ccall(:jl_exit, Cvoid, (Int32,), 86)
    exit(86)
  end 

  #def __init__(self, h, tera_w,tera_a,ro_s,dz,lf,v_e,delta_gw,W):
  element=Leach(h,tera_w,tera_a,kd,ief,args.density,args.depth_source,
                args.depth,args.darcy,args.depth_mixing,args.extent)  



  kw = calc_kw!(element)
  ldf = calc_ldf!(element)
  sam = calc_sam!(element)
  LF = calc_LF!(element)




  # fix_print_with_import
  println("kw = $kw")

  # fix_print_with_import
  println("ldf = $ldf")

  # fix_print_with_import
  println("sam = $sam")

  # fix_print_with_import
  println("fattore di lisciviazione : $LF ")
end

end # module