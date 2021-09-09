module Scs

using ArchGDAL: isempty
export ScsClass, #struct
       calc_R, calc_S, #functions associated with ScsClass
       fuori #function

#!/usr/bin/python
# -*- coding: utf-8 -*-

#importazione moduli
""" imports
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

mutable struct ScsClass
  """docstring for element"""

  c::Float64
  cn::Int

  ScsClass( c, cn ) = new( float(c), int(cn) )
end


function calc_R( scs::ScsClass, S )

  #import ipdb; ipdb.set_trace()

  R = (scs.c^2) / (scs.c + S)
  I = scs.c - R

  return R,I
end


function calc_S(scs::ScsClass)
    
  #import ipdb; ipdb.set_trace()

  S = 25.4 * ( ( 1000 / scs.cn ) - 10 )

  return S
end

####### end class daf 


function fuori()
  # fix_print_with_import
  print("hello")
end


if abspath(PROGRAM_FILE) == @__FILE__

  #input variables


  aps = ArgParseSettings( description = "SCS algorithm" )
  @add_arg_table! aps begin
    "-s","--soil"
      help = "Id of land cover class"
      required = true
    "-c","--concentration"
      help = "Pollution substance concentration"
      required = true
    "-i","--infiltration"
      help = "Infilration level (a,b,c,d)"
      default = "a"
      required = false
  end
  args = parse_args(aps)

  ## show values ##
  #print ("Sostanza: %s" % args.substance )
  #import ipdb; ipdb.set_trace()

  res_fields = Functions.cn_extract( args.infiltration, args.soil )

  cn = res_fields[0]

  #controllo se i valori sono presenti
  if isempty(res_fields)
    # fix_print_with_import
    print("scs not computable")
    exit()
  end

  #def __init__(self, c0,x,y,alfa_x,alfa_y,lamba,v_e,kd,ro_s,tera_e,s_w):
  element = ScsClass( args.concentration, cn )

  S = calc_S(element)

  R, I = calc_R( element, S )

  # fix_print_with_import
  print("runoff calcolato: %f " % R)
  # fix_print_with_import
  print("infiltrazione: %f " % I)


  #kw = calc_kw( element, args.density, tera_w, kd, h, tera_a )
end
end # module