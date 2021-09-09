module Rivers

export River, #struct
       calc_concentration!, #function associated with River
       fuori #function

#!/usr/bin/python
# -*- coding: utf-8 -*-


#importazione moduli
""" imports
  from __future__ import print_function
  from builtins import object
  import os, osr, sys, argparse, math

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

#fine importazione moduli


mutable struct River
  """docstring for element"""

  ma::Float64
  t::Float64
  x::Float64
  dl::Float64
  v::Float64
  w::Float64
  k::Float64

  C

  # element=river(args.concentration,args.time,args.distance,args.fickian,args.velocity)
  River( ma, t, x, dl, v, w, k ) = new( float(ma), float(t), float(x), float(dl), float(v), float(w), float(k) )
end

function calc_concentration!( r::River )   
  c1 = r.x - (r.v * r.t)
  c1_1 = -(c1^2)

  #import ipdb; ipdb.set_trace()
  c2 = c1_1 / ( 4 * r.dl * r.t )
  c2_1 = exp( -r.k * r.t )
  c3 = exp(c2) * c2_1

  c4 = ( r.ma / r.w ) / ( √( 4 * π * r.dl * r.t ) )

  r.C = c4 * c3

  #import ipdb; ipdb.set_trace()

  return r.C
end


function fuori()
 # fix_print_with_import
 println("hello")
end




if abspath(PROGRAM_FILE) == @__FILE__
  aps = ArgParseSettings(description = "Evaluation of the chemical concentration in the river due to an injection of chemical mass")
  @add_arg_table! aps begin
    "-Cs","--concentration"
      help = "Concentration of substance"
      arg_type = Float64
      required = true
    "-t","--time"
      help = "Time"
      arg_type = Float64
      required = true
    "-k","--decadiment"
      help = "First order of decadiment"
      default = 0
      required = false
    "-Dl","--fickian"
      help = "Longitudinal Fickian transport coefficient"
      default = 0.05
      required = false
    "-x","--distance"
      help = "Distance from secondary source"
      arg_type = Int
      required = true
    "-v","--velocity"
      help = "Mean velocity of river"
      arg_type = Float64
      required = true
    "-w","--width"
      help = "Mean width of the river"
      default = 0.05
      required = false
  end
  args = parse_args( aps )

  element = River( args.concentration, args.time, args.distance, args.fickian, args.velocity, args.width, args.decadiment )  

  C = calc_concentration!(element)


  # fix_print_with_import
  println( "Concentrazione (x,t) : $C " )
end
end # module