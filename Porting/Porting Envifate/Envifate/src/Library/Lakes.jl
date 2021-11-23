module Lakes

export Lake, #struct
       calc_concentration!, #functions associated with "Lake"
       fuori #function


#!/usr/bin/python
# -*- coding: utf-8 -*-

#importazione moduli
#=imports
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
=#

using ArchGDAL
using ArgParse

#fine importazione moduli

mutable struct Lake
  """docstring for element"""

  ma::Float64
  t::Float64
  x::Float64
  y::Float64
  d_x::Float64
  d_y::Float64
  v_x::Float64
  v_y::Float64
  k_d::Float64

  C


  #=
    element=Lake(args.concentration,args.time,args.distance_x,args.distance_y,
              args.fickian_x,args.fickian_y,args.velocity_x,args.velocity_y,
              args.lambdak)  
  =#
  Lake( ma, t, x, y, d_x, d_y, v_x, v_y, k_d ) = new( ma, t, x, y, d_x, d_y, v_x, v_y, k_d )
end


function calc_concentration!( l::Lake )
  #c1 = (l.x - (l.v_x * l.t))^2 / ( 4 * l.d_x * l.t )
  #c2 = (l.y - (l.v_y * l.t))^2 / ( 4 * l.d_y * l.t )

  c1_1 = l.x - ( l.v_x * l.t )
  c1_2 = c1_1^2
  c1_3 = c1_2 / ( 4 * l.d_x * l.t )

  c2_1 = l.y - ( l.v_y * l.t )
  c2_2 = c2_1^2
  c2_3 = c2_2 / ( 4 * l.d_y * l.t )

  #c3 = exp( -(c1 + c2) )
  c3_2 = exp( -(c1_3 + c2_3) )
  c3_1 = exp( -l.k_d * l.t )

  #c4 = c3 * c3_1
  c4 = c3_2 * c3_1

  c5 = l.ma / ( 4 * π * l.t * √(l.d_x * l.d_y) )

  #import ipdb; ipdb.set_trace()  
  l.C = c4 * c5
  return l.C
end


function fuori()
  # fix_print_with_import
  println("hello")
end



#if __name__ == "__main__":
if abspath(PROGRAM_FILE) == @__FILE__
  aps = ArgParseSettings(description="Two-dimensional model for concentration of a chemical introduced as a pulse over the depth of a vertically mixed lake.")
  @add_arg_table! aps begin
    "-Cs","--concentration"
      help = "Concentration of substance per depth of water"
      arg_type = Float64
      required = true
    "-t","--time"
      help = "Time"
      arg_type = Float64
      required = true
    # "-k","--decadiment" help = "First order of decadiment" required = false default=1.7
    "-Dx","--fickian_x"
      help = "X direction Fickian transport coefficient"
      default = 0.05
      required = false
    "-Dy","--fickian_y"
      help = "Y direction Fickian transport coefficient"
      default = 0.05
      required = false
    #attenzione!!!! i valori d_x e d_y sono troppo bassi
    "-x","--distance_x"
      help = "Distance from secondary source in the X direction"
      arg_type = Int
      required = true
    "-y","--distance_y"
      help = "Distance from secondary source in the Y direction"
      arg_type = Int
      required = true
    "-Vx","--velocity_x"
      help = "Mean velocity of river in the X direction"
      arg_type = Float64
      required = true
    "-Vy","--velocity_y"
      help = "Mean velocity of river in the Y direction"
      arg_type = Float64
      required = true
    "-l","--lambdak"
      help = "First order decadiment"
      arg_type = Float64
      default=0.0
      required = false
  end
  args = parse_args( aps )

  element = Lake( args.concentration, args.time, args.distance_x, args.distance_y,
                  args.fickian_x, args.fickian_y, args.velocity_x ,args.velocity_y,
                  args.lambdak )  

  C = calc_concentration!( element )


  # fix_print_with_import
  println("Concentrazione (x,y,t) : $C ")
end

end # module