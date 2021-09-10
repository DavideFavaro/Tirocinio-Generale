module Sediments

export Sediment, #struct
       calc_q, calcolo_e!, #functions associated with "Sediment"
       fuori #function

#!/usr/bin/python
# -*- coding: utf-8 -*-

#importazione moduli
#= imports
  from __future__ import print_function
  from builtins import object
  import os, osr, sys, argparse, math,pdb

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
=#

using ArchGDAL
using ArgParse
using Parameters

#fine importazione moduli

#classe sediment
@with_kw mutable struct Sediment
  """docstring for element"""

  # Q: input sedimento (kg/sec)
  # t: tempo finale
  # h: profondità metri
  # Dx,Dy: coefficienti di diffusione
  # x0,y0: coordinate sorgente
  # x,y: coordinate target point
  # V: velocità media corrente
  # w: velocità di sedimentazione
  # dt: delta t per la discretizzazione dell'integrale
  # U: amplitude of the oscillatory current
  # tide: tidal cycle es. 12 (ore)

  Q::Float64
  t
  h::Int
  dx::Float64
  dy::Float64
  x::Int
  y::Int
  V::Float64
  w::Float64
  dt::Int
  U::Float64 = 0.0
  tide::Int = 0

  ω = 0
  ew

  #costruttore ( metodo __init__ )
  function Sediment(Q,t,h,dx,dy,x,y,V,w,dt,U,tide)
    if U > 0 && tide > 0
      ω = 2*π/tide
      return new(float(Q),t,int(h),float(dx),float(dy),int(x),int(y),float(V),float(w),int(dt),float(U),int(tide),ω)
    end
  end
end

#Class methods
#metodo calc_q
function calc_q( s::Sediment )
  return s.Q / ( 4 * π *s.h * isqrt(s.dx * s.dy) )
end

#metodo calcolo_e
function calcolo_e!( s::Sediment, i )
  s.ew = s.ω > 0 ? s.U/(s.ω*(cos(deg2rad(s.ω))) - (cos(deg2rad(s.ω*i*s.dt))) ) : 0

  e1 = exp(-(( s.x - s.V*(s.t-i*s.dt) + s.ew ) / ( 4*s.dx * (s.t-i*s.dt) ) ))

  e2 = exp(-( (s.y^2) / ( 4*s.dy * (s.t-i*s.dt) ) ) - ( (s.w*(s.t-i*s.dt))/s.h ) )

  return e1*e2
end

#funzione fuori
function fuori()
  # fix_print_with_import
  print("out");
end 



#if __name__ == "__main__":
if abspath(PROGRAM_FILE) == @__FILE__

  # Q: input sedimento (kg/sec)
  # t: tempo finale
  # h: profondità metri
  # Dx,Dy: coefficienti di diffusione
  # x0,y0: coordinate sorgente
  # x,y: coordinate target point
  # V: velocità media corrente
  # w: velocità di sedimentazione
  # dt: delta t per la discretizzazione dell'integrale
  # U: amplitude of the oscillatory current
  # tide: tidal cycle es. 12 (ore)


  aps = ArgParseSettings(description="Mathematical model to evaluate temporal evolution of dredging-induced turbidity plumes in the far field under oscillatory tidal currents")
  @add_arg_table! aps begin
    "-Q", "--input"
      help = "Sediment input rate (kg/sec)"
      arg_type = Float64
      required = true
    "-t", "--final_time"
      help = "Final temporal step"
      arg_type = Int
      required = true
    "-hh", "--depth"
      help = "Depth"
      arg_type = Int
      required = true
    "-dt", "--delta"
      help = "Interval time of analysis"
      arg_type = Int
      required = true
    "-dx", "--dx"
      help = "Diffusion X component"
      arg_type = Float64
      required = true
    "-dy", "--dy"
      help = "Diffusion Y component"
      arg_type = Float64
      required = true
    "-x", "--xcoord"
      help = "Target coordinate X"
      arg_type = Float64
      required = true
    "-y", "--ycoord"
      help = "Target coordinate Y"
      arg_type = Float64
      required = true
    "-V", "--current"
      help = "Average speed current"
      arg_type = Float64
      required = true
    "-w", "--speed_sediment"
      help = "Average sediment velocity"
      arg_type = Float64
      required = true
    "-U", "--amplitude"
      help = "Amplitude of the oscillatory current"
      arg_type = Int
      default = 0
      required = false
    "-tide", "--tide"
      help = "Tidal cycle"
      arg_type = Int
      default = 0
      required = false
  end
  args = parse_args( aps )

  # (self,t,Q,h,dx,dy,x0,y0,x,y,V,w,dt,U=0,tide=0):

  element = Sediment( args.final_time, args.input, args.depth,
                      args.dx, args.dy, args.xcoord, args.ycoord,
                      args.current, args.speed_sediment, args.delta,
                      args.amplitude, args.tide )  

  q = calc_q(element)

  n = int(element.t/element.dt)

  csum = 0

  for i in range( 1, n, step = 1 )
    e = calcolo_e!(element, i)
    csum += e * (1 / (element.t-i*element.dt) )
  end

  cfinal = q*csum*element.dt

  # fix_print_with_import
  println("Sedimento finale (x,y,t) : %f " % cfinal)
end

end # module