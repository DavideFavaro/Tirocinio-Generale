module PlumeLoop

export Plume, #struct
       calc_h!, calc_sigma!, calc_g!, calc_C!, #functions associated with "Plume"
       fuori #function

#!/usr/bin/python
# -*- coding: utf-8 -*-

#importazione moduli
"""imports
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

include("Functions.jl")

#fine importazione moduli

mutable struct Plume
  """docstring for element"""
  # q: concentrazione inquinante m3/sec
  # d: distanza (coordinata x)
  # y,z coordinate in metri
  # c_stability: classe di stabilità
  # u: velocità del vento nella direzione x
  # h: altezza camino
  # v_s: velocità gas
  # d_s: diametro camino
  # t_s: temperatura gas all'uscita del camino
  # t_a: temperatura ambiente

  q::Float64
  d::Int
  y::Int
  z::Int
  c_stability
  u::Float64
  h::Float64
  v_s::Float64
  d_s::Float64
  t_s::Float64
  t_a::Float64
  x_w::Int
  outdoor

  H
  sigmay
  sigmaz
  g1
  g2

  Plume(q,d,y,z,c_stability,u,h,v_s,d_s,t_s,t_a,x_w,outdoor) = new(float(q),int(d),int(y),int(z),c_stability,outdoor,float(u),float(v_s),float(d_s),float(t_s),float(t_a),float(h),int(x_w))
end


function calc_h!( p::Plume )
  try
    fb = 9.81 * ( (p.d_s * p.v_s) / 4 ) * ( ( p.t_s / p.t_a ) / p.t_s )

    delta_h = 1.6 * fb^0.333333 * p.d^0.666667

    p.H = p.h + delta_h
  catch
    p.H = p.h
  end
  # import ipdb; ipdb.set_trace()
  return p.H
end


function calc_sigma!( p::Plume )
    sigma_values = Functions.air_extract( p.c_stability, p.outdoor )
    sigmay1 = sigma_values[0]
    sigmay2 = sigma_values[1]
    sigmayexp = sigma_values[2]
    sigmaz1 = sigma_values[3]
    sigmaz2 = sigma_values[4]
    sigmazexp = sigma_values[5]

    p.sigmay = ( sigmay1 * p.d ) / (1 + sigmay2 * p.d )^sigmayexp
    p.sigmaz = ( sigmaz1 * p.d ) / (1 + sigmaz2 * p.d)^sigmazexp

    return p.sigmay, p.sigmaz
end


function calc_g!( p::Plume )
  p.g1 = exp( ( -0.5 * p.y^2 ) / p.sigmay^2 )

  p.g2 = exp( ( -0.5 * (p.z - p.h)^2 ) / p.sigmaz^2 ) + exp( ( -0.5 * (p.z + p.h)^2 ) / p.sigmaz^2 )

  return p.g1, p.g2
end


function calc_C!( p::Plume )
  p.C = ( p.q / (p.u * 3600) ) * ( (p.g1 * p.g2) / ( 2 * π * p.sigmay * p.sigmaz ) )

  return p.C
end


function fuori()
  # fix_print_with_import
  println("out")
end


if abspath(PROGRAM_FILE) == @__FILE__
  # q: concentrazione inquinante m3/sec
  # d: distanza (coordinata x)
  # y,z coordinate in metri
  # c_stability: classe di stabilità
  # u: velocità del vento nella direzione x
  # h: altezza camino
  # v_s: velocità gas
  # d_s: diametro camino
  # t_s: temperatura gas all'uscita del camino
  # t_a: temperatura ambiente


  aps = ArgParseSettings(description="Evaluation of atmospheric chemical concentrations over level terrain resulting from steady emission of a chemical from a smokestack")
  @add_arg_table! aps begin
    "-q", "--concentration"
      help = "Rate of chemical emission"
      arg_type = Float64
      required = true
    "-u", "--wind_spped"
      help = "Wind speed in the main direction"
      arg_type = Float64
      required = true
    "-x", "--distance"
      help = "Distance in wind direction"
      arg_type = Int
      required = true
    "-x_w", "--wind_direction"
      help = "Main wind direction (degree)"
      arg_type = Int
      required = true
    "-y", "--y"
      help = "Y coordinate of target point"
      arg_type = Int
      required = true
    "-z", "--z"
      help = "Height of target point"
      arg_type = Int
      required = true
    "-c", "--class_stability"
      help = "Atmosphere Pasquill class stability"
      required = true
    "-o", "--outdoor"
      help = "Environment (c=country, u=urban)"
      required = true
    "-h_s", "--height"
      help = "Height of stack"
      arg_type = Int
      required = true
    "-v_s", "--gas_velocity"
      help = "Stack gas velocity"
      default = 0.00
      required = false
    "-d_s", "--stack_diameter"
      help = "Diameter of stack"
      default=0.00
      required = false
    "-t_s", "--stack_temperature"
      help = "Absolute stack gas temperature"
      default = 0.00
      required = false
    "-t_a", "--env_temperature"
      help = "Absolute ambient air temperature"
      default = 0.00
      required = false
  end

  args = parse_args( aps )


  # (self,q,d,y,z,c_stability,outdoor,u,h,v_s,d_s,t_s,t_a):


  altezza = 100
  distanza = 10
  tvalue = 0
  tlist = []
  while distanza <= 1000
    altezza += 10
    distanza += 10

    element = Plume( args.concentration, args.distance, distanza, args.z,
                     args.class_stability, args.outdoor, args.wind_spped, args.height,
                     args.gas_velocity, args.stack_diameter, args.stack_temperature,
                     args.env_temperature, args.wind_direction )


    sigmay, sigmaz = calc_sigma!(element)

    g1, g2 = calc_g!(element)

    hvero = calc_h!(element)

    cfinal = calc_C!(element)

    # fix_print_with_import
    println( distanza, cfinal )
  end
end

end # module