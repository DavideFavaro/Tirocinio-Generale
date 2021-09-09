#!/usr/bin/python
# -*- coding: utf-8 -*-

#importazione moduli

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

#fine importazione moduli


class plume(object):
  """docstring for element"""


  def __init__(self,q,d,y,z,c_stability,outdoor,u,h,v_s,d_s,t_s,t_a,x_w):
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

    self.q=float(q)
    self.d=int(d)
    self.y=int(y)
    self.z=int(z)
    self.c_stability=c_stability
    self.outdoor=outdoor
    self.u=float(u)
    self.v_s=float(v_s)
    self.d_s=float(d_s)
    self.t_s=float(t_s)
    self.t_a=float(t_a)
    self.h=float(h)
    self.x_w=int(x_w)




  def calc_h(self):

    try:
      fb=9.81*((self.d_s*self.v_s)/4)*((self.t_s/self.t_a)/self.t_s)

      delta_h=1.6*math.pow(fb,0.333333)*(math.pow(self.d,0.666667))

      self.H=self.h+delta_h
    except:
      self.H=self.h

    # import ipdb; ipdb.set_trace()
    return self.H

  def calc_sigma(self):
    sigma_values=functions.air_extract(self.c_stability,self.outdoor)
    sigmay1=sigma_values[0]
    sigmay2=sigma_values[1]
    sigmayexp=sigma_values[2]
    sigmaz1=sigma_values[3]
    sigmaz2=sigma_values[4]
    sigmazexp=sigma_values[5]



    self.sigmay=(sigmay1*self.d)/math.pow((1+sigmay2*self.d),sigmayexp)
    self.sigmaz=(sigmaz1*self.d)/math.pow((1+sigmaz2*self.d),sigmazexp)

    return self.sigmay, self.sigmaz

  def calc_g(self):
    self.g1=math.exp((-0.5*math.pow(self.y,2))/math.pow(self.sigmay,2))

    self.g2=math.exp((-0.5*math.pow(self.z-self.h,2))/math.pow(self.sigmaz,2))+\
            math.exp((-0.5*math.pow(self.z+self.h,2))/math.pow(self.sigmaz,2))

    return self.g1,self.g2

  def calc_C(self):

    self.C=(self.q*100/(self.u*3600))*((self.g1*self.g2)/(2*math.pi*self.sigmay*self.sigmaz))


    return self.C

def fuori():
  # fix_print_with_import
  print("out");


if __name__ == "__main__":

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


  parser = argparse.ArgumentParser(description='Evaluation of atmospheric chemical concentrations over level terrain resulting from steady emission of a chemical from a smokestack')
  parser.add_argument('-q','--concentration', help='Rate of chemical emission',type=float,required=True)
  parser.add_argument('-u','--wind_spped',help='Wind speed in the main direction', type=float,required=True)
  parser.add_argument('-x','--distance',help='Distance in wind direction', required=True, type=int)
  parser.add_argument('-x_w','--wind_direction',help='Main wind direction (degree)', required=True, type=int)
  parser.add_argument('-y','--y',help='Y coordinate of target point', required=True, type=int)
  parser.add_argument('-z','--z',help='Height of target point', required=True, type=int)
  parser.add_argument('-c','--class_stability',help='Atmosphere Pasquill class stability', required=True)
  parser.add_argument('-o','--outdoor',help='Environment (c=country, u=urban)', required=True)
  parser.add_argument('-h_s','--height',help='Height of stack', required=True, type=int)
  parser.add_argument('-v_s','--gas_velocity',help='Stack gas velocity', required=False, default=0.00)
  parser.add_argument('-d_s','--stack_diameter',help='Diameter of stack', required=False, default=0.00)
  parser.add_argument('-t_s','--stack_temperature',help='Absolute stack gas temperature', required=False, default=0.00)
  parser.add_argument('-t_a','--env_temperature',help='Absolute ambient air temperature', required=False, default=0.00)

  args = parser.parse_args()


  # (self,q,d,y,z,c_stability,outdoor,u,h,v_s,d_s,t_s,t_a):

  element=plume(args.concentration,args.distance,args.y,args.z,
                args.class_stability,args.outdoor,args.wind_spped,args.height,
                args.gas_velocity,args.stack_diameter,args.stack_temperature,
                args.env_temperature,args.wind_direction)



  sigmay,sigmaz=element.calc_sigma()

  #print "sigmay e sigmaz %f e %f " % (sigmay,sigmaz)

  g1,g2=element.calc_g()

  #print "g1 e g2 %f %f " % (g1,g2)

  hvero=element.calc_h()

  #print "h stack %f mentre h vero è %f " % (args.height,hvero)


  cfinal=element.calc_C()

  # fix_print_with_import
  print("Concentrazione finale (x,y,z) : %f " % cfinal)
