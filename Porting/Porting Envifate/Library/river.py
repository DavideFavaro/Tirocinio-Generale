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


class river(object):
  """docstring for element"""


  def __init__(self,ma,t,x,dl,v,w,k):
    # element=river(args.concentration,args.time,args.distance,args.fickian,args.velocity)
    self.ma=float(ma)
    self.t=float(t)
    self.x=float(x)
    self.dl=float(dl)
    self.v=float(v)
    self.w=float(w)
    self.k=float(k)




  def calc_concentration(self):
    
    c1=self.x-(self.v*self.t)
    c1_1=-math.pow(c1,2)

    #import ipdb; ipdb.set_trace()
    c2=c1_1/(4*self.dl*self.t)
    c2_1=math.exp(-self.k*self.t)
    c3=math.exp(c2)*c2_1

    c4=(self.ma/self.w)/(math.sqrt(4*math.pi*self.dl*self.t))

    self.C=c4*c3

    #import ipdb; ipdb.set_trace()

    return self.C



def fuori():
  # fix_print_with_import
  print("hello");


if __name__ == "__main__":




  parser = argparse.ArgumentParser(description='Evaluation of the chemical concentration in the river due to an injection of chemical mass')
  parser.add_argument('-Cs','--concentration', help='Concentration of substance',type=float,required=True)
  parser.add_argument('-t','--time',help='Time', type=float,required=True)
  parser.add_argument('-k','--decadiment',help='First order of decadiment', required=False, default=0)
  parser.add_argument('-Dl','--fickian',help='Longitudinal Fickian transport coefficient', required=False, default=0.05)
  parser.add_argument('-x','--distance',help='Distance from secondary source', required=True,type=int)
  parser.add_argument('-v','--velocity',help='Mean velocity of river', required=True, type=float)
  parser.add_argument('-w','--width',help='Mean width of the river', required=False, default=0.05)

  args = parser.parse_args()
   

  element=river(args.concentration,args.time,args.distance,args.fickian,args.velocity,args.width,args.decadiment)  



  C=element.calc_concentration()


  # fix_print_with_import
  print("Concentrazione (x,t) : %f " % C)








