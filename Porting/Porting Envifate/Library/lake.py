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


class lake(object):
  """docstring for element"""


  def __init__(self,ma,t,x,y,d_x,d_y,v_x,v_y,k_d):
    """
    element=lake(args.concentration,args.time,args.distance_x,args.distance_y,
               args.fickian_x,args.fickian_y,args.velocity_x,args.velocity_y,
               args.lambdak)  
    """
    self.ma=float(ma)
    self.t=float(t)
    self.x=float(x)
    self.y=float(y)
    self.d_x=float(d_x)
    self.d_y=float(d_y)
    self.v_x=float(v_x)
    self.v_y=float(v_y)
    self.k_d=float(k_d)



  def calc_concentration(self):
    
    # c1=(math.pow(self.x-(self.v_x*self.t),2))/(4*self.d_x*self.t)
    # c2=(math.pow(self.y-(self.v_y*self.t),2))/(4*self.d_y*self.t)

    c1_1=self.x-(self.v_x*self.t)
    c1_2=math.pow(c1_1,2)
    c1_3=c1_2/(4*self.d_x*self.t)


    c2_1=self.y-(self.v_y*self.t)
    c2_2=math.pow(c2_1,2)
    c2_3=c2_2/(4*self.d_y*self.t)

    #c3=math.exp(-(c1+c2))

    c3_2=math.exp(-(c1_3+c2_3))

    c3_1=math.exp(-self.k_d*self.t)

    #c4=c3*c3_1
    c4=c3_2*c3_1

    c5=self.ma/(4*math.pi*self.t*math.sqrt(self.d_x*self.d_y))

    #import ipdb; ipdb.set_trace()
    
    self.C=c4*c5

    return self.C



def fuori():
  # fix_print_with_import
  print("hello");


if __name__ == "__main__":




  parser = argparse.ArgumentParser(description='Two-dimensional model for concentration of a chemical introduced as a pulse over the depth of a vertically mixed lake.')
  parser.add_argument('-Cs','--concentration', help='Concentration of substance per depth of water',type=float,required=True)
  parser.add_argument('-t','--time',help='Time', type=float,required=True)
  # parser.add_argument('-k','--decadiment',help='First order of decadiment', required=False, default=1.7)
  parser.add_argument('-Dx','--fickian_x',help='X direction Fickian transport coefficient', required=False, default=0.05)
  parser.add_argument('-Dy','--fickian_y',help='Y direction Fickian transport coefficient', required=False, default=0.05)
  #attenzione!!!! i valori d_x e d_y sono troppo bassi
  parser.add_argument('-x','--distance_x',help='Distance from secondary source in the X direction', required=True,type=int)
  parser.add_argument('-y','--distance_y',help='Distance from secondary source in the Y direction', required=True,type=int)
  parser.add_argument('-Vx','--velocity_x',help='Mean velocity of river in the X direction', required=True, type=float)
  parser.add_argument('-Vy','--velocity_y',help='Mean velocity of river in the Y direction', required=True, type=float)
  parser.add_argument('-l','--lambdak',help='First order decadiment', type=float,required=False,default=0.0)

  args = parser.parse_args()
   

  element=lake(args.concentration,args.time,args.distance_x,args.distance_y,
               args.fickian_x,args.fickian_y,args.velocity_x,args.velocity_y,
               args.lambdak)  



  C=element.calc_concentration()


  # fix_print_with_import
  print("Concentrazione (x,y,t) : %f " % C)








