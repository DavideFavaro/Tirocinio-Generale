#!/usr/bin/python
# -*- coding: utf-8 -*-

#importazione moduli

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

#fine importazione moduli


class scs(object):
  """docstring for element"""
  def __init__(self,cn,c):
    self.c=float(c)
    self.cn=int(cn)

   


  def calc_R(self,S):
    
    #import ipdb; ipdb.set_trace()
    
    R=math.pow(self.c, 2)/(self.c+S)
    I=self.c-R


    return R,I


  def calc_S(self):
    
    #import ipdb; ipdb.set_trace()
    
    S=25.4*((1000/int(self.cn))-10)

    return S


  ####### end class daf 

def fuori():
  # fix_print_with_import
  print("hello");


if __name__ == "__main__":

  #input variables

  parser = argparse.ArgumentParser(description='SCS algorithm')
  parser.add_argument('-s','--soil', help='Id of land cover class',required=True)
  parser.add_argument('-c','--concentration',help='Pollution substance concentration', required=True)
  parser.add_argument('-i','--infiltration',help='Infilration level (a,b,c,d)', required=False, default="a")
  
  args = parser.parse_args()
   
  ## show values ##
  #print ("Sostanza: %s" % args.substance )
  #import ipdb; ipdb.set_trace()

  res_fields=functions.cn_extract(args.infiltration,args.soil)

  cn=res_fields[0]

  #controllo se i valori sono presenti
  if '' in res_fields:
    # fix_print_with_import
    print("scs not computable")
    sys.exit()


  #def __init__(self, c0,x,y,alfa_x,alfa_y,lamba,v_e,kd,ro_s,tera_e,s_w):
  element=scs(cn,args.concentration)

  S=element.calc_S()

  R,I=element.calc_R(S)

  # fix_print_with_import
  print("runoff calcolato: %f " % R)
  # fix_print_with_import
  print("infiltrazione: %f " % I)


  #kw=element.calc_kw(args.density,tera_w,kd,h,tera_a)










