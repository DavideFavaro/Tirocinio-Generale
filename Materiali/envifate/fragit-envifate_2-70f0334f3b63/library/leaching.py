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


class leach(object):
  """docstring for element"""


  def __init__(self,h,tera_w,tera_a,kd,ief,ro_s,dz,lf,v_e,delta_gw,W):
    self.h=h
    self.tera_w=tera_w
    self.tera_a=tera_a
    self.ro_s=ro_s
    self.dz=dz
    self.kd=kd
    self.ief=float(ief)
    self.v_e=v_e
    self.delta_gw=float(delta_gw)
    self.W=float(W)
    self.lf=lf



  def calc_kw(self):
    foc=0.01
    self.koc=self.kd
    self.kd=self.koc*foc
    self.kw=self.ro_s/(self.tera_w+(self.kd*self.ro_s)+(self.h*self.tera_a))
    

    return self.kw


  def calc_ldf(self):
    

    v_darcy=float(self.v_e)*100*86400*365
    self.ldf=1+(v_darcy*(self.delta_gw/(self.ief*self.W)))
    #import ipdb; ipdb.set_trace()
    #import pdb; pdb.set_trace()
    return self.ldf


  def calc_sam(self):
    #import ipdb; ipdb.set_trace()
    
    self.sam=float(self.dz)/float(self.lf)

    return self.sam    


  def calc_LF(self):
    
    #self.LF=(self.kw)/self.ldf
    self.LF=(self.kw*self.sam)/self.ldf

    return self.LF 


def fuori():
  # fix_print_with_import
  print("hello");


if __name__ == "__main__":

  #check input variables
  # -s : id inquinante
  # -ro : densità suolo
  # -t : tessitura
  # -dz : spessore sorgente
  # -lf : profondità falda
  # -v : velocità Darcy
  # -dgw : spessore zona miscelazione falda
  # -W : estensione sorgente


  parser = argparse.ArgumentParser(description='One-dimensional leaching script')
  parser.add_argument('-s','--substance', help='Id of pollution element',required=True)
  parser.add_argument('-t','--texture',help='Soil texture', required=True)
  parser.add_argument('-ro','--density',help='Soil density', required=False, default=1.7)
  parser.add_argument('-dz','--depth_source',help='Source depth', required=False, default=1)
  parser.add_argument('-lf','--depth',help='Aquifer depth', required=True)
  parser.add_argument('-v','--darcy',help='Darcy velocity', required=False, default=0.000025)
  parser.add_argument('-dgw','--depth_mixing',help='Depth of mixing zone', required=False, default=1)
  parser.add_argument('-w','--extent',help='Source extent', required=False, default=10000)
  parser.add_argument('-p','--rainfall',help='Annual rainfall average (mm/d)', required=True)
  args = parser.parse_args()
   
  ####### importante########
  # inserire pioggia (attenzione unità di misura)


  ## show values ##
  #print ("Sostanza: %s" % args.substance )
  
  pioggia=float(args.rainfall)/10

  # import ipdb; ipdb.set_trace()


  lst_fields=['c_henry','koc_kd']
  res_fields=functions.substance_extract(args.substance,lst_fields)

  h=res_fields[0]
  kd=res_fields[1]

  lst_fields_t=['tot_por','c_water_avg','ief']
  res_texture=functions.texture_extract(args.texture,lst_fields_t)
    
  #inserire pioggia Ief=ief_tabulato*pioggia(al quadrato)


  tera_a=res_texture[0]
  tera_w=res_texture[1]
  #import ipdb; ipdb.set_trace()
  ief=res_texture[2]*math.pow(pioggia,2)

  

  #controllo se i valori sono presenti
  if '' in res_fields or '' in res_texture:
    # fix_print_with_import
    print("leaching not computable")
    sys.exit()

  #def __init__(self, h, tera_w,tera_a,ro_s,dz,lf,v_e,delta_gw,W):
  element=leach(h,tera_w,tera_a,kd,ief,args.density,args.depth_source,
                args.depth,args.darcy,args.depth_mixing,args.extent)  



  kw=element.calc_kw()
  ldf=element.calc_ldf()
  sam=element.calc_sam()

  LF=element.calc_LF()




  # fix_print_with_import
  print("kw = %f" % kw)

  # fix_print_with_import
  print("ldf = %f" % ldf)

  # fix_print_with_import
  print("sam = %f" % sam)

  # fix_print_with_import
  print("fattore di lisciviazione : %f " % LF)








