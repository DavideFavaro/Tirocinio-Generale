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


class class_daf(object):
  """docstring for element"""
  def __init__(self, c0,x,y,alfa_x,alfa_y,lambda1,v_e,kd,ro_s,tera_e,s_w,T):
    self.c0=c0
    self.x=x
    self.y=y
    self.alfa_x=float(alfa_x)
    self.alfa_y=float(alfa_y)
    self.lambda1=lambda1
    self.v_e=float(v_e)
    self.kd=kd
    self.ro_s=ro_s
    self.tera_e=tera_e
    self.s_w=s_w
    self.T=T

   



  def calc_R(self):
    
    self.R=1+(self.kd*(self.ro_s/self.tera_e))

    return self.R


  def calc_DAF_ispra2(self):

    #import ipdb; ipdb.set_trace()

    if self.alfa_x==0:
      self.alfa_x=self.x*0.1
    if self.alfa_y==0:
      self.alfa_y=self.alfa_x/3
    
    #daf1=(self.x/2*self.alfa_x)*(1-math.sqrt(1+((4*self.lambda1*self.alfa_x*self.R)/self.v_e)))
    daf1=math.exp(self.x/(2*self.alfa_x)*0)
    #daf1e=math.exp(daf1)
    daf2=math.erf(self.s_w/(4*math.sqrt(self.alfa_y*self.x)))
    
    #import pdb; pdb.set_trace() 
    self.DAF=daf1*daf2

    return self.DAF


  def calc_DAF(self):

    #import ipdb; ipdb.set_trace()

    if self.alfa_x==0:
      self.alfa_x=(self.x/100)*0.1
    if self.alfa_y==0:
      self.alfa_y=self.alfa_x/3

    dx=self.alfa_x*self.v_e
    dy=self.alfa_y*self.v_e

    daf_a=self.c0/(4*self.tera_e*math.pi*self.T*(math.sqrt(dx*dy)))

    #import ipdb; ipdb.set_trace()

    daf_b=math.exp( -(( math.pow((self.x-(self.v_e*self.T)),2)/(4*dx*self.T) ) + ( math.pow(self.y,2)/(4*dy*self.T) ) ) )

           
     
    self.DAF=daf_a*daf_b

    return self.DAF


  def calc_DAF_uni(self):

    #import ipdb; ipdb.set_trace()

    if self.alfa_x==0:
      self.alfa_x=self.x*0.1

    dx=self.alfa_x*self.v_e

    daf_a=self.c0/(2*self.tera_e*math.sqrt(4*dx*math.pi*self.T)) 


    #import ipdb; ipdb.set_trace()

    daf_b=math.exp( -(( math.pow((self.x-(self.v_e*self.T)),2)/(4*dx*self.T) ) ) )
       
     
    self.DAF=daf_a*daf_b

    return self.DAF


  def calc_DAF_c(self):
    #continuous

    #import ipdb; ipdb.set_trace()

    if self.alfa_x==0:
      self.alfa_x=self.x*0.1
    if self.alfa_y==0:
      self.alfa_y=self.alfa_x/3

    dx=self.alfa_x*self.v_e
    dy=self.alfa_x*self.v_e

    r=math.sqrt(math.pow(self.x,2)+(math.pow(self.y,2)*(dx/dy)) )

    daf_a=self.c0/(4*self.tera_e*math.sqrt(math.pi)*math.sqrt(self.v_e*r)*math.sqrt(dy))


    #import ipdb; ipdb.set_trace()

    daf_b=math.exp( ((self.x-r)*self.v_e)/(2*dx) ) 
       
     
    self.DAF=daf_a*daf_b

    return self.DAF


  def calc_DAF_ispra(self):

    ######################## modello di domenico ###########################
    # vedere appendice C pagina 2 del documento Criteri metodologici per l'applicazione dell'analisi assoluta di rischio ai siti contaminati
    # la formula originale prevede la produttoria delle 3 componenti x,y,z moltiplicata per 1/4
    # eliminando la terza componente dell'asse z ?? necessario moltplicare per 1/2 (quindi 0.5)
    # per verifica vedere Domenico P.A. e Schwartz F.W. (1998), Physical and Chemical Hydrogeology, John Wiley and Sons, New York.
    # da pagina 642 a pag 644
    
    if self.alfa_x==0:
      self.alfa_x=self.x*0.1
    if self.alfa_y==0:
      self.alfa_y=self.alfa_x/3

    R=1+(self.kd*(self.ro_s/self.tera_e))

    
    daf1=0.50*math.exp((self.x/2*self.alfa_x)*(1-math.sqrt(1+((4*self.lambda1*self.alfa_x*R)/self.v_e))))
    #daf1=math.exp((self.x/(2*self.alfa_x)))
    #daf2=math.erf(self.s_w/(4*math.sqrt(self.alfa_y*self.x)))
    daf21=math.erf((self.y+0.5*self.s_w)/(2*math.sqrt(self.alfa_y*self.x)))
    daf22=math.erf((self.y-0.5*self.s_w)/(2*math.sqrt(self.alfa_y*self.x)))
    #daf_prova=math.erf((self.y+0.5*self.s_w)/(2*math.sqrt(self.alfa_y*self.x)))
    daf3=daf21-daf22
    
    #import ipdb; ipdb.set_trace()

    DAF_tot=daf1*daf3

    return DAF_tot


  ####### end class daf 

def fuori():
  # fix_print_with_import
  print("hello");


if __name__ == "__main__":

  #input variables

  parser = argparse.ArgumentParser(description='Two-dimensional DAF (dilution attenuation factor) model.')
  parser.add_argument('-s','--substance', help='Id of pollution element',required=True)
  parser.add_argument('-t','--texture',help='Soil texture', required=True)
  parser.add_argument('-ro','--density',help='Soil density', required=False, default=1.7)
  parser.add_argument('-v','--darcy',help='Darcy velocity', required=False, default=0.000025)
  parser.add_argument('-ax','--alfa_x',help='Longitudinal dispersion index', required=False, default=0.00)
  parser.add_argument('-ay','--alfa_y',help='Trasversal dispersion index', required=False, default=0.00)
  parser.add_argument('-x','--x',help='Coordinate x', required=True, type=int)
  parser.add_argument('-y','--y',help='Coordinate y', required=True,type=int)
  parser.add_argument('-c0','--concentration_source',help='Source concentration', required=True, type=float)
  parser.add_argument('-l','--lambda1',help='First order decadiment index', required=False, default=0.0, type=float)
  parser.add_argument('-sw','--sw',help='Extension x', required=False, default=1,type=float)
  parser.add_argument('-T','--Time',help='Time', required=True, type=float)
  parser.add_argument('-opz','--option',help='0: Domenico-Schwartz; 1: Fickian', required=False, default=0)
  parser.add_argument('-lev','--level',help='0: input mass; 1: continuos; 2: one-dimensional', required=False, default=0)
  args = parser.parse_args()
   
  ## show values ##
  #print ("Sostanza: %s" % args.substance )
  


  

  lst_fields=['koc_kd']
  res_fields=functions.substance_extract(args.substance,lst_fields)

  kd=res_fields[0]

  lst_fields_t=['por_eff','grain']
  res_texture=functions.texture_extract(args.texture,lst_fields_t)



  tera_e=res_texture[0]
  grain=res_texture[1]

  opz=args.option
  lev=args.level

  #controllo se i valori sono presenti
  if '' in res_fields or '' in res_texture:
    # fix_print_with_import
    print("daf not computable")
    sys.exit()


  #def __init__(self, c0,x,y,alfa_x,alfa_y,lamba,v_e,kd,ro_s,tera_e,s_w):
  element=class_daf(args.concentration_source,args.x,args.y,args.alfa_x,args.alfa_y,
              args.lambda1,args.darcy,kd,args.density,tera_e,args.sw,args.Time)

  

  if opz=="1":  
    # fix_print_with_import
    print("Algoritmo: Fickian model")
    if lev=="0":
      daf=element.calc_DAF()

    elif lev=="1":
      daf=element.calc_DAF_c()

    else:
      daf=element.calc_DAF_uni()
  else:
    # fix_print_with_import
    print("Algoritmo: Domenico-Schwartz")
    daf=element.calc_DAF_ispra()

  
  #import ipdb; ipdb.set_trace()

  # fix_print_with_import
  print("concentrazione punto conformit??: %f " % (daf*element.c0))


  #kw=element.calc_kw(args.density,tera_w,kd,h,tera_a)










