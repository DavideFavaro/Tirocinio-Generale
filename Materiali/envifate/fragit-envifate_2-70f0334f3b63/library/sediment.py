#!/usr/bin/python
# -*- coding: utf-8 -*-

#importazione moduli

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

#fine importazione moduli


class sediment(object):
  """docstring for element"""


  def __init__(self,t,Q,h,dx,dy,x,y,V,w,dt,U=0,tide=0):
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

    self.Q=float(Q)
    self.h=int(h)
    self.dx=float(dx)
    self.dy=float(dy)
    self.x=int(x)
    self.y=int(y)
    self.V=float(V)
    self.w=float(w)
    self.dt=int(dt)
    self.tide=int(tide)
    self.U=float(U)
    self.t=t

    if U>0 and tide>0:
      self.omega=2*math.pi/self.tide
    else:
      self.omega=0





  def calc_q(self):

    return self.Q/(4*math.pi*self.h*math.sqrt(self.dx*self.dy))

  def calcolo_e(self,i):
    if self.omega>0:
      self.ew=self.U/(self.omega*(math.cos(math.radians(self.omega))) - (math.cos(math.radians(self.omega*i*self.dt))) )
    else:
      self.ew=0
    

    e1=math.exp(-(( self.x - self.V*(self.t-i*self.dt) + self.ew ) / ( 4*self.dx * (self.t-i*self.dt) ) ))


    e2=math.exp(-( math.pow(self.y,2) / ( 4*self.dy * (self.t-i*self.dt) ) ) - ( (self.w*(self.t-i*self.dt))/self.h ) )



    return e1*e2

    


def fuori():
  # fix_print_with_import
  print("out");


if __name__ == "__main__":

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


  parser = argparse.ArgumentParser(description='Mathematical model to evaluate temporal evolution of dredging-induced turbidity plumes in the far field under oscillatory tidal currents')
  parser.add_argument('-Q','--input', help='Sediment input rate (kg/sec)',type=float,required=True)
  parser.add_argument('-t','--final_time',help='Final temporal step', type=int,required=True)
  parser.add_argument('-hh','--depth',help='Depth', required=True, type=int)
  parser.add_argument('-dt','--delta',help='Interval time of analysis', required=True, type=int)
  parser.add_argument('-dx','--dx',help='Diffusion X component', required=True, type=float)
  parser.add_argument('-dy','--dy',help='Diffusion Y component', required=True, type=float)
  parser.add_argument('-x','--xcoord',help='Target coordinate X', required=True, type=float)
  parser.add_argument('-y','--ycoord',help='Target coordinate Y', required=True, type=float)
  parser.add_argument('-V','--current',help='Average speed current', required=True, type=float)
  parser.add_argument('-w','--speed_sediment',help='Average sediment velocity', required=True, type=float)
  parser.add_argument('-U','--amplitude',help='Amplitude of the oscillatory current', required=False, type=int, default=0)
  parser.add_argument('-tide','--tide',help='Tidal cycle', required=False, type=int, default=0)


  args = parser.parse_args()
   

  # (self,t,Q,h,dx,dy,x0,y0,x,y,V,w,dt,U=0,tide=0):

  element=sediment(args.final_time,args.input,args.depth,
                args.dx,args.dy,args.xcoord,args.ycoord,
                args.current,args.speed_sediment,args.delta,
                args.amplitude,args.tide)  




  q=element.calc_q()


  n=int(element.t/element.dt)

  csum=0

  for i in range(n):
    e=element.calcolo_e(i)
    csum+=e*(1 / (element.t-i*element.dt) )





  cfinal=q*csum*element.dt

  # fix_print_with_import
  print("Sedimento finale (x,y,t) : %f " % cfinal)








