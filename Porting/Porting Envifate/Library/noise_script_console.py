from __future__ import print_function
import processing

r_layer=v_layer=None
for lyr in list(QgsMapLayerRegistry.instance().mapLayers().values()):
    if lyr.name() == "dddd":
        r_layer = lyr
        # fix_print_with_import
        print("trovato")
        break
        
for lyr in list(QgsMapLayerRegistry.instance().mapLayers().values()):
    if lyr.name() == "sorgente":
        v_layer = lyr
        # fix_print_with_import
        print("trovato")
        break


for lyr in list(QgsMapLayerRegistry.instance().mapLayers().values()):
    if lyr.name() == "es_veg":
        veg_layer = lyr
        # fix_print_with_import
        print("trovato")
        break

#target_ds=processing.runalg("grass7:r.viewshed",r_layer,"1742773,5138047","1.75","10000",False,"1736295.0,1751615.0,5132537.0,5143122.0",0,None)
#interlyr = iface.addRasterLayer(target_ds['output'], 'intervisibilita')

veg_ds=processing.runalg("grass7:v.to.rast.attribute",veg_layer,0,"categ","1739030.54285,1747656.9909,5135774.94093,5141558.95782",0,-1,0.0001,None)
#veglyr = iface.addRasterLayer(veg_ds['output'], 'rasterizzazione')

veg_ds_layer=QgsRasterLayer(veg_ds['output'],'friction')

cost_ds=processing.runalg("grass7:r.walk",r_layer,veg_ds_layer,v_layer,"0","100","0.72,6.0,1.9998,-1.9998","1.0","-0.2125",True,True,True,"1739000.0,1747700.0,5135700.0,5141600.0",0,-1,0.0001,None)
cost_lyr = iface.addRasterLayer(cost_ds['output'], 'cost distance')

