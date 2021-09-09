module Noises

export find_distances, find_coords, seq #functions

#!/usr/bin/python
# -*- coding: utf-8 -*-  
# tutto ruota attorno al cost distance
# in pratica, dopo numerosi test e vedendo come affronta il problema spreadGIS
# si vede che la chiave è la costa distance con il friction surface calcolata come 
# riduzione del livello sonoro.
# spreadGIS fa così: riclassifica la mappa dell'uso del suolo utilizzando dei valori che sono i coefficienti di riduzione
# poi fa la sottrazione tra cost distance pulita e costa distance con le barriere    
# ALGORITHM|2017-04-26 14:49:11|processing.runalg("grass7:v.to.rast.attribute","/home/francescogeri/Desktop/rumore/es_veg.shp",0,"categ","1739030.54285,1747656.9909,5135774.94093,5141558.95782",0,-1,0.0001,None)
# ALGORITHM|2017-04-26 14:53:28|processing.runalg("grass7:v.to.rast.attribute","/home/francescogeri/Desktop/rumore/es_veg.shp",0,"categ","1739030.54285,1747656.9909,5135774.94093,5141558.95782",0,-1,0.0001,"/home/francescogeri/Desktop/rrr.tif")
# ALGORITHM|2017-04-26 14:53:46|processing.runalg("grass7:r.walk","/home/francescogeri/Desktop/rumore/dddd.tif","/home/francescogeri/Desktop/rrr.tif","/home/francescogeri/Desktop/rumore/sorg_buf.shp","0","100","0.72,6.0,1.9998,-1.9998","1.0","-0.2125",True,True,True,"1739000.0,1747700.0,5135700.0,5141600.0",0,-1,0.0001,None)   

#importazione moduli
"""imports
    from __future__ import print_function
    from builtins import range
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

#fine importazione moduli

function find_distances(max_dist, npts, n_digits = 8) #offset = 0,
    increment = max_dist * (npts - 1)^-1 # - 1 accounts for the inclusion of the end points, **-1 does division
    increment = round(increment, n_digits)   
    # If the increment is 0, then just create a vector of 0's
    #if round(increment, n_digits) == 0:
    if increment == 0
        dist_vec = [0] * npts
    else
        # create dist_vec
        dist_vec = [ 0, max_dist, increment ]
        #**# Below is no longer necessary with change in where rounding occurs
        #for i in range(len(dist_vec)):
        #    dist_vec[i] = round(dist_vec[i], n_digits)
    
        # check that the final point was added to dist_vec (if dist_vec is one short, add the final point)
        # Sometimes there is a mismatch due to rounding, and the last point is not added.
        if len(dist_vec) == (npts - 1)
            append!(dist_vec, max_dist)
        end
    end     
    return dist_vec
end

function seq(start, stop, by)
    #**# For now, this only works for ascending sequences    
    # Test if stop is greater than or less than start -     
    if stop < start
        throw( DomainError("'stop' must be greater than 'start'") )
    end
    # Check if by is positive
    if by <= 0
        throw( DomainError("'by' must be positive") )
    end  
    my_lst = [ value for value in range(start, stop, step = by) ]    
    return my_lst
end

function find_coords(source, receiver, npts)
    direction = source > receiver ? -1 : 1
    
    abs_diff = abs(receiver - source)
    # Get distances between source & receiver in this dimension    
    coords_1d = find_distances(abs_diff, npts)  # 0 corresponds to the offset input.     
    # Convert distances to coordinates        
    for i in range( 1, length(coords_1d), step = 1 )
        this_dist = coords_1d[i]
        # Start at source location, then add distance in the appropriate direction
        new_dist = source + this_dist * direction
        coords_1d[i] = new_dist
    end  
    return coords_1d
end

end # module