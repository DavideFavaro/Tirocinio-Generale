import ArchGDAL as agd

using CombinedParsers
using CombinedParsers.Regexp
using DataFrames
using HTTP

include("..\\Library\\Functions.jl")


sat_file = split( @__DIR__ , "\\Porting\\")[1] * "\\Mappe\\sat\\sette_sorelle.shp"
sat = agd.read(sat_file)

layer = agd.getlayer(sat, 0)
geometry = agd.getgeom(collect(layer)[1])

src_crs = agd.getspatialref(geometry)
trg_crs = agd.importEPSG(4326)

#= A QUANTO PER QUESTA SCRITTURA NON FUNZIONA MENTRE QUELLA SOTTO SI
coord_transform = agd.createcoordtrans(x->x, src_crs, trg_crs)
agd.transform!(geometry, coord_transform)
=#
agd.createcoordtrans(src_crs, trg_crs) do transform
    agd.transform!(geometry, transform)
end

#=
geom = agd.toWKT(geometry)[11:end-2]
coords_str = split.( split( geom, ",", keepempty=false ), " ", keepempty=false )
coords = map( xy -> Tuple(tryparse.(Float64, xy)), coords_str )
=#
geom = agd.toWKT(geometry)
aoi = geom[1:7] * replace( geom[9:end], " " => "%20" )




function getAoi( path::AbstractString )::AbstractString
    data = agd.read(path)
    layer = agd.getlayer(data, 0)
    geometry = agd.getgeom(collect(layer)[1])
    src_crs = agd.getspatialref(geometry)
    trg_crs = agd.importEPSG(4326)
 #= A QUANTO PER QUESTA SCRITTURA NON FUNZIONA MENTRE QUELLA SOTTO SI
    coord_transform = agd.createcoordtrans(x->x, src_crs, trg_crs)
    agd.transform!(geometry, coord_transform)
 =#
    agd.createcoordtrans(src_crs, trg_crs) do transform
        agd.transform!(geometry, transform)
    end 
 #=
    geom = agd.toWKT(geometry)[11:end-2]
    coords_str = split.( split( geom, ",", keepempty=false ), " ", keepempty=false )
    coords = map( xy -> Tuple(tryparse.(Float64, xy)), coords_str )
 =#
    geom = agd.toWKT(geometry)
    return geom[1:7] * replace( geom[9:end], " " => "%20" )
end


aoi = getAoi(sat_file)