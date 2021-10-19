module GroundData



using CSV

include("./Global.jl")
include("./GroundDataV.jl")
include("./GroundDataL.jl")
include("./GroundDataER.jl")
include("./GroundDataFVG.jl")
include("./GroundDataTAA.jl")


export getGroundData


@enum Data_Type METEO=1 AIRQUALITY=2 
@enum Data_Source STATIONS=1 SENSORS=2
@enum Region AA=1, ER=2 FVG=3 L=4 T=5 V=6


functions = [ getDataAA,
              getDataER,
              getDataFVG,
              getDataL,
              getDataT,
              getDataV ]



"""
    getGroundData( filePath::AbstractString="."; regions::Region..., type::Data_Type=METEO, source::Data_Source=STATIONS )

Obtain `type` data of `regions` from `source` and save it as `filePath` 
"""
function getGroundData( filePath::AbstractString="."; regions::Region..., type::Data_Type=METEO, source::Data_Source=STATIONS )
    for region in regions
        res = functions[Integer(region)]( type, source )
        CSV.write( filePath, res )
    end
end

end # module