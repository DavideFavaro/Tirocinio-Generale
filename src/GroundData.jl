module GroundData




#=
colonna	     |   descrizione
------------------------------------------------------------------------------------
uiid	     |   identificatore univoco delle stazioni di qualunque genere di misura
------------------------------------------------------------------------------------
parameter	 |   tipo di parametro misurato
------------------------------------------------------------------------------------
unit	     |   unita di misura (possibilmente SI)
------------------------------------------------------------------------------------
value	     |   valore misurato
------------------------------------------------------------------------------------
date	     |   anno mese giorno ora (UTM)
------------------------------------------------------------------------------------
longitude	 |   longitudine della stazione
------------------------------------------------------------------------------------
latitude	 |   latitudine della stazione
------------------------------------------------------------------------------------
quote	     |   quota stazione || quota stazione più quota misura per il vento
------------------------------------------------------------------------------------
validation	 |   bool (già segnalato dalla stazione)
------------------------------------------------------------------------------------
note	     |   errori e outlayers? altro?
------------------------------------------------------------------------------------
=#



using CSV
using DataFrames

include("./GroundDataAA.jl")
include("./GroundDataER.jl")
include("./GroundDataFVG.jl")
include("./GroundDataL.jl")
include("./GroundDataT.jl")
include("./GroundDataV.jl")


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

columns = [ 
            :uiid,                      # identificatore univoco delle stazioni di qualunque genere di misura
            :parameter,                 # tipo di parametro misurato
            :unit,                      # unita di misura (possibilmente SI)
            :value,                     # valore misurato
            :date,                      # anno mese giorno ora (UTM)
            :longitude,                 # longitudine della stazione
            :latitude,                  # latitudine della stazione
            :quote,                     # quota stazione || quota stazione più quota misura per il vento
            :rel_measurement_height     # quota relativa della misurazioni
            :validation,                # bool (già segnalato dalla stazione)
            :note                       # errori e outlayers? altro?
          ]


function standardize( df::DataFrame )







"""
    getGroundData( filePath::AbstractString="."; regions::Region..., type::Data_Type=METEO, source::Data_Source=STATIONS )

Obtain `type` data of `regions` from `source` and save it as `filePath` 
"""
function getGroundData( targetDirectory::AbstractString="."; regions::Region..., type::Data_Type=METEO, source::Data_Source=STATIONS )
    for region in regions
        res = functions[ Integer(region) ]( type, source )
        



    end
end













end # module