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


str = occursin( "GroundData.jl", @__FILE__ ) ? "" : "src\\"
include("$(@__DIR__)\\$(str)GroundDataAA.jl")
include("$(@__DIR__)\\$(str)GroundDataER.jl")
include("$(@__DIR__)\\$(str)GroundDataFVG.jl")
include("$(@__DIR__)\\$(str)GroundDataL.jl")
include("$(@__DIR__)\\$(str)GroundDataT.jl")
include("$(@__DIR__)\\$(str)GroundDataV.jl")


export getGroundData


@enum Region AA=1 ER=2 FVG=3 L=4 T=5 V=6


const regions_modules = [
                      GroundDataAA,
                      GroundDataER,
                      GroundDataFVG,
                      GroundDataL,
                      GroundDataT,
                      GroundDataV
                  ]

"""
"""
function createMap( attributes::AbstractVector, destinations::AbstractVector; stop::Int64=0 )
  min_len = min( length.( [attributes, destinations] )... ) - stop
  return [ attributes[i] => destinations[i] for i in 1:min_len ]
end



"""
"""
function standardize( map::AbstractVector, bridge::Union{ Nothing, Symbol, Pair{Symbol, Symbol} }, dfSta::DataFrame, dfSen::Union{DataFrame, Nothing}=nothing )

  # Separate the column that have a mapping from the others
  complete_map = []
  missing_map = []
  for (x,y) in map
    if isnothing(x)
      push!( missing_map, y )
    else
      push!( complete_map, x => y )
    end
  end

  # Check wether there is a second dataframe or there is only one containing all the needed informations
  if isnothing(dfSen)
      dataframe = select( dfSta, complete_map... )
  else
      dataframe = innerjoin( dfSta, dfSen, on=bridge )
      select!( dataframe, complete_map... )
  end

  # For each column that doesn't have a mapping insert the column as an array of missings
  insertcols!( dataframe, ( missing_map .=> Ref(missings( nrow(dataframe) )) )... )

  return dataframe
end



"""
"""
function generateUuidsTable()
  columns = [ :local_id, :name, :longitude, :latitude ]

  ress = []
  for r in instances(Region)
    for t in [:METEO, :AIRQUALITY]
      rgn = regions_modules[Integer(r)]
      station = rgn.getData( type=t )
      map = createMap( rgn.stat_info[t], columns )
      res = standardize( map, nothing, station )
      push!( ress, res... )
    end
  end

  return ress
end



"""
    getGroundData( filePath::AbstractString="."; regions::Region..., type::Symbol=METEO, source::Symbol=STATIONS )

Obtain `type` data of `regions` from `source` and save it as `filePath` 
"""
function getGroundData( type::Symbol=:METEO, regions::Region... )
  columns = [
              :parameter,                 # tipo di parametro misurato
              :unit,                      # unita di misura (possibilmente SI)
              :value,                     # valore misurato
              :freqency,                  # frequency of measurements
              :date,                      # anno mese giorno ora (UTM)
              :longitude,                 # longitudine della stazione
              :latitude,                  # latitudine della stazione
              :height,                    # quota stazione || quota stazione più quota misura per il vento
              :rel_measurement_height,    # quota relativa della misurazioni
              :validation,                # bool (già segnalato dalla stazione)
              :note                       # errori e outlayers? altro?
            ]

  ress = []
  for region in regions
      println("Region: $region\n")
      rnum = Integer(region)
      rgn = regions_modules[rnum]

      resSta = rgn.getData( type=type, source=:STATIONS )
      resSen = rgn.getData( type=type, source=:SENSORS )

      s = type == :METEO ? 0 : 1 
      map = createMap( rgn.attributes[type], columns, stop=s )
      bridge = rgn.ids[type]
      res = standardize( map, bridge, resSta, resSen )

      push!( ress, res )
  end

  return ress
end

#   res = getGroundData( :METEO, AA, FVG, L, T, V )
#   uuids = generateUuidsTable()


end # module