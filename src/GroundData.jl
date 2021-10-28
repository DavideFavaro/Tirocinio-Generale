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
using Dates
using UUIDs


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




Base.convert(::Type{Int64}, s::AbstractString ) = parse( Int64, s )
Base.convert(::Type{Float64}, s::AbstractString ) = isempty(s) ? missing : parse( Float64, s )


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
  if isnothing(dfSen) || isnothing(bridge)
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
  df = DataFrame( 
         :local_id => Union{Missing, Any}[],
         :name => String[],
         :longitude => Union{Missing, Float64}[],
         :latitude => Union{Missing, Float64}[],
         :type => Symbol[],
         :region => Region[]
       )
  for r in instances(Region)
    if r != ER
      for t in [:METEO, :AIRQUALITY]
        rgn = regions_modules[Integer(r)]
        station = rgn.getData( type=t )
        stat_info = rgn.getRegionStationsInfo(t)
        map = createMap( stat_info, columns )
        res = standardize( map, nothing, station )
        insertcols!( res, :type => t, :region => r )
        append!( df, res )
      end
    end
  end

  # AA e V hanno delle stazioni ( 7 e 2 rispettivamente, che non hanno latitudine e longitudine )
  dropmissing!( df, [:longitude, :latitude], disallowmissing=true )

  dict = Dict( (point[:longitude], point[:latitude]) => uuid4() for point in eachrow(unique( df[:, [:longitude, :latitude]] ) ) )
  df[!, :uuid] = getindex.( Ref(dict), (df[!,:longitude], df[!,:latitude])... )

  return df
end

#   uuids = generateUuidsTable()
#   CSV.write( "./Dati stazioni/stazioni.csv", uuids )



"""
    getGroundData( filePath::AbstractString="."; regions::Region..., type::Symbol=METEO, source::Symbol=STATIONS )

Obtain `type` data of `regions` from `source` and save it as `filePath` 
"""
function getGroundData( type::Symbol=:METEO, regions::Region... )
  columns = [
              :parameter,                 # tipo di parametro misurato
              :unit,                      # unita di misura (possibilmente SI)
              :value,                     # valore misurato
              :frequency,                 # frequency of measurements
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
      rnum = Integer(region)
      rgn = regions_modules[rnum]

      resSta = rgn.getData( type=type, source=:STATIONS )
      resSen = rgn.getData( type=type, source=:SENSORS )

      s = type == :METEO ? 0 : 1
      attributes = rgn.getRegionAttributes(type)
      map = createMap( attributes, columns, stop=s )
      bridge = rgn.getRegionIds(type)
      res = standardize( map, bridge, resSta, resSen )

      # res[!, :date] = DateTime.( res[!, :date] )

      push!( ress, res )
  end

  return ress
end

#   res = getGroundData( :METEO, AA, FVG, L, T, V )




end # module