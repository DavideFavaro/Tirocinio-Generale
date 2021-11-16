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

using Revise


str = occursin( "GroundData.jl", @__FILE__ ) ? "" : "src\\"
include("$(@__DIR__)\\$(str)GroundDataRegions\\GroundDataAA.jl")
include("$(@__DIR__)\\$(str)GroundDataRegions\\GroundDataER.jl")
include("$(@__DIR__)\\$(str)GroundDataRegions\\GroundDataFVG.jl")
include("$(@__DIR__)\\$(str)GroundDataRegions\\GroundDataL.jl")
include("$(@__DIR__)\\$(str)GroundDataRegions\\GroundDataT.jl")
include("$(@__DIR__)\\$(str)GroundDataRegions\\GroundDataV.jl")


export getGroundData


@enum Region AA=1 FVG=2 L=3 T=4 V=5


const regions_modules = [
                      GroundDataAA,
                      GroundDataFVG,
                      GroundDataL,
                      GroundDataT,
                      GroundDataV
                  ]



Base.convert(::Type{Int64}, s::AbstractString ) = parse( Int64, s )
Base.convert(::Type{Float64}, s::AbstractString ) = isempty(s) || s == "NULL" ? missing : parse( Float64, s )
Base.convert(::Type{DateTime}, s::AbstractString ) = isempty(s) ? missing : DateTime(s)

"""
    pop!( df::DataFrame, columns, condition::Function  )

Delete rows from `df` where `condition` on the elements of `columns` is true and return the deleted rows as a new `DataFrame`
"""
function pop!( df::DataFrame, columns, condition::Function  )
  if columns isa AbstractVector
    res = filter( x -> any( condition, x[columns] ), df )
    filter!( x -> !any( condition, x[columns] ), df )
  else
    res = filter( x -> condition(x[columns]) == true, df )
    filter!( x -> condition(x[columns]) == false, df )
  end
  return res
end



"""
    createMap( attributes::AbstractVector, destinations::AbstractVector[; n::Int64=0] )

Create an array of pairs from the tow input vectors, skipping the last `n` pairs
"""
function createMap( attributes::AbstractVector, destinations::AbstractVector )
  min_len = min( length.( [attributes, destinations] )... )
  return [ attributes[i] => destinations[i] for i in 1:min_len ]
end



"""
    standarize(
    
Generate a dataframe in standard format from `dfSta` and `dfSen` using `bridge` to join the two dataframes and `map` to select the desired columns
and map them to the standard ones.
If `dfSen` is not provided the selection of the column will be done on `dfSta`. 
"""
function standardize( map::AbstractVector, dfSta::DataFrame, dfSen::Union{DataFrame, Nothing}=nothing, bridge::Union{ Nothing, Symbol, Pair{Symbol, Symbol} }=nothing )

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
  insertcols!( dataframe, ( missing_map .=> Ref( missings( nrow(dataframe) ) ) )... )

  return dataframe
end



"""
    generateUuidsTable()

Generate a CSV file containing informations on all the stations of all the regions and their uuids, a second file containing the information on the stations
that lack spatial coordinates and return the dataframe of the first CSV.
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
    for t in [:METEO, :AIRQUALITY]
      rgn = regions_modules[Integer(r)]
      station = try
                  rgn.getData( type=t )
                catch e
                  println( "Skipped $r $t, due to:\n", e )
                  continue
                end
      stat_info = rgn.getRegionStationsInfo(t)
      map = createMap( stat_info, columns )
      res = standardize( map, station )
      insertcols!( res, :type => t, :region => r )
      append!( df, res )
    end
  end

  # AA e V hanno delle stazioni ( 7 e 2 rispettivamente, che non hanno latitudine e longitudine )
  mdf = pop!( df, [:latitude, :longitude], ismissing )
  CSV.write( "./Dati stazioni/missing_stazioni.csv", mdf )
  disallowmissing!(df)
  unique!( df, [:longitude, :latitude] )
  insertcols!( df, :uuid => [ uuid4() for i in 1:nrow(df) ] )
  CSV.write( "./Dati stazioni/stazioni.csv", df )
  return df
end

#   uuids = generateUuidsTable()



"""
    getGroundData( filePath::AbstractString="."; regions::Region..., type::Symbol=METEO, source::Symbol=STATIONS )

Obtain `type` data of `regions` from `source` and save it as `filePath` 
"""
function getGroundData( type::Symbol=:METEO, regions::Region... )

  # `s` will be used to regulate the number of columns of the df based on the data type, the conditional assignment
    # checks also wether `type` is valid
  s = type == :METEO ? 0 :
          type == :AIRQUALITY ? 1 : throw( DomainError( type, "`type` must be either `:METEO` OR `:AIRQUALITY`" ) )

  columns = Pair{Symbol, Vector}[
              :uuid                   => String[],
              :parameter              => String[],                    # tipo di parametro misurato
              :unit                   => Union{String, Missing}[],    # unita di misura (possibilmente SI)
              :value                  => Union{Float64, Missing}[],   # valore misurato
              :frequency              => Union{String, Missing}[],    # frequency of measurements
              :date                   => Union{DateTime, Missing}[],  # anno mese giorno ora (UTM)
              :longitude              => Float64[],                   # longitudine della stazione
              :latitude               => Float64[],                   # latitudine della stazione
              :height                 => Union{Float64, Missing}[],   # quota stazione || quota stazione più quota misura per il vento
              :validation             => Any[],                       # bool (già segnalato dalla stazione)
              :note                   => Any[],                       # errori e outlayers? altro?
              :rel_measurement_height => String[]                     # quota relativa della misurazioni
            ]

  # There are no `:rel_measurement_height` columns for the airquality data so if that is the type of data
    # the created dataframe will lack the specific column
  df = DataFrame( columns[1:end-s] )
  mdf = DataFrame( columns[2:end-s] )
  allowmissing!( mdf, [:longitude, :latitude] )
  for region in regions
      # Use the Region enum number to obtain the corresponding module
      rnum = Integer(region)
      rgn = regions_modules[rnum]
      # Obtain the dataframes that will be used to create the standard format dataframe, while checking for
       # possible problems
      resSta = try 
                  rgn.getData( type=type, source=:STATIONS )
               catch e
                  println( "Skipped $region, due to:\n", e )
                  continue
               end
      resSen = try
                  rgn.getData( type=type, source=:SENSORS )
               catch e
                  println( "Skipped $region, due to:\n", e )
                  continue
               end
      # Create the standard dataframe
      attributes = rgn.getRegionAttributes(type)
      map = createMap( attributes, Symbol.(names(df)[2:end]) )
      bridge = rgn.getRegionIds(type)
      res = standardize( map, resSta, resSen, bridge )
      mres = pop!( res, [:latitude, :longitude], ismissing )
      uuids = isfile(".\\Dati stazioni\\stazioni.csv") ? CSV.read( ".\\Dati stazioni\\stazioni.csv", DataFrame ) : generateUuidsTable()
      res = innerjoin( uuids[:, [:uuid, :longitude, :latitude]], res, on=[:longitude, :latitude] )
      append!( df, res )
      append!( mdf, mres )
  end
  return df, mdf
end

#   resmt, mresmt = getGroundData( :METEO, AA, FVG, L, T, V )
#   resaq, mresaq = getGroundData( :AIRQUALITY, AA, FVG, L, T, V )



"""
"""
function saveGroundData( path::AbstractString, data::DataFrame; overwrite::Bool=false )
  condition = !overwrite && basename(path) in readdir(dirname(path))
  if condition
      # Attributes to confront in order to see if the entry is a duplicate
      attrs = [:uuid, :value, :date, :height, :rel_measurement_height]
      # Preexisting data
      old_data = CSV.read( path, DataFrame )
      filter!( row -> row[attrs] in old_data[:, attrs], data )
  end
  CSV.write( path, data, append=condition )
end

#    path = *( @__DIR__, "\\..\\Dati stazioni")
#    saveGroundData( path*"\\data.csv", resmt )
#    saveGroundData( path*"\\missing_data.csv", mresmt )



end # module