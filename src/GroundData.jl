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


include("./src/Global.jl")
include("./GroundDataAA.jl")
include("./GroundDataER.jl")
include("./GroundDataFVG.jl")
include("./GroundDataL.jl")
include("./GroundDataT.jl")
include("./GroundDataV.jl")


export getGroundData


@enum Region AA=1 ER=2 FVG=3 L=4 T=5 V=6


const functions = [
                      GroundDataAA.getDataAA,
                      GroundDataER.getDataER,
                      GroundDataFVG.getDataFVG,
                      GroundDataL.getDataL,
                      GroundDataT.getDataT,
                      GroundDataV.getDataV
                  ]

#= POSSIBILE VARIANTE PER EFFETUARE IL JOIN DEI DF   
bridges = [
            [ # METEO
                :SCODE,                     # Alto Adige
                nothing,                    # Emilia Romagna
                nothing,                    # Friuli Venezia Giulia
                :idsensore,                 # Lombardia
                :codice, => :station_id,    # Trentino
                :idstaz => :station_id      # Veneto
            ]
            [ # AIRQUALITY
                :SCODE,                     # Alto Adige
                nothing,                    # Emilia Romagna
                nothing,                    # Friuli Venezia Giulia
                :idsensore, # Lombardia
                nothing,                    # Trentino
                :codseqst => :station_id    # Veneto
            ]
          ]
=#
const bridges = [ # ( METEO, AIRQUALITY )
                  ( :SCODE, :SCODE ), # Alto Adige
                  ( nothing, nothing ), # Emilia Romagna
                  ( nothing, nothing ), # Friuli Venezia Giulia
                  ( :idsensore, :idsensore ), # Lombardia
                  ( :codice, nothing ), # Trentino
                  ( :idstaz, :codseqst ) # Veneto
                ]


#= POSSIBILE VARIANTE PER IL MAPPING DAL COLONNE DEI VARI DF A QUELLI IN FORMA STANDARD
region_maps = ( # 1: Alto Adige;  2: Emilia Romagna;  3: Friuli Venezia Giulia;  4: Lombardia;  5: Trentino;  6: Veneto;  7: Column Name
                [ # METEO
                    [ :DESC_I, missing, :param, :tipologia, :info, :param,
                      :parameter ],
                    [ :UNIT, missing, :unit, :unit_dimisura, :unit, :unit,
                      :unit ],
                    [ :VALUE, missing, :value, :valore, :value, :value,
                      :value ],
                    [ :DATE, missing, :observation_time, :data, :date, :instant,
                      :date ],
                    [ :LONG, missing, missing, :lng, :longitude, :x,
                      :longitude ],
                    [ :LAT, missing, missing, :lat, :latitudine, :y,
                      :latitude ],
                    [ :ALT, missing, :station_altitude, :quota, :quota, :quota,
                      :height ],
                    [ :ALT, missing, :rel_measure_height, :quota, :quota, :quota,
                      :rel_measurement_height ]#=,
                    [ missing, missing, missing, :stato, missing, missing,
                      :validation ],
                    [ missing, missing, missing, missing, missing, missing,
                     :note ]=#
                ],
                [ # AIRQUALITY
                    [ :MCODE, missing, :parametro, :nometiposensore, :Inquinante, :param,
                      :parameter],
                    [ missing, missing, :unita_misura, :unitamisura, Symbol("Unita di misura"), missing,
                      :unit ],
                    [ :VALUE, missing, missing, :valore, :Valore, :value,
                      :value ],
                    [ :DATE, missing, :data_misura, :data, :Data_Ora, :date,
                      :date ],
                    [ :LONG, missing, :longitudine, :lng, missing, :lon,
                      :longitude ],
                    [ :LAT, missing, :latitudine, :lat, missing, :lat,
                      :latitude ],
                    [ nothing, missing, missing, :quota, missing, missing,
                      :height ]#=,
                    [ :FLAGS, missing, :dati_insuff, :stato, missing, missing,
                      :validation ],
                    [ :NOTE, missing, missing, missing, missing, tipozona,
                      :note ]=#
                ]
              )
=#
const maps = [ # Array of couples of arrays containing for each region the mapping of the meteo data (both from stations and sensors) to the standardized form
                ( # Alto Adige
                    [ # METEO
                        # E' in italiano, non sembra esserci in inglese
                        :DESC_I => :parameter,
                        :UNIT => :unit,
                        :VALUE => :value,
                        :DATE => :date,
                        :LONG => :longitude,
                        :LAT => :latitude,
                        :ALT => :height,
                        :ALT => :rel_measurement_height#=,
                        nothing => :validation ),
                        nothing => :note	),=#
                    ],
                    [ # AIRQUALITY
                        :MCODE => :parameter,
                        nothing => :unit,
                        :VALUE => :value,
                        :DATE => :date,
                        :LONG => :longitude,
                        :LAT => :latitude,
                        nothing => :heigth#=,
                        # Sono tutti missing quindi non so cosa rappresenti
                        :FLAGS => :validation,
                        # :VALUE è -1 quando mancante
                        nothing  => :note=#
                    ]
                ),
                ( # Emilia Romagna
                    [ # METEO
                        nothing => :parameter,
                        nothing => :unit,
                        nothing => :value,
                        nothing => :date,
                        nothing => :longitude,
                        nothing => :latitude,
                        nothing => :height,
                        nothing => :rel_measurement_height#=,
                        nothing => :validation,
                        nothing => :note=#
                    ],
                    [ # AIRQUALITY
                        ( nothing => :parameter ),
                        ( nothing => :unit ),
                        ( nothing => :value ),
                        ( nothing => :date ),
                        ( nothing => :longitude ),
                        ( nothing => :latitude ),
                        ( nothing => :heigth )#=,
                        ( nothing => :validation ),
                        ( nothing => :note )=#
                    ]
                ),
                ( # Friuli Venezia Giulia
                    [ # METEO
                        :param => :parameter,
                        :unit => :unit,
                        :value => :value,
                        :observation_time, :date,
                        nothing => :longitude,
                        nothing => :latitude,
                        :station_altitude => :height,
                        :rel_measure_height => :rel_measurement_height#=,
                        nothing => :validation ),
                        nothing => :note )=#
                    ],
                    [ # AIRQUALITY
                        :parametro => :parameter,
                        :unita_misura => :unit,
                        # Ogni dataframe ha valori diversi ( alcuni hanno media giornaliera altri oraria altri altro )
                        :value => :value,
                        :data_misura => :date,
                        :longitudine => :longitude,
                        :latitudine => :latitude,
                        nothing => :heigth#=,
                        :dati_insuff => :validation,
                        nothing => :note=#
                    ]
                ),
                ( # Lombardia
                    [ # METEO
                        :tipologia => :parameter,
                        :unit_dimisura, :unit,
                        :valore => :value,
                        :data => :date,
                        :lng => :longitude,
                        :lat => :latitude,
                        :quota => :height,
                        :quota => :rel_measurement_height#=,
                        # non so cosa rappresenti
                        :stato => :validation ),
                        :note => nothing )=#
                    ],
                    [ # AIRQUALITY
                        :nometiposensore => :parameter,
                        :unitamisura => :unit,
                        :valore => :value,
                        :data => :date,
                        :lng => :longitude,
                        :lat => :latitude,
                        :quota => :heigth#=,
                        # non so cosa rappresenti | :valore a -9999.0 rappresenta un valore mancante e corrisponde a stao NA apparentemente
                        :stato => :validation,
                        nothing => :note=#
                    ]
                ),
                ( # Trentino
                    [ # METEO
                        # o :attribute ?
                        :info => :parameter,
                        :unit => :unit,
                        :value => :value,
                        :date => :date,
                        :longitudine => :longitude,
                        :latitudine => :latitude,
                        :quota => :height,
                        :quota => :rel_measurement_height#=,
                        nothing => :validation ),
                        nothing => :note )=#
                    ],
                    [ # AIRQUALITY
                        :Inquinante => :parameter,
                        :Unita_di_misura => :unit,
                        :Valore => :value,
                        :Data_Ora => :date,
                        nothing => :longitude,
                        nothing => :latitude,
                        nothing => :heigth#=,
                        nothing => :validation ),
                        nothing => :note )=#
                    ]
                ),
                ( # Veneto
                    [ # METEO
                        :param => :parameter,
                        :unit => :unit,
                        :value => :value,
                        :instant => :date,
                        :x => :longitude,
                        :y => :latitude,
                        :quota => :height,
                        :quota => :rel_measurement_height#=,
                        nothing => :validation,
                        nothing => :note=#
                    ],
                    [ # AIRQUALITY
                        :param => :parameter,
                        nothing => :unit,
                        :value => :value,
                        :date => :date,
                        :lon => :longitude,
                        :lat => :latitude,
                        nothing => :heigth#=,
                        nothing => :validation,
                        :tipozona => :note=#
                    ]
                )
             ]


const columns = [ 
                  :uiid,                      # identificatore univoco delle stazioni di qualunque genere di misura
                  :parameter,                 # tipo di parametro misurato
                  :unit,                      # unita di misura (possibilmente SI)
                  :value,                     # valore misurato
                  :date,                      # anno mese giorno ora (UTM)
                  :longitude,                 # longitudine della stazione
                  :latitude,                  # latitudine della stazione
                  :quote,                     # quota stazione || quota stazione più quota misura per il vento
                  :rel_measurement_height,    # quota relativa della misurazioni
                  :validation,                # bool (già segnalato dalla stazione)
                  :note                       # errori e outlayers? altro?
                ]



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

  println( complete_map )
  println( missing_map )

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
    getGroundData( filePath::AbstractString="."; regions::Region..., type::Data_Type=METEO, source::Data_Source=STATIONS )

Obtain `type` data of `regions` from `source` and save it as `filePath` 
"""
function getGroundData( type::Global.Data_Type=METEO, regions::Region... )
  ress = []
  for region in regions
      rnum = Integer(region)
      fun = functions[ rnum ]

      resSta = fun( type=type, source=STATIONS )
      resSen = fun( type=type, source=SENSORS )

      tnum = integer(type)
      res = standardize( resSta, resSen, maps[rnum][tnum], bridges[rnum][tnum] )

      push!( ress, res )
  end
  return ress
end

#   res = getGroundData( Global.METEO, AA, FVG, L, T, V )











end # module