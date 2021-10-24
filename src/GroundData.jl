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


@enum Data_Type     METEO=1 AIRQUALITY=2 
@enum Data_Source   STATIONS=1 SENSORS=2
@enum Region        AA=1 ER=2 FVG=3 L=4 T=5 V=6


functions = [
                getDataAA,
                getDataER,
                getDataFVG,
                getDataL,
                getDataT,
                getDataV
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
bridges = [ # ( METEO, AIRQUALITY )
            ( :SCODE, :SCODE ), # Alto Adige
            ( nothing, nothing ), # Emilia Romagna
            ( nothing, nothing ), # Friuli Venezia Giulia
            ( :idsensore, :idsensore ), # Lombardia
            ( :codice, => :station_id, nothing ), # Trentino
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
maps = [ # Array of couples of arrays containing for each region the mapping of the meteo data (both from stations and sensors) to the standardized form
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
                    nothing => :value,
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
                    :longitude => :longitude,
                    :latitudine => :latitude,
                    :quota => :height,
                    :quota => :rel_measurement_height#=,
                    nothing => :validation ),
                    nothing => :note )=#
                ],
                [ # AIRQUALITY
                    :Inquinante => :parameter,
                    Symbol("Unita di misura") => :unit,
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



"""
"""
function standardize( dfSta::DataFrame, dfSen::DataFrame=nothing, map::AbstractVector{Pair}, bridge::Union{ Pair, Tuple{Pair}} )
    if isnothing(dfSen)
        return select( dfSta, map )
    else
        dataframe = innerjoin( dfs[1], dfs[2], on=bridge )
        select!( dataframe, map )
    return dataframe
end






"""
    getGroundData( filePath::AbstractString="."; regions::Region..., type::Data_Type=METEO, source::Data_Source=STATIONS )

Obtain `type` data of `regions` from `source` and save it as `filePath` 
"""
function getGroundData( targetDirectory::AbstractString="."; regions::Region..., type::Data_Type=METEO )
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













end # module