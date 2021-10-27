module GroundDataER
"""
Module for the download and processing of atmospheric data gathered by measuring stations located in Emilia Romagna, Italy
"""

#=
Link ai dati:
    Stazioni qualità aria - Link ottenuto dal bottone per il download (CSV):
        https://docs.google.com/spreadsheets/d/1UFdHVhzULdE-F4o1CFZ7JJ4hC7xJB1IBt_LUCB4K1Uk/export?format=csv
    Dati AQ:
        https://sdati-test.datamb.it/arex/


Per gli altri dati non ci sono link espliciti e nemmeno bottoni o opzioni per il download 
=#



export getDataER


const attributes = Dict(
                      :METEO      => [ nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing ],
                      :AIRQUALITY => [ nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing ]
                   )

const ids = Dict(
                :METEO      => nothing,
                :AIRQUALITY => nothing
            )

const stat_info = Dict(
                      :METEO      => [ nothing, nothing, nothing, nothing ],
                      :AIRQUALITY => [ nothing, nothing, nothing, nothing ]
                  )



"""
"""
function getDataER(; type::Symbol=:METEO, source::Symbol=:STATIONS )
    println("Non ancora implementato")
    return nothing
end



end # module