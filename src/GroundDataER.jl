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


str = occursin( "GroundDataER.jl", @__FILE__ ) ? "" : "src\\"
include("$(@__DIR__)\\$(str)Global.jl")


export getDataER



function getDataER(; type::Data_Type=meteo, source::Data_Source=stations )
    println("Non ancora implementato")
    return nothing
end
end # module