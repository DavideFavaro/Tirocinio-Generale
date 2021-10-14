module GroundDataEmiliaRomagna

#=
Link ai dati:
    Stazioni qualità aria - Link ottenuto dal bottone per il download (CSV):
        https://docs.google.com/spreadsheets/d/1UFdHVhzULdE-F4o1CFZ7JJ4hC7xJB1IBt_LUCB4K1Uk/export?format=csv

Per gli altri dati non ci sono link espliciti e nemmeno bottoni o opzioni per il download 
=#


@enum Data_Type meteo=1 airquality=2 
@enum Data_Source stations=1 sensors=2

function getDataER( targetDirectory::AbstractString; type::Data_Type=meteo, source::Data_Source=stations )
    println("Non ancora implementato")
    return nothing
end
end # module