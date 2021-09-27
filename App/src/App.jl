module App





using Genie
using SQLite
using CSV
using DataFrames

Genie.newapp_mvc("CopernicusData")


include(joinpath("config", "initializers", "searchlight.jl"))
SQLite.DB("db/data.sqlite")



Genie.newresource("product")





df = CSV.read( "D:\\Documents and Settings\\DAVIDE-FAVARO\\My Documents\\GitHub\\Tirocinio\\Xml di prova\\data.csv", DataFrame )
cols = names(df)
types = eltype.( eachcol(df) )
# SQLite.createtable!( SQLite.DB("db/data_db.sqlite"), "products", Tables.Schema(cols, types) )




io = open( "D:\\Documents and Settings\\DAVIDE-FAVARO\\My Documents\\GitHub\\Tirocinio\\App\\colonne.txt", "w" )
for i in 1:82
    println( io, "column( :$(cols[i]), $(types[i]) )" )
end
close( io )






end # module