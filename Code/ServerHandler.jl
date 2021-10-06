module ServerHandler

using CSV
using DataFrames
using LibPQ

#include( *( @__DIR__, "\\Code\\data.sql" ) )

params = Dict(
    :dbname => "data",
    :host => "localhost",
    :port => "5432",
    :user => "postgress",
    :password => "password",
    :connect_timeout => "10"
)

con = LibPQ.Connection( "dbname=$(params[:dbname]) host=$(params[:host]) port=$(params[:port]) user=$(params[:user]) password=$(params[:password])" )








"""
    updateProductsVal( key::AbstractString, table::AbstractString, attribute::AbstractString, newval::Any )

For the row in table `table` of the database that has `key` as its pk, update "attribute" to have value `newval`
"""
function updateProductsVal( key::AbstractString, table::AbstractString, attribute::AbstractString, newval )
    # Connettiti al db
    # Effettua l'update dell'attributo nella tabella del db indicata alla riga indicata da key 

    println("Valore modificato")
end

"""
    updateProductsVal( keys::AbstractVector{AbstractString}, table::AbstractString, attribute::AbstractString, newval::Any )

For all the rows in table `table` of the database that have one of `keys` as their pk, update `attribute` to have value `newval`
"""
function updateProductsVal( keys::AbstractVector{AbstractString}, table::AbstractString, attribute::AbstractString, newval )
    # Connettiti al db
    # Effettua l'update dell'attributo nella tabella del db indicata alla riga indicata da key 

    println("Valore modificato")
end



"""
    updateProductsVal( condition::Function, table::AbstractString, attribute::AbstractString, newval::Any )

For all the rows in table `table` of the database that verify `condition`, update `attribute` to have value `newval`
"""
function updateProductsVal( condition::Function, table::AbstractString, attribute::AbstractString, newval )
    # Connettiti al db
    # Effettua l'update dell'attributo nella tabella del db indicata per tutte le righe che rispettano la condizione

    println("Valore modificato")
end







end # module