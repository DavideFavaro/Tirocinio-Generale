module TestCn

#= imports
	import sqlite3
	import ipdb
=#

import DBInterface as dbi
import SQLite as sql


db = sql.DB("substance.db")

query_cn = sql.Stmt( db, "SELECT * FROM cn" )
results = dbi.execute( query_cn )

#listaclc = Dict()
#for row in results
#	lista_soil = [ x for x in row ]
#	listaclc[ row[5] ] = lista_soil	
#end
listaclc = Dict( row[5] => [ x for x in row ] for row in results )

println(listaclc)

end # module