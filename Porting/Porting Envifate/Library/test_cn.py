
import sqlite3
import ipdb
conn = sqlite3.connect("substance.db")
cursor=conn.cursor()

listaclc={}

query_cn="select * from cn"
#import pdb; pdb.set_trace()
cursor.execute(query_cn)

sql_fetch=cursor.fetchall()


for row in sql_fetch:
	lista_soil=[]
	for x in row:
		lista_soil.append(x)
	#ipdb.set_trace()
	listaclc[row[5]]=lista_soil	
	#listaclc.append(lista_soil)

conn.close()

print(listaclc)